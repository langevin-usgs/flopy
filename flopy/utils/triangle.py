import os
import subprocess

import matplotlib.pyplot as plt
import numpy as np

from ..mbase import resolve_exe
from ..utils.cvfdutil import centroid_of_polygon
from ..utils.geospatial_utils import GeoSpatialUtil


class Triangle:
    """
    Class to work with the triangle program to unstructured triangular grids.
    Information on the triangle program can be found at
    https://www.cs.cmu.edu/~quake/triangle.html

    Parameters
    ----------
    model_ws : str
        workspace location for creating triangle files (default is '.')
    exe_name : str
        path and name of the triangle program. (default is triangle, which
        means that the triangle program must be in your path)
    maximum_area : float
        the maximum area for any triangle.  The default value is None, which
        means that the user must specify maximum areas for each region.
    angle : float
        Triangle will continue to add vertices until no angle is less than
        this specified value.  (default is 20 degrees)
    nodes : ndarray
        Two dimensional array of shape (npoints, 2) with x and y positions
        of fixed node locations to include in the resulting triangular mesh.
        (default is None)
    additional_args : list
        list of additional command line switches to pass to triangle

    Returns
    -------
    None

    """

    def __init__(
        self,
        model_ws=".",
        exe_name="triangle",
        maximum_area=None,
        angle=20.0,
        nodes=None,
        additional_args=None,
    ):
        self.model_ws = model_ws
        self.exe_name = resolve_exe(exe_name)
        self.angle = angle
        self.maximum_area = maximum_area
        self._nodes = nodes
        self.additional_args = additional_args
        self._initialize_vars()

    def add_linestring(self, linestring):
        """
        Add a line to the triangle mesh.  The line will be preserved
        in the resulting mesh.

        Parameters
        ----------
        linestring : list, geojson, shapely.geometry, shapefile.Shape
            add linestring method accepts any of these geometries:

            a list of (x, y) points
            geojson LineString object
            shapely LineString object
            shapefile LineString shape
            flopy.utils.geometry.LineString object


        Returns
        -------
        None

        """
        if isinstance(linestring, (list, tuple, np.ndarray)):
            linestring = [linestring]

        geom = GeoSpatialUtil(linestring, shapetype="LineString")
        linestring = geom.points

        self._linestrings.append(linestring)

    def add_polygon(self, polygon, ignore_holes=False):
        """
        Add a polygon

        Parameters
        ----------
        polygon : list, geojson, shapely.geometry, shapefile.Shape
            add polygon method accepts any of these geometries:

            a list of (x, y) points
            geojson Polygon object
            shapely Polygon object
            shapefile Polygon shape
            flopy.utils.geometry.Polygon object
        ignore_holes : bool
            method to ignore holes in polygon and only use the exterior
            coordinates

        Returns
        -------
        None

        """
        if isinstance(polygon, (list, tuple, np.ndarray)):
            polygon = [polygon]

        geom = GeoSpatialUtil(polygon, shapetype="Polygon")
        polygon = geom.points
        if polygon[0][0] == polygon[0][-1]:
            polygon[0] = polygon[0][:-1]
        self._polygons.append(polygon[0])
        if not ignore_holes:
            if len(polygon) > 1:
                for hole in polygon[1:]:
                    self.add_hole(hole)

    def add_hole(self, hole):
        """
        Add a point that will turn enclosing polygon into a hole

        Parameters
        ----------
        hole : tuple
            (x, y)

        Returns
        -------
        None

        """
        self._holes.append(hole)

    def add_region(self, point, attribute=0, maximum_area=None):
        """
        Add a point that will become a region with a maximum area, if
        specified.

        Parameters
        ----------
        point : tuple
            (x, y)

        attribute : integer or float
            integer value assigned to output elements

        maximum_area : float
            maximum area of elements in region

        Returns
        -------
        None

        """
        self._regions.append([point, attribute, maximum_area])

    def build(self, verbose=False):
        """
        Build the triangular mesh

        Parameters
        ----------
        verbose : bool
            If true, print the results of the triangle command to the terminal
            (default is False)

        Returns
        -------
        None

        """

        # provide some protection by removing existing files
        self.clean()

        # write the active domain to a file
        fname = os.path.join(self.model_ws, f"{self.file_prefix}.0.node")
        
        if self._linestrings is not None:
            ls_nds = self._write_nodefile(fname)
        else:
            self._write_nodefile(fname)

        # poly file
        fname = os.path.join(self.model_ws, f"{self.file_prefix}.0.poly")
        if self._linestrings is not None:
            self._write_polyfile(fname,ls_nds)
        else:
            self._write_polyfile(fname)

        # Construct the triangle command
        cmds = [self.exe_name]
        if self.maximum_area is not None:
            cmds.append(f"-a{self.maximum_area}")
        else:
            cmds.append("-a")
        if self.angle is not None:
            cmds.append(f"-q{self.angle}")
        if self.additional_args is not None:
            cmds += self.additional_args
        cmds.append("-A")  # assign attributes
        cmds.append("-p")  # triangulate .poly file
        cmds.append("-V")  # verbose
        cmds.append("-D")  # delaunay triangles for finite volume
        cmds.append("-e")  # edge file
        cmds.append("-n")  # neighbor file
        cmds.append(f"{self.file_prefix}.0")  # output file name

        # run Triangle
        buff = subprocess.check_output(cmds, cwd=self.model_ws)
        buff = buff.decode()
        if verbose:
            print(buff)

        # load the results
        self._load_results()
        self.ncpl = self.ele.shape[0]
        self.nvert = self.node.shape[0]

        # create verts and iverts
        self.verts = self.node[["x", "y"]]
        self.verts = np.array(self.verts.tolist(), float)
        self.iverts = []
        for row in self.ele:
            self.iverts.append([row[1], row[2], row[3]])


    def refine_nds(self, nd_numbers, area, iteration=1, verbose=False):
        # refine existing mesh at specified nodes (refines all elements touching the node)
        c1 = np.isin(self.ele['iv1'], nd_numbers)
        c2 = np.isin(self.ele['iv2'], nd_numbers)
        c3 = np.isin(self.ele['iv3'], nd_numbers)
        ele2refine = np.where(c1 | c2 | c3)[0] # returns the element numbers of all triangles touching the specified nodes
        self.refine_ele(ele2refine, area, iteration, verbose)


    def refine_ele(self, ele_num, area, iteration=1, verbose=False):
        """
        Refine the triangular mesh with specified minimum area at specified elements

        Parameters
        ----------
        ele_num : list or array of element numbers to refine
        area: maximum area constraint to be applied to triangles specified in ele_num
        iteration: iteration number of mesh generation to apply refinement to 
        verbose : bool
            If true, print the results of the triangle command to the terminal
            (default is False)

        Returns
        -------
        None

        """
        fname = os.path.join(self.model_ws, f"{self.file_prefix}.{iteration}.area")
        self._write_areafile(fname, ele_num, area)
        # # provide some protection by removing existing files

        # Construct the triangle command
        cmds = [self.exe_name]
        cmds.append("-r")  # Refine existing mesh
        cmds.append("-a") # apply area constraint 
        if self.maximum_area is not None: # apply area constraint a second time for max area
            cmds.append(f"-a{self.maximum_area}")
        else:
            cmds.append("-a")
        if self.angle is not None:
            cmds.append(f"-q{self.angle}")
        if self.additional_args is not None:
            cmds += self.additional_args
        
        cmds.append("-A")  # assign attributes
        cmds.append("-p")  # triangulate .poly file
        cmds.append("-V")  # verbose
        cmds.append("-D")  # delaunay triangles for finite volume
        cmds.append("-e")  # edge file
        cmds.append("-n")  # neighbor file
        cmds.append(f"{self.file_prefix}.{iteration}")  # output file name

        # run Triangle
        buff = subprocess.check_output(cmds, cwd=self.model_ws)
        buff = buff.decode()
        if verbose:
            print(buff)

        # load the results
        self._load_results(iteration=iteration+1)
        self.ncpl = self.ele.shape[0]
        self.nvert = self.node.shape[0]

        # create verts and iverts
        self.verts = self.node[["x", "y"]]
        self.verts = np.array(self.verts.tolist(), float)
        self.iverts = []
        for row in self.ele:
            self.iverts.append([row[1], row[2], row[3]])


    def plot(
        self,
        ax=None,
        layer=0,
        edgecolor="k",
        facecolor="none",
        cmap="Dark2",
        a=None,
        masked_values=None,
        **kwargs,
    ):
        """
        Plot the grid.  This method will plot the grid using the shapefile
        that was created as part of the build method.

        Note that the layer option is not working yet.

        Parameters
        ----------
        ax : matplotlib.pyplot axis
            The plot axis.  If not provided it, plt.gca() will be used.
            If there is not a current axis then a new one will be created.
        layer : int
            Layer number to plot
        cmap : string
            Name of colormap to use for polygon shading (default is 'Dark2')
        edgecolor : string
            Color name.  (Default is 'scaled' to scale the edge colors.)
        facecolor : string
            Color name.  (Default is 'scaled' to scale the face colors.)
        a : numpy.ndarray
            Array to plot.
        masked_values : iterable of floats, ints
            Values to mask.
        kwargs : dictionary
            Keyword arguments that are passed to
            PatchCollection.set(``**kwargs``).  Some common kwargs would be
            'linewidths', 'linestyles', 'alpha', etc.

        Returns
        -------
        None

        """
        from ..discretization import VertexGrid
        from ..plot import PlotMapView

        cell2d = self.get_cell2d()
        vertices = self.get_vertices()
        ncpl = len(cell2d)

        modelgrid = VertexGrid(
            vertices=vertices, cell2d=cell2d, ncpl=ncpl, nlay=1
        )

        pmv = PlotMapView(modelgrid=modelgrid, ax=ax, layer=layer)
        if a is None:
            pc = pmv.plot_grid(
                facecolor=facecolor, edgecolor=edgecolor, **kwargs
            )
        else:
            pc = pmv.plot_array(
                a,
                masked_values=masked_values,
                cmap=cmap,
                edgecolor=edgecolor,
                **kwargs,
            )

        return pc

    def get_boundary_marker_array(self):
        """
        Get an integer array that has boundary markers

        Returns
        -------
        iedge : ndarray
            integer array of size ncpl containing a boundary ids.  The array
            contains zeros for cells that do not touch a boundary.  The
            boundary ids are the segment numbers for each segment in each
            polygon that is added with the add_polygon method.

        """
        iedge = np.zeros((self.ncpl), dtype=int)
        boundary_markers = np.unique(self.edge["boundary_marker"])
        for ibm in boundary_markers:
            icells = self.get_edge_cells(ibm)
            iedge[icells] = ibm
        return iedge

    def plot_boundary(self, ibm, ax=None, **kwargs):
        """
        Plot a line and vertices for the specified boundary marker

        Parameters
        ----------
        ibm : integer
            plot the boundary for this boundary marker

        ax : matplotlib.pyplot.Axes
           axis to add the plot to.  (default is plt.gca())

        kwargs : dictionary
            dictionary of arguments to pass to ax.plot()

        Returns
        -------
        None

        """
        if ax is None:
            ax = plt.gca()
        idx = np.where(self.edge["boundary_marker"] == ibm)[0]
        for i in idx:
            iv1 = self.edge["endpoint1"][i]
            iv2 = self.edge["endpoint2"][i]
            x1 = self.node["x"][iv1]
            x2 = self.node["x"][iv2]
            y1 = self.node["y"][iv1]
            y2 = self.node["y"][iv2]
            ax.plot([x1, x2], [y1, y2], **kwargs)

    def plot_vertices(self, ax=None, **kwargs):
        """
        Plot the mesh vertices

        Parameters
        ----------
        ax : matplotlib.pyplot.Axes
           axis to add the plot to.  (default is plt.gca())

        kwargs : dictionary
            dictionary of arguments to pass to ax.plot()

        Returns
        -------
        None

        """
        if ax is None:
            ax = plt.gca()
        ax.plot(self.node["x"], self.node["y"], lw=0, **kwargs)

    def label_vertices(self, ax=None, onebased=True, **kwargs):
        """
        Label the mesh vertices with their vertex numbers

        Parameters
        ----------
        ax : matplotlib.pyplot.Axes
           axis to add the plot to.  (default is plt.gca())

        onebased : bool
            Make the labels one-based if True so that they correspond to
            what would be written to MODFLOW.

        kwargs : dictionary
            dictionary of arguments to pass to ax.text()

        Returns
        -------
        None

        """
        if ax is None:
            ax = plt.gca()
        for i in range(self.verts.shape[0]):
            x = self.verts[i, 0]
            y = self.verts[i, 1]
            s = i
            if onebased:
                s += 1
            ax.text(x, y, str(s), **kwargs)

    def plot_centroids(self, ax=None, **kwargs):
        """
        Plot the cell centroids

        Parameters
        ----------
        ax : matplotlib.pyplot.Axes
           axis to add the plot to.  (default is plt.gca())

        kwargs : dictionary
            dictionary of arguments to pass to ax.plot()

        Returns
        -------
        None

        """
        if ax is None:
            ax = plt.gca()
        xcyc = self.get_xcyc()
        ax.plot(xcyc[:, 0], xcyc[:, 1], lw=0, **kwargs)

    def label_cells(self, ax=None, onebased=True, **kwargs):
        """
        Label the cells with their cell numbers

        Parameters
        ----------
        ax : matplotlib.pyplot.Axes
           axis to add the plot to.  (default is plt.gca())

        onebased : bool
            Make the labels one-based if True so that they correspond to
            what would be written to MODFLOW.

        kwargs : dictionary
            dictionary of arguments to pass to ax.text()

        Returns
        -------
        None

        """
        if ax is None:
            ax = plt.gca()
        xcyc = self.get_xcyc()
        for i in range(xcyc.shape[0]):
            x = xcyc[i, 0]
            y = xcyc[i, 1]
            s = i
            if onebased:
                s += 1
            ax.text(x, y, str(s), **kwargs)

    def get_xcyc(self):
        """
        Get a 2-dimensional array of x and y cell center coordinates.

        Returns
        -------
        xcyc : ndarray
            column 0 contains the x coordinates and column 1 contains the
            y coordinates

        """
        ncpl = len(self.iverts)
        xcyc = np.empty((ncpl, 2), dtype=float)
        for i, icell2d in enumerate(self.iverts):
            points = []
            for iv in icell2d:
                x = self.verts[iv, 0]
                y = self.verts[iv, 1]
                points.append((x, y))
            xc, yc = centroid_of_polygon(points)
            xcyc[i, 0] = xc
            xcyc[i, 1] = yc
        return xcyc

    def get_cell2d(self):
        """
        Get a list of the information needed for the MODFLOW DISV Package.

        Returns
        -------
        cell2d : list (of lists)
            innermost list contains cell number, x, y, number of vertices, and
            then the vertex numbers comprising the cell.

        """
        cell2d = []
        xcyc = self.get_xcyc()
        for i, icell2d in enumerate(self.iverts):
            ic2dr = icell2d[::-1]
            cell2d.append([i, xcyc[i, 0], xcyc[i, 1], len(icell2d)] + ic2dr)
        return cell2d

    def get_vertices(self):
        """
        Get a list of vertices in the form needed for the MODFLOW DISV Package.

        Returns
        -------
        vertices : list (of lists)
            innermost list contains vertex number, x, and y

        """
        vertices = []
        for i, row in enumerate(self.verts):
            vertices.append([i, row[0], row[1]])
        return vertices

    def get_edge_cells(self, ibm):
        """
        Get a list of cell numbers that correspond to the specified boundary
        marker.

        Parameters
        ----------
        ibm : integer
            boundary marker value

        Returns
        -------
        cell_list : list
            list of zero-based cell numbers

        """
        # Create the edge dictionary if it doesn't exist
        if self.edgedict is None:
            self._create_edge_dict()

        # Create a list of cells for boundary marker ibm
        cell_list = []
        edgedict = self.edgedict
        for n, ivlist in enumerate(self.iverts):
            itmp = ivlist + [ivlist[0]]
            for i in range(len(ivlist)):
                ie = (itmp[i], itmp[i + 1])
                if ie in edgedict:
                    if edgedict[ie] == ibm:
                        cell_list.append(n)

        return cell_list

    def get_cell_edge_length(self, n, ibm):
        """
        Get the length of the edge for cell n that corresponds to
        boundary marker ibm

        Parameters
        ----------
        n : int
            cell number.  0 <= n < self.ncpl

        ibm : integer
            boundary marker number

        Returns
        -------
        length : float
            Length of the edge along that boundary marker.  Will
            return None if cell n does not touch boundary marker.

        """

        assert 0 <= n < self.ncpl, "Not a valid cell number"

        # Create the edge dictionary if it doesn't exist
        if self.edgedict is None:
            self._create_edge_dict()

        ivlist = self.iverts[n]
        itmp = ivlist + [ivlist[0]]
        d = None
        for i in range(len(ivlist)):
            iv1 = itmp[i]
            iv2 = itmp[i + 1]
            ie = (itmp[i], itmp[i + 1])
            if ie in self.edgedict:
                if self.edgedict[ie] == ibm:
                    x1, y1 = self.verts[iv1]
                    x2, y2 = self.verts[iv2]
                    d = ((x1 - x2) ** 2 + (y1 - y2) ** 2) ** 0.5
                    return d
        return d

    def get_attribute_array(self):
        """
        Return an array containing the attribute value for each cell.  These
        are the attribute values that are passed into the add_region() method.

        Returns
        -------
        attribute_array : ndarray

        """
        return self.ele["attribute"]


    def clean(self):
        """
        Remove the input and output files created by this class and by the
        Triangle program

        Returns
        -------
        None

        """
        # remove input files
        for ext in ["poly", "node"]:
            fname = os.path.join(self.model_ws, f"{self.file_prefix}0.{ext}")
            if os.path.isfile(fname):
                os.remove(fname)
                if os.path.isfile(fname):
                    print(f"Could not remove: {fname}")
        # remove output files
        for ext in ["poly", "ele", "node", "neigh", "edge"]:
            fname = os.path.join(self.model_ws, f"{self.file_prefix}1.{ext}")
            if os.path.isfile(fname):
                os.remove(fname)
                if os.path.isfile(fname):
                    print(f"Could not remove: {fname}")


    def _initialize_vars(self):
        self.file_prefix = "_triangle"
        self.ncpl = 0
        self.nvert = 0
        self._active_domain = None
        self._polygons = []
        self._linestrings = []
        self._holes = []
        self._regions = []
        self.verts = None
        self.iverts = None
        self.edgedict = None

    def _load_results(self,iteration=1):
        # node file
        ext = "node"
        dt = [("ivert", int), ("x", float), ("y", float)]
        fname = os.path.join(self.model_ws, f"{self.file_prefix}.{iteration}.{ext}")
        setattr(self, ext, None)
        with open(fname, "r") as f:
            line = f.readline()
            f.close()
            ll = line.strip().split()
            nvert = int(ll[0])
            ndim = int(ll[1])
            assert ndim == 2, "Dimensions in node file is not 2"
            iattribute = int(ll[2])
            if iattribute == 1:
                dt.append(("attribute", int))
            ibm = int(ll[3])
            if ibm == 1:
                dt.append(("boundary_marker", int))
            a = np.loadtxt(fname, skiprows=1, comments="#", dtype=dt)
            assert a.shape[0] == nvert
            setattr(self, ext, a)

        # ele file
        ext = "ele"
        dt = [("icell", int), ("iv1", int), ("iv2", int), ("iv3", int)]
        fname = os.path.join(self.model_ws, f"{self.file_prefix}.{iteration}.{ext}")
        setattr(self, ext, None)
        with open(fname, "r") as f:
            line = f.readline()
            f.close()
            ll = line.strip().split()
            ncells = int(ll[0])
            npt = int(ll[1])
            assert npt == 3, "Nodes per triangle in ele file is not 3"
            iattribute = int(ll[2])
            if iattribute == 1:
                dt.append(("attribute", int))
            a = np.loadtxt(fname, skiprows=1, comments="#", dtype=dt)
            assert a.shape[0] == ncells
            setattr(self, ext, a)

        # edge file
        ext = "edge"
        dt = [("iedge", int), ("endpoint1", int), ("endpoint2", int)]
        fname = os.path.join(self.model_ws, f"{self.file_prefix}.{iteration}.{ext}")
        setattr(self, ext, None)
        with open(fname, "r") as f:
            line = f.readline()
            f.close()
            ll = line.strip().split()
            nedges = int(ll[0])
            ibm = int(ll[1])
            if ibm == 1:
                dt.append(("boundary_marker", int))
            a = np.loadtxt(fname, skiprows=1, comments="#", dtype=dt)
            assert a.shape[0] == nedges
            setattr(self, ext, a)

        # neighbor file
        ext = "neigh"
        dt = [
            ("icell", int),
            ("neighbor1", int),
            ("neighbor2", int),
            ("neighbor3", int),
        ]
        fname = os.path.join(self.model_ws, f"{self.file_prefix}.{iteration}.{ext}")
        setattr(self, ext, None)
        with open(fname, "r") as f:
            line = f.readline()
            f.close()
            ll = line.strip().split()
            ncells = int(ll[0])
            nnpt = int(ll[1])
            assert nnpt == 3, "Neighbors per triangle in neigh file is not 3"
            a = np.loadtxt(fname, skiprows=1, comments="#", dtype=dt)
            assert a.shape[0] == ncells
            setattr(self, ext, a)

    def _write_nodefile(self, fname):
        # some code to account for possibility of overlaping nodes in input polygons, lines and points
        if self._linestrings is not None:
            verts = [] # list of verticies
            ls_nds =  [] # list to hold node index for each line
            for p in self._polygons:
                for vertex in p:
                    verts.append(vertex)

            nverts=len(verts)
            for j,l in enumerate(self._linestrings):
                ls_nds.append([])
                for v in l:
                    if type(v)==np.ndarray:
                        v=tuple(v)
                    if v not in verts: #if the linestring vertex is not already in the list add it
                        verts.append(v)
                        ls_nds[j].append(nverts)
                        nverts=nverts+1
                    else: #vertex is already in the list
                        # get the vertex index
                        vert_ind=np.where((np.array(verts)==np.array(v)).all(axis=1))[0][0]
                        # add the index of that node to the list for that line
                        ls_nds[j].append(vert_ind)
                        
            if self._nodes is not None:
                for i in range(self._nodes.shape[0]):
                    v=tuple(self._nodes[i])
                    if v not in verts:
                        verts.append(tuple(self._nodes[i]))
                        nverts=nverts+1
            
        f = open(fname, "w")
        nvert = 0
        for p in self._polygons:
            nvert += len(p)
        if self._nodes is not None:
            nvert += self._nodes.shape[0]
        if self._linestrings is not None:
            nvert=len(verts)
        s = f"{nvert} 2 0 0\n"
        f.write(s)

        if self._linestrings is not None:
            for i,v in enumerate(verts):
                    s = f"{i} {v[0]} {v[1]}\n"  # {'  '} {j+1}\n"
                    f.write(s)
            f.close()           
        else:
            ip = 0
            for p in self._polygons:
                for vertex in p:
                    s = f"{ip} {vertex[0]} {vertex[1]}\n"
                    f.write(s)
                    ip += 1
      
            if self._nodes is not None:
                for i in range(self._nodes.shape[0]):
                    s = f"{ip} {self._nodes[i, 0]} {self._nodes[i, 1]}\n"
                    f.write(s)
                    ip += 1
    
            f.close()
        if self._linestrings is not None:
            return ls_nds
        

    def _write_areafile(self, fname, ele_num, area):
        if type(area) in (int, float, np.float64):
            area = np.ones(len(ele_num)) * area
        with open(fname, "w") as f:
            nele = self.ele.shape[0]
            s = f"{nele}\n"
            f.write(s)
            j = 0
            for e in self.ele:
                if e[0] in ele_num:
                    s = f"{e[0]} {area[j]}\n"
                    f.write(s)
                    j += 1
                else:
                    s = f"{e[0]} {-1}\n"
                    f.write(s)                
        return
    
        
    def _write_polyfile(self, fname, ls_nds=None):
        f = open(fname, "w")

        # vertices, write zero to indicate read from node file
        s = "0 0 0 0\n"
        f.write(s)

        # segments are sum of polygons and linestrings
        nseg = 0
        for p in self._polygons:
            nseg += len(p)
        for l in self._linestrings:
            # linestrings are a coordinate list so # of segments is # of coordinates -1
            nseg += len(l) - 1 
        
        bm = 1
        s = f"{nseg} {bm}\n"
        f.write(s)

        iseg = 0
        ipstart = 0
        for p in self._polygons:
            nseg = len(p)
            # number of segments is equal to the number of 
            # points.  This way, there is a final segment that
            # closes the polygon
            for i in range(nseg):
                ep1 = i
                ep2 = i + 1
                if ep2 > nseg - 1:
                    ep2 = 0
                ep1 += ipstart
                ep2 += ipstart
                s = f"{iseg} {ep1} {ep2} {iseg + 1}\n"
                f.write(s)
                iseg += 1
            ipstart += len(p)

        # linestrings
        if self._linestrings is not None:
            poly_iseg=iseg
            for j,l in enumerate(self._linestrings):
                nseg=len(l)-1
                for i in range(nseg):
                    ep1 = ls_nds[j][i]
                    ep2 = ls_nds[j][i + 1]
                    s = f"{iseg} {ep1} {ep2} {'  '} {poly_iseg + j + 1}\n"
                    f.write(s)
                    iseg += 1
                ipstart += len(p)
                    
        
        # holes
        nholes = len(self._holes)
        s = f"{nholes}\n"
        f.write(s)
        for i, hole in enumerate(self._holes):
            s = f"{i} {hole[0]} {hole[1]}\n"
            f.write(s)

        # regions
        nregions = len(self._regions)
        s = f"{nregions}\n"
        f.write(s)
        for i, region in enumerate(self._regions):
            pt = region[0]
            attribute = region[1]
            maxarea = region[2]
            if maxarea is None:
                maxarea = -1.0
            s = f"{i} {pt[0]} {pt[1]} {attribute} {maxarea}\n"
            f.write(s)

        f.close()

    def _create_edge_dict(self):
        """
        Create the edge dictionary

        """
        edgedict = {}
        for _, iv1, iv2, iseg in self.edge:
            if iseg != 0:
                edgedict[(iv1, iv2)] = iseg
                edgedict[(iv2, iv1)] = iseg
        self.edgedict = edgedict

    def unique_vertices(xy_vertarray_list):
        """
        This routine looks for duplicate points in one or more
        x,y vertex arrays provided as input.  The routine returns
        a new vertex array with only the unique points and a separate
        index array equal in size to the sum of all the rows provided
        in xy_vertarray_list.  The index array contains the index number
        of the original x,y row in the new unique array.
        """
        # stack the vertex arrays
        combined_array = np.vstack(xy_vertarray_list)

        # find unique x,y pairs in the combined list
        unq = np.unique(combined_array, axis=0)

        # find the index of each combined_array row in the unique array
        nrows = combined_array.shape[0]
        index_in_unq = np.empty(nrows, dtype=int)
        for i, row in enumerate(combined_array):
            indices = np.argwhere(np.all(combined_array == row, axis=1))
            index_in_unq[i] = indices.ravel().min()

        return unq, index_in_unq