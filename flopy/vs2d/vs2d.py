import os

import numpy as np

from ..mbase import BaseModel


class Vs2d(BaseModel):
    def __init__(
        self,
        modelname="vs2dt",
        namefile_ext="fil",
        exe_name="vs2d.exe",
        model_ws=".",
        verbose=False,
        load=True,
        silent=0,
    ):
        # Call constructor for parent object
        BaseModel.__init__(self, modelname, namefile_ext, exe_name, model_ws)
        self.verbose = verbose
        self.load = load
        self.silent = silent
        return

    def write_name_file(self):
        """
        Write vs2dt.fil

        vs2dt.dat        file05.dat
        vs2dt.out        file06.dat
        file07.out       file07.dat
        variables.out    file08.dat
        balance.out      file09.dat
        obsPoints.out    file11.dat
        # vs2dt3.3       file0x.dat

        """
        fn_path = os.path.join(self.model_ws, self.namefile)
        f_nam = open(fn_path, "w")
        files = [
            "vs2dt.dat",
            "vs2dt.out",
            "file07.out",
            "variables.out",
            "balance.out",
            "obspoints.out",
            "# vs2dt3.3",
        ]
        for f in files:
            f_nam.write("{}\n".format(f))
        f_nam.close()
        return

    def load_balance(self):
        """

        Returns
        -------
        balance : ndarray
            structured array of balance terms versus time.

        """

        nly = self.vs2dt.nly
        nxr = self.vs2dt.nxr

        fname = os.path.join(self.model_ws, "balance.out")
        column_width = 21
        if self.vs2dt.numt > 0:
            raise Exception(
                "Cannot load balance.out unless extra precision"
                " is used (numt is negative)"
            )
        if os.path.isfile(fname):
            with open(fname, "r") as f:
                line = f.readline()
                l1 = [
                    line[i : i + column_width]
                    for i in range(0, len(line), column_width)
                ]
                line = f.readline()
                l2 = [
                    line[i : i + column_width]
                    for i in range(0, len(line), column_width)
                ]
                line = f.readline()
                l3 = [
                    line[i : i + column_width]
                    for i in range(0, len(line), column_width)
                ]
                column_names = []
                for i in range(len(l1)):
                    name = l1[i].strip()
                    l2string = l2[i].strip()
                    l3string = l3[i].strip()
                    if len(l2string) > 0:
                        if name.endswith("-"):
                            name = name[:-1] + l2string
                        else:
                            name += "_" + l2string
                    if len(l3string) > 0:
                        name += "_" + l3string
                    name = name.replace(" ", "_")
                    name = name.replace("-", "")
                    column_names.append(name)
                dt = [(name, float) for name in column_names]
                balance = np.loadtxt(f, dt)

        return balance

    def load_obs(self):
        """

        Returns
        -------
        results : tuple
            obs_results is a numpy structured array

        """

        obs_results = None
        nly = self.vs2dt.nly
        nxr = self.vs2dt.nxr

        fname = os.path.join(self.model_ws, "obsPoints.out")
        if os.path.isfile(fname) and os.stat(fname).st_size > 0:
            obs_dtype = [
                ("time", float),
                ("node", int),
                ("x", float),
                ("z", float),
                ("head", float),
                ("pressure_head", float),
                ("theta", float),
                ("saturation", float),
                ("vx", float),
                ("vz", float),
                ("et", float),
            ]
            obs_results = np.loadtxt(fname, dtype=obs_dtype, skiprows=3)

        else:
            print("Could not find obs file {}".format(fname))

        return obs_results

    def load_variables(self):
        """

        Returns
        -------
        results : ndarray
            is an array of the dependent variable (head) for each time they
            were saved.

        """

        variables = []
        times = []
        nly = self.vs2dt.nly
        nxr = self.vs2dt.nxr

        fname = os.path.join(self.model_ws, "variables.out")
        if os.path.isfile(fname):
            with open(fname, "r") as f:
                line = f.readline()  # read blank line
                while True:
                    line = f.readline()  # read the time line
                    if not line:
                        break
                    linelist = line.strip().split()
                    time_in_file = float(linelist[2])
                    times.append(time_in_file)
                    # assert np.allclose(np.float(time_in_file), np.float(tm))
                    line = f.readline()  # read blank line
                    if not line:
                        break
                    a = np.fromfile(
                        f, dtype=np.float, count=nly * nxr, sep=" "
                    )
                    a = a.reshape((nly, nxr))
                    variables.append(a)
        else:
            print("Could not find variables.out {}".format(fname))

        return variables, times

    def load_binary_results(self, fname=None):
        from ..utils.binaryfile import binaryread

        if fname is None:
            fname = os.path.join(self.model_ws, "fort.12")
        phead_list = []
        times = []
        if os.path.isfile(fname) and os.stat(fname).st_size > 0:
            nxr = self.vs2dt.nxr
            nly = self.vs2dt.nly
            with open(fname) as f:
                while True:
                    try:
                        stim = binaryread(f, np.float64)[0]
                    except:
                        break
                    if not stim:
                        break
                    times.append(stim)
                    phead = binaryread(f, np.float64, shape=(nxr, nly))
                    phead = np.transpose(phead)
                    phead_list.append(phead)
        return phead_list, times

    def plot_obs(self, axes=None, obs=None, figsize=(10, 10), **kwargs):
        """

        Parameters
        ----------
        axes : list of matplotlib.axes.Axes
            List of axes onto which the observation plots will be made.  There
            must be enough axes to make the plots.  A separate plot will be
            made for each dependent variable in the obsPoints.out file.  For
            a vs2dt flow simulation, there will be 7 plots (head, pressure
            head, saturation, theta, vx, vz, et).  The default is None,
            which means the axes will be automatically generated.
        obs : list of (layer, column) observation tuples
            This is typically the list of observations locations that is
            passed to the vs2dt package constructor.  The default is None,
            which means it will use the vs2dt.obs.
        figsize : tuple
            Size of the figure.  Default is (10, 10)
        kwargs : dictionary
            Keyword arguments that are passed into ax.plot() method that is
            used to create the plots.

        Returns
        -------
        axes : list of matplotlib.axes.Axes

        """
        import matplotlib.pyplot as plt

        # load the results and determine number of plots
        obs_results = self.load_obs()
        if obs_results is None:
            print("No observations to plot")
            return
        if obs is None:
            obs = self.vs2dt.obs
        else:
            assert isinstance(obs, list), "obs must be a list of (lay, col)"
        colnames = obs_results.dtype.names[4:]
        nplots = len(colnames)

        # if list of axes not passed in, then create them
        if axes is None:
            fig, axes = plt.subplots(nplots, 1, sharex=True, figsize=figsize)
        assert nplots <= len(
            axes
        ), "axes must be have at least {} axes".format(nplots)

        # plotting keywords
        marker = "o"
        if "marker" in kwargs:
            marker = kwargs.pop("marker")

        for layer, row in obs:
            node = (
                row * self.vs2dt.nly + layer
            )  # node numbering NOT like MODFLOW
            node += 1
            obsn = obs_results[obs_results["node"] == node]

            for i, ax in enumerate(axes):
                colname = obs_results.dtype.names[i + 4]
                ax.plot(
                    obsn["time"],
                    obsn[colname],
                    marker=marker,
                    label="CELL ({},{})".format(layer + 1, row + 1),
                    **kwargs,
                )

        # add ylabel and legend and label x axis for last plot
        for i, ax in enumerate(axes):
            colname = obs_results.dtype.names[i + 4]
            ax.set_ylabel(colname.upper())
            ax.legend()
        axes[-1].set_xlabel("TIME ({})".format(self.vs2dt.tunit))

        return axes

    def plot_variables(
        self,
        axes=None,
        figsize=(10, 10),
        variables=None,
        plot_titles=None,
        **kwargs,
    ):
        import matplotlib.pyplot as plt

        from ..utils import SpatialReference

        # load the results and determine number of plots
        if variables is None:
            variables, times = self.load_variables()
            if plot_titles is None:
                plot_titles = []
                for t in times:
                    plot_titles.append("PRESSURE HEAD (time={})".format(t))
        nplots = len(variables)

        # if list of axes not passed in, then create them
        if axes is None:
            fig, axes = plt.subplots(nplots, 1, figsize=figsize)
            if nplots == 1:
                axes = [axes]
        assert nplots <= len(
            axes
        ), "axes must be have at least {} axes".format(nplots)

        delr = self.vs2dt.get_delr()
        delc = self.vs2dt.get_dell()
        sr = SpatialReference(delr, delc)
        for iplot, a in enumerate(variables):
            atoplot = variables[iplot]
            atoplot[0, :] = 1.0e30
            atoplot[-1, :] = 1.0e30
            atoplot[:, 0] = 1.0e30
            atoplot[:, -1] = 1.0e30
            atoplot = np.ma.masked_equal(atoplot, 1.0e30)
            qm = sr.plot_array(atoplot, ax=axes[iplot], **kwargs)
            plt.colorbar(qm, shrink=0.5, ax=axes[iplot])

        if plot_titles is not None:
            for iplot, plot_title in enumerate(plot_titles):
                axes[iplot].set_title(plot_title)
        plt.tight_layout()

        return axes

    # todo def plot_balance(self):
