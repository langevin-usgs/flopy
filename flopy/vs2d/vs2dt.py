"""
vs2dt module. Contains the vs2dt package class. Note that the user can access
the Vs2dt class as `flopy.vs2d.Vs2dt`.

Additional information for this Vs2dt package can be found in the Vs2dt
User's Manual.

"""

import numpy as np

from ..pakbase import Package
from ..utils.flopy_io import line_parse
from ..utils.util_array import read1d


def load_record(f, recname, varnames, vartypes, vardict, verbose):
    if verbose:
        print("   loading {}...".format(recname))

    line = f.readline()
    linelist = line_parse(line)

    for i, (varname, vartype) in enumerate(zip(varnames, vartypes)):
        if verbose:
            print("{}:".format(recname))

        if vartype is bool:
            val = linelist[i] == "T"
        else:
            val = vartype(linelist[i])

        # if name ends in tr, then it is a transient variable
        if varname.endswith("_tr"):
            if varname in vardict:
                if isinstance(vardict[varname], list):
                    vardict[varname].append(val)
                else:
                    val0 = vardict[varname]
                    vardict[varname] = [val0, val]
            else:
                vardict[varname] = val
        else:
            vardict[varname] = val

        if verbose:
            print("    {} = {}".format(varname, val))

    return vardict


def write_record(f, vars, recname, irech=None, verbose=False):
    if verbose:
        print("... Writing record {}".format(recname))
    if irech is None:
        varlist = [str(var) for var in vars]
    else:
        varlist = []
        for var in vars:
            if isinstance(var, list):
                v = var[irech]
                if isinstance(v, bool):
                    v = bool2str(v)
                varlist.append(str(v))
            else:
                v = var
                if isinstance(v, bool):
                    v = bool2str(v)
                varlist.append(str(v))
    s = " ".join(varlist)
    if recname is not None:
        s += "  Record {}".format(recname)
    if irech is not None:
        s += "  Period {}".format(irech + 1)
    s += "\n"
    f.write(s)
    return


def bool2str(b):
    if b:
        return "T"
    else:
        return "F"


class Vs2dt(Package):
    """
    vs2dt package

    Parameters
    ----------
    model : model object
        The model object (of type :class:`flopy.vs2d.Vs2d`) to which
        this package will be added.
    titl : str
        80-character problem description.
    tmax : float
        Maximum simulation time.  Default is 1.0.
    stim : float
        Initial time.  Default is 0.0.
    ang : float
        Angle by which grid is to be tilted (Must be between -90 and +90
        degrees, ANG = 0 for no tilting, see Healy (1990) for further
        discussion), degrees. Default is 0.
    zunit : str
        Units used for length (A4), 'm' for meters. Default is m.
    tunit : str
        Units used for time (A4), 'day' for days. Default is day.
    cunx : str
        Units used for mass (A4), 'kg' for kilograms. Default is kg.
    nxr : int
        Number of cells in horizontal or radial direction.  Default is 1.
    nly : int
        Number of cells in vertical direction.  Default is 10.
    nrech : int
        Number of recharge periods. (NOTE: set NRECH to a negative number
        (-1 times actual number of recharge periods) to output binary values
        of head and concentration at selected observation times to file
        fort.12. Selecting this option allows the simulation to be restarted
        at any observation time; however, it may require a large amount of
        disk storage space.)  Default is -1.
    numt : int
        Maximum number of time steps. (NOTE: if enhanced precision in print
        out to file 9 and file 11 is desired set NUMT equal to a negative
        number. That is, multiply actual maximum number of time steps by –1)
        Default is -100.
    rad : bool
        Logical variable = T if radial coordinates are used; otherwise = F.
        Default is False.
    itstop : bool
        Logical variable = T if simulation is to terminate after ITMAX
        iterations in one time step; otherwise = F.
        Default is True.
    trans : bool
        Logical variable = T if solute transport is to be simulated.
        Default is False.
    cis : bool
        Logical variable = T if centered-in-space differencing is to be
        used; = F if backward-in-space differencing is to be used for
        transport equation. Default is True.
    cit : bool
        Logical variable = T if centered-in-time differencing is to be
        used; = F if backward-in-time or fully implicit differencing is to
        be used. Default is False.
    sorp : bool
        Unknown required variable read by program but not documented.
        Default is False.
    f6p : bool
        Logical variable = T if mass balance is to be written to file 6 for
        each time step; = F if mass balance is to be written to file 6 only
        at observation times and ends of recharge periods. Default is True.
    thpt : bool
        Logical variable = T if volumetric moisture contents are to be written
        to file 6; otherwise = F.  Default is True.
    spnt : bool
        Logical variable = T if saturations are to be written to file 6;
        otherwise = F.  Default is True.
    ppnt : bool
        Logical variable = T if pressure heads are to be written to file 6;
        otherwise F.  Default is True.
    hpnt : bool
        Logical variable = T if total heads are to be written to file 6;
        otherwise = F.  Default is True.
    vpnt : bool
        Logical variable = T if velocities are to be written to file 6;
        otherwise = F. Default is True.
    ifac : int
        = 0 if grid spacing in horizontal (or radial) direction is to be read
        in for each column and multiplied by FACX.
        = 1 if all horizontal grid spacing is to be constant and equal to FACX.
        = 2 if horizontal grid spacing is variable, with spacing for the first
        two columns equal to FACX and the spacing for each subsequent column
        equal to XMULT times the spacing of the previous column, until the
        spacing equals XMAX, where- upon spacing becomes constant at XMAX.
        Default is 1.
    facx : float
        Constant grid spacing in horizontal (or radial) direction (if IFAC = 1);
        constant multiplier for all spacing (if IFAC = 0); or initial spacing
        (if IFAC = 2), L. Default is 1.0.
    dxr : float
        Grid spacing in horizontal or radial direc tion. Number of entries
        must equal NXR, L.  Default is None.
    xmult : float
        Multiplier by which the width of each cell is increased from that of
        the previous cell. Default is None.
    xmax : float
        Maximum allowed horizontal or radial spacing, L. Default is None.
    jfac : int
        = 0 if grid spacing in vertical direction is to be read in for each
        row and multiplied by FACZ.
        = 1 if all vertical grid spacing is to be constant and equal to FACZ.
        = 2 if vertical grid spacing is variable, with spacing for the first
        two rows equal to FACZ and the spacing for each subsequent row equal
        to ZMULT times the spacing at the previous row, until spacing equals
        ZMAX, whereupon spacing becomes constant at ZMAX.  Default is 1.
    facz : float
        Constant grid spacing in vertical direction (if JFAC=1); constant
        multiplier for all spacing (if JFAC=0); or initial vertical spacing
        (if JFAC=2), L. Default is 1.0
    delz : float
        Grid spacing in vertical direction; number of entries must equal
        NLY, L. Default is None.
    zmult : float
        Multiplier by which each cell is increased from that of previous cell.
        Default is None.
    zmax : float
        Maximum allowed vertical spacing, L. Default is None.
    pltim : list
        Elapsed times at which pressure heads and concentrations are written
        to file 8, and heads, concentrations, saturations, velocties, and/or
        moisture contents to file 6, T. Default is None.
    obs : list
        List of zero-based (layer, column) observations. nobs is calculated
        as the length of this list.  Default is None.
    mb9 : list
        The index number of each mass balance component to be written to
        file 9 (balance.out).  (See table 7, from p. 66, in Healy (1990) listed
        at end of these instructions.). nmb9 is calculated as the length of
        this list. Default is a list from 1 to 72 with every third entry.
        This restricts the output to just the simulation cumulative totals
        and does not include the time step totals and time step rates.  Thus,
        the default is range(1, 73, 3).
    eps : float
        Head closure criterion for iterative solution of flow equation, L.
        Default is 0.01.
    hmax : float
        Relaxation parameter for iterative solution. See discussion in Lappala
        and others (1987) for more detail. Value is generally in the range of
        0.4 to 1.2. Default is 0.7.
    wus : float
        Weighting option for intercell relative hydraulic conductivity:
        WUS = 1 for full upstream weighting. WUS = 0.5 for arithmetic mean.
        WUS = 0.0 for geometric mean. Default is 1.0.
    eps1 : float
        Concentration closure criterion for iterative solution of transport
        equation, M/L3. Default is 0.01.
    minit : int
        Minimum number of iterations per time step. Default is 1.
    itmax : int
        Maximum number of iterations per time step. Default is 100.
    phrd : bool
        Logical variable = T if initial conditions are read in as
        pressure heads; = F if initial conditions are read in as moisture
        contents. Default is True.
    ntex : int
        Number of textural classes or lithologies having different values of
        hydraulic conductivity, specific storage, and/or constants in the
        functional relations among pressure head, relative conductivity, and
        moisture content, must be less than 11. Default is 1.
    nprop : int
        Number of flow properties to be read in for each textural class. When
        using Brooks and Corey, van Genuchten or Nimmo- Rossi functions, set
        NPROP = 6; when using Haverkamp functions, set NPROP =8. When using
        tabulated data, set NPROP = 6 plus number of data points in table.
        [For example, if the number of pressure Heads in the table is equal
        to N1, then set NPROP=3*(N1+1)+3]. Default is 6.
    nprop1 : int
        Number of transport properties to be read in for each textural class.
        For no adsorption set NPROP1 = 6. For Langmuir or Freundlich isotherm
        set NPROP1 = 7. For ion exchange set NPROP1 = 8. Present only if
        TRANS = T. Default is 6.
    hft : int
        Hydraulic function type, 0 for Brooks- Corey; 1 for van Genuchten;
        2 for Haverkamp; 3 for tabular data; and 4 for Rossi-Nimmo.
        Default is 1.
    adt : int
        Adsorption type, 1 for linear adsorption; 2 Langmuir isotherm;
        3 for Freundlich isotherm; 4 for mono-monovalent ion ex- change;
        5 for mono-divalent ion exchange; 6 for di-monovalent ion exchange;
        and 7 di-divalent exchange. Default is 1.
    aniz : tuple
        Ratio of hydraulic conductivity in the z-coordinate direction to that
        in the x-coordinate direction for textural class ITEX.
        Default is (1.0,).
    hk : tuple
        hk(itex, 1) is saturated hydraulic conductivity (K) in the x-coordinate
        direction for class itex, L/T.
        hk(itex 2) is specific storage (Ss) for class itex, L-1.
        hk(itex, 3) is porosity (f) for class itex.
        Definitions for the remaining sequential values on this line are
        dependent upon which functional relation is selected to represent the
        nonlinear coefficients. Five different functional relations are
        allowed: (0) Brooks and Corey, (1) van Genuchten, (2) Haverkamp,
        (3) tabular data, and (4) Rossi-Nimmo. In the following descriptions,
        definitions for the different functional relations are indexed by the
        above numbers. For tabular data, all pressure heads are input first
        (in decreasing order from the largest to the smallest), all relative
        hydraulic conductivities are then input in the same order, followed by
        all for additional details.
        hk(itex, 4) is
        (0) hb, Brooks-Corey bubbling pressure head (must be less than 0), L;
        (1) alpha, van Genuchten alpha. NOTE: alpha is as defined by van
        Genuchten (1980) and is the negative reciprocal of alpha' used in
        earlier versions (prior to version 3.0) of VS2DT, L;
        (2) A', Haverkamp parameter (must be less than 0.0), L;
        (3) Largest pressure head in table.
        (4) psi0, Rossi-Nimmo parameter.
        hk(itex, 5) is
        (0) Residual moisture content;
        (1) Residual moisture content;
        (2) Residual moisture content;
        (3) Second largest pressure head in table;
        (4) psiD, Rossi-Nimmo parameter.
        hk(itex, 6) is
        (0) Brooks-Corey pore-size distribution index;
        (1) n, van Genuchten parameter, beta' in Healy (1990) and Lappala and
        others (1987);
        (2) Beta', Haverkamp parameter;
        (3) Third largest pressure head in table;
        (4) lambda, Rossi-Nimmo parameter.
        hk(itex, 7) is
        (0) not used;
        (1) note used;
        (2) alpha, Haverkamp parameter (must be less than 0.0), L;
        (3) Fourth largest pressure head in table;
        (4) not used.
        hk(itex, 8) is
        (0) not used;
        (1) not used;
        (2) beta, Haverkamp parameter;
        (3) Fifth largest pressure head in table;
        (4) not used.
        For functional relations (0), (1), (2), and (4) no further values are
        required on this line for this textural class. For tabular data (3),
        data input continues as follows:
        hk(itex, 9) is next largest pressure head in table.
        hk(itex, n1 + 3) is minimum pressure head in table.  Here n1 = number
        of pressure heads in table; nprop = 3 * (n1 + 1) + 3)
        hk(itex, n1 + 4) is always input a value of 99.
        hk(itex, n1 + 5) is relative hydraulic conductivity corresponding to
        first pressure head.
        hk(itex, n1 + 6) is relative hydraulic conductivity corresponding to
        second pressure head.
        hk(itex, 2 * n1 + 4) is relative hydraulic conductivity corresponding
        smallest pressure head.
        hk(itex, 2 * n1 + 5) is always input a value of 99.
        hk(itex, 2 * n1 + 6) is moisture content corresponding to first
        pressure head.
        hk(itex, 2 * n1 + 7) is moisture content corresponding to second
        pressure head.
        hk(itex, 3 * n1 + 5) is moisture content corresponding to smallest
        pressure head.
        hk(itex, 3 * n1 + 6) is always input a value of 99.
        Regardless of which functional relation is selected there must be
        nprop + 1 values on line B-7.
        Default is ((1.0, 1.e-5, 0.2, 0.0, -0.8, 2.0), )
    ht : tuple
        ht(itex, 1) is longitudinal dispersivity.
        ht(itex, 2) is transverse dispersivity.
        ht(itex, 3) is molecular diffusion coefficient.
        ht(itex, 4) is decay constant.
        ht(itex, 5) is bulk density (set to zero for no adsoprtion or ion
        exchange).
        ht(itex, 6) is 0 for no adsorption or ion exchange; is Kd for linear
        adsorption isotherm; is Kl for Langmuir isotherm; is Kf for Freundlich
        isotherm; or is Km for ion exchange.
        ht(itex, 7) is Q for Langmuir isotherm; n for Freundlich isotherm; or
        Q for ion exchange. Not used when adsorption or exchange is not
        simulated.
        ht(itex, 8) is C0 for ion exchange; only used for ion exchange.
        Default is None.
    jtex : int or ndarray of shape (nly, nxr)
        Indices (itex) for textural class for each node, read in row by row.
        There must be nly * nxr entries. Default is 1.
    iread : int
        If IREAD = 0, all initial conditions in terms of pressure head or
        moisture content as determined by the value of PHRD are set equal to
        FACTOR. If IREAD = 1, all initial conditions are read from file IU in
        user-designated format and multiplied by FACTOR. If IREAD = 2 initial
        conditions are defined in terms of pressure head, and an equilibrium
        profile is specified above a free-water surface at a depth of DWTX
        until a pressure head of HMIN is reached. all pressure heads above
        this are set to HMIN. If IREAD = 3 initial heads and concentrations are
        read unformatted from file fort.13 for continuation of a previous
        simulation beginning at time STIM (line A-2). Default is 0.
    factor : float
        Multiplier or constant value, depending on value of IREAD, for initial
        conditions. Default is 0.0.
    dwtx : float
        Depth to free-water surface above which an equilibrium profile is
        computed, L. Default is None.
    hmin : float
        Minimum pressure head to limit height of equilibrium profile, L.
        Must be negative. Default is None.
    iu : int
        Unit number from which initial head or moisture content values are to
        be read. Default is None.
    ifmt : str
        Format to be used in reading initial values from unit IU. Must be
        enclosed in quotation marks, for example '(10X,E10.3)'.
        Default is 'free'.
    pinit : float or ndarray of size (nly, nxr)
        Initial conditions.  If phrd is True, then pinit is initial pressure
        head.  Otherwise pinit is initial moisture content. Default is None.
    bcit : bool
        Logical variable = T if evaporation is to be simulated at any time
        during the simulation; otherwise = F. Default is False.
    etsim : bool
        Logical variable = T if evapotranspiration (plant-root extraction) is
        to be simulated at any time during the simulation. Default is False.
    npv : int
        Number of ET periods to be simulated. NPV values for each variable
        required for the evaporation and/or evapotranspiration options must be
        entered on the following lines. If ET variables are held constant
        throughout the simulation code, NPV = 1. (NOTE: Set NPV equal to a
        negative number {-1 times number of ET periods} if solute uptake by
        plant roots is not allowed; otherwise, solute is removed from the domain
        by root uptake.). Default is None.
    etcyc : float
        Length of each ET period, T. Default is None.
    peval : float
        Potential evaporation rate (PEV) at beginning of each ET period.
        Number of entries must equal NPV, L/T. Default is None.
    rdc : float
        is an array of shape (2, npv). To conform with the sign convention
        used in most existing equations for potential evaporation, all entries
        must be greater than or equal to 0. The program multiplies all nonzero
        entries by -1 so that the evaporative flux is treated as a sink rather
        than a source.  Default is None.
        rdc(1, :) is surface resistance to evaporation (SRES) at beginning of
        ET period, L-1. For a uniform soil, SRES is equal to the reciprocal of
        the distance from the top active node to land surface, or 2/DELZ(2).
        If a surface crust is present, SRES may be decreased to account for
        the added resistance to water movement through the crust. Number of
        entries must equal NPV.
        rdc(2, :) is pressure potential of the atmosphere (HA) at beginning of
        each ET period; may be estimated using equation 6 of Lappala and
        others (1987), L. Number of entries must equal NPV.
        rdc(3, :) is Rooting depth at beginning of each ET period, L.
        Number of entries must equal NPV .
        rdc(4, :) is Root activity at base of root zone at beginning of each
        ET period, L-2. Number of entries must equal NPV.
        rdc(5, :) is Root activity at top of root zone at beginning of each
        ET period, L-2. Number of entries must equal NPV.
        Note: Values for root activity generally are determined empirically,
        but typically range from 0 to 3x104 m/m3. As programmed, root activity
        varies linearly from land surface to the base of the root zone, and
        its distribution with depth at any time is represented by a trapezoid.
        In general, root activities will be greater at land surface than at
        the base of the root zone.
        rdc(6, :) is Pressure head in roots (HROOT) at beginning of each ET
        period, L. Number of entries must equal NPV.
    ptval : float
        Potential evapotranspiration rate (PET) at beginning of each ET period,
        L/T. Number of entries must equal NPV. As with PEV, all values must be
        greater than or equal to 0. Default is None.
    iread_conc : int
        If IREAD = 0, all initial concentrations are set equal to FACTOR.
        If READ =1, all initial concentrations are read from file IU in user
        designated format and multiplied by FACTOR.. Default is None.
    factor_conc : float
        Multiplier or constant value, depending on value of IREAD, for initial
        concentrations. Default is None.
    iu_conc : int
        Unit number from which initial concentrations are to be read.
        Default is None.
    ifmt_conc : int
        Format to be used in reading initial concentrations from unit iu_conc.
        Must be enclosed in quotation marks, for example '(10X,E10.3)'.
        Default is 'free'.
    cinit : float or ndarray of size (nly, nxr)
        Initial concentrations. Default is None.
    boundary_faces : list
        boundary_faces is a list of boundary faces.  It has the following form,
        [[(l1, r1), (l2, r2), ...], [(l1, r1), (l2, r2), ...], ... ], where
        each inner list is a boundary face that itself is a list of
        (layer, column) cells that comprise that boundary face. numbf is
        equal to len(boundary_faces).
        Default is None.
    tper_tr : float or list of floats
        Length of the recharge period, T. Default is 1.0.
    delt_tr : float or list of floats
        Length of initial time step for this period, T. Default is 1.0.
    tmlt_tr : float or list of floats
        Multiplier for time step length. Default is 1.0.
    dltmx_tr : float or list of floats
        Maximum allowed length of time step, T. Default is 1.e6.
    dltmin_tr : float or list of floats
        Minimum allowed length of time step, T. Default is 1.e-5.
    tred_tr : float or list of floats
        Factor by which time-step length is reduced if convergence is not
        obtained in ITMAX iterations. Values usually should be in the range
        0.1 to 0.5. If no reduction of time-step length is desired, input a
        value of 0.0. Default is 0.1.
    dsmax_tr : float or list of floats
        Maximum allowed change in head per time step for this period, L.
        Default is 0.01.
    sterr_tr : float or list of floats
        Steady-state head criterion; when the maximum change in head between
        successive time steps is less than STERR, the program assumes that
        steady state has been reached for this period and advances to next
        recharge period, L. Default is 0.0.
    pond_tr : float or list of floats
        Maximum allowed height of ponded water for constant flux nodes. See
        Lappala and other (1987) for detailed discussion of POND, L.
        Default is 0.0.
    prnt_tr : bool or list of bools
        Logical variable = T if heads, concentration, moisture contents, and/or
        saturations are to be printed to file 6 after each time step; = F if
        they are to be written to file 6 only at observation times and ends
        of recharge periods. Default is True.
    bcit_tr : bool or list of bools
        Logical variable = T if evaporation is to be simulated for this
        recharge period; otherwise = F. Default is False.
    etsim_tr : bool or list of bools
        Logical variable = T if evapotranspiration (plant-root extraction) is
        to be simulated for this recharge period; otherwise = F.
        Default is False.
    seep_tr : bool or list of bools
        Logical variable = T if seepage faces are to be simulated for this
        recharge period; otherwise = F. Default is False.
    seepage_faces_tr: list
        seepage_faces is a list of seepage faces.  Because this is a _tr
        variable there must be nrech values in this list.  This means that the
        list of seepage faces for the first period is obtained as
        seepage_faces0 = seepage_faces_tr[0].  There may be multiple seepage
        faces for an individual period, so seepage_faces0, for example, is
        a list of seepage faces. A seepage face is (jj, jlast, seepjn), where
        jj is the number of nodes on the possible seepage face; jlast is the
        number of the node which initially represents the highest node of the
        seep; jlast can range from 0 (bottom of the face) up to JJ (top of the
        face).  seepjn is a list of (layer, column) tuples that comprise the
        seepage face.  The (layer, column) tuples must be in order from the
        lowest to the highest elevation; jj pairs of (layer, column) values
        are required.
        Default is None.
    boundary_conditions_tr : list
        boundary_conditions_tr is a list that must have nrech entries.  Each
        entry is a list of boundary conditions for that recharge period.  The
        boundary conditions for an individual recharge period is defined as
        a list boundary tuples for each boundary cell.  The tuple is defined as
        (layer, column, ntx, pdfnum, ntc, cf, ). layer and column are the
        zero-based layer and column numbers for the boundary.
        ntx is node type identifier for boundary conditions. = 0 for no
        specified boundary (needed for resetting some nodes after initial
        recharge period); = 1 for specified pressure head; = 2 for specified
        flux per unit horizontal surface area in units of L/T; = 3 for possible
        seepage face; = 4 for specified total head; = 5 for evaporation;
        = 6 for specified volumetric flow in units of L3/T; = 7 for gravity
        drain. (The gravity drain boundary condition allows gravity driven
        vertical fow out of the domain assuming a unit vertical hydraulic
        gradient. Flow into the domain cannot occur.)
        pdfnum is specified head for NTX = 1 or 4 or specified flux for
        NTX = 2 or 6. If codes 0, 3, 5, or 7 are specified, the line should
        contain a dummy value for PFDUM.
        ntc is node type identifier for transport boundary conditions.
        = 0 for no specified boundary; = 1 for specified concentration;
        cf is specified concentration for NTC = 1 or NTX = 1, 2, 4, 6, or 7.
        Present only if trans = True.
        Default is None.

    Notes
    -----
    Arguments that end with _tr are transient arguments, which mean they can
    have a different value for each recharge period.  If a _tr variable is
    entered as a scalar and has a single value, then that value is used for
    all recharge periods.

    The C records start with tper_tr.  A full set of C records are required
    for each recharge period.

    """

    def __init__(
        self,
        model,
        titl="Flopy created vs2d model",
        tmax=1.0,
        stim=0.0,
        ang=0.0,
        zunit="m",
        tunit="day",
        cunx="kg",
        nxr=1,
        nly=10,
        nrech=-1,
        numt=-100,
        rad=False,
        itstop=True,
        trans=False,
        cis=True,
        cit=False,
        sorp=False,
        f6p=True,
        thpt=True,
        spnt=True,
        ppnt=True,
        hpnt=True,
        vpnt=True,
        ifac=1,
        facx=1.0,
        dxr=None,
        xmult=None,
        xmax=None,
        jfac=1,
        facz=1.0,
        delz=None,
        zmult=None,
        zmax=None,
        pltim=None,
        obs=None,
        mb9=range(1, 73),
        eps=0.01,
        hmax=0.7,
        wus=1.0,
        eps1=0.01,
        minit=1,
        itmax=100,
        phrd=True,
        ntex=1,
        nprop=6,
        nprop1=6,
        hft=0,
        adt=1,
        aniz=(1.0,),
        hk=((1.0, 1.0e-5, 0.2, 0.0, -0.8, 2.0),),
        ht=None,
        jtex=1,
        iread=0,
        factor=0.0,
        dwtx=None,
        hmin=None,
        iu=None,
        ifmt="free",
        pinit=None,
        bcit=False,
        etsim=False,
        npv=None,
        etcyc=None,
        peval=None,
        rdc=None,
        ptval=None,
        iread_conc=None,
        factor_conc=None,
        iu_conc=None,
        ifmt_conc="free",
        cinit=None,
        boundary_faces=None,
        tper_tr=1.0,
        delt_tr=1.0,
        tmlt_tr=1.0,
        dltmx_tr=1.0e6,
        dltmin_tr=1.0e-5,
        tred_tr=0.1,
        dsmax_tr=0.01,
        sterr_tr=0.0,
        pond_tr=0.0,
        prnt_tr=True,
        bcit_tr=False,
        etsim_tr=False,
        seep_tr=False,
        seepage_faces_tr=None,
        boundary_conditions_tr=None,
        extension="dat",
        unitnumber=None,
        filenames=None,
    ):
        # set default unit number of one if not specified
        if unitnumber is None:
            unitnumber = 10

        # set filenames
        if filenames is None:
            filenames = [None]
        elif isinstance(filenames, str):
            filenames = [filenames]

        # Fill namefile items
        name = "vs2dt"
        units = None
        extra = [""]

        # set package name
        fname = [filenames[0]]

        # Call ancestor's init to set self.parent, extension, name and unit number
        Package.__init__(
            self,
            model,
            extension="dat",
            name=name,
            unit_number=units,
            extra=extra,
            filenames=fname,
        )

        self.titl = titl

        self.tmax = tmax
        self.stim = stim
        self.ang = ang

        self.zunit = zunit
        self.tunit = tunit
        self.cunx = cunx

        self.nxr = nxr
        self.nly = nly

        self.nrech = nrech
        self.numt = int(numt)

        self.rad = rad
        self.itstop = itstop
        self.trans = trans

        self.cis = cis
        self.cit = cit
        self.sorp = sorp

        self.f6p = f6p

        self.thpt = thpt
        self.spnt = spnt
        self.ppnt = ppnt
        self.hpnt = hpnt
        self.vpnt = vpnt

        self.ifac = ifac
        self.facx = facx

        self.xmult = xmult
        self.xmax = xmax
        self.dxr = dxr

        self.jfac = jfac
        self.facz = facz

        self.delz = delz
        self.zmult = zmult
        self.zmax = zmax

        self.nplt = None
        self.f8p = False
        if pltim is not None:
            self.f8p = True
            self.nplt = len(pltim)
            if self.nplt == 0:
                self.f8p = "F"
        self.pltim = pltim

        self.nobs = None
        self.f11p = False
        if obs is not None:
            self.f11p = True
            self.nobs = len(obs)
        self.obs = obs

        self.nmb9 = None
        self.f9p = False
        if mb9 is not None:
            self.f9p = True
            self.nmb9 = len(mb9)
        self.mb9 = mb9

        self.eps = eps
        self.hmax = hmax
        self.wus = wus
        self.eps1 = eps1

        self.minit = minit
        self.itmax = itmax

        self.phrd = phrd

        self.ntex = ntex
        self.nprop = nprop
        self.nprop1 = nprop1

        self.hft = hft
        self.adt = adt

        self.aniz = aniz
        self.hk = hk
        self.ht = ht

        if np.isscalar(jtex):
            jtex = np.ones((nly, nxr), dtype=np.int) * jtex
        assert jtex.shape == (nly, nxr), "{} /= {}".format(
            jtex.shape, (nly, nxr)
        )
        self.jtex = jtex

        self.iread = iread
        self.factor = factor

        self.dwtx = dwtx
        self.hmin = hmin

        self.iu = iu
        self.ifmt = ifmt
        self.pinit = pinit

        self.bcit = bcit
        self.etsim = etsim

        self.npv = npv
        self.etcyc = etcyc

        self.peval = peval
        self.rdc = rdc

        self.ptval = ptval

        self.iread_conc = iread_conc
        self.factor_conc = factor_conc

        self.iu_conc = iu_conc
        self.ifmt_conc = ifmt_conc
        self.cinit = cinit

        self.f7p = False
        if boundary_faces is not None:
            self.f7p = True
        self.boundary_faces = boundary_faces

        self.tper_tr = tper_tr
        self.delt_tr = delt_tr

        self.tmlt_tr = tmlt_tr
        self.dltmx_tr = dltmx_tr
        self.dltmin_tr = dltmin_tr
        self.tred_tr = tred_tr

        self.dsmax_tr = dsmax_tr
        self.sterr_tr = sterr_tr

        self.pond_tr = pond_tr

        self.prnt_tr = prnt_tr

        self.bcit_tr = bcit_tr
        self.etsim_tr = etsim_tr
        self.seep_tr = seep_tr

        self.seepage_faces_tr = seepage_faces_tr

        self.boundary_conditions_tr = boundary_conditions_tr
        self.parent.add_package(self)
        return

    def write_file(self):
        f = open(self.fn_path, "w")

        # A-1
        rec = "A-1 TITL"
        if self.parent.verbose:
            print("... Writing record {}...".format(rec))
        f.write("{}\n".format(self.titl))

        # A-2
        write_record(
            f,
            [self.tmax, self.stim, self.ang],
            "A-2 TMAX STIM ANG",
            verbose=self.parent.verbose,
        )

        # A-3
        rec = "A-3 ZUNIT TUNIT CUNX"
        if self.parent.verbose:
            print("... Writing record {}...".format(rec))
        f.write(
            "{:>4s}{:>4s}{:>4s}\n".format(self.zunit, self.tunit, self.cunx)
        )

        # A-4
        write_record(
            f, [self.nxr, self.nly], "A-4 NXR NLY", verbose=self.parent.verbose
        )

        # A-5
        write_record(
            f,
            [self.nrech, self.numt],
            "A-5 NRECH NUMT",
            verbose=self.parent.verbose,
        )

        # A-6
        rec = "A-6 RAD ITSTOP TRANS"
        if self.parent.verbose:
            print("... Writing record {}...".format(rec))
        fmt = "{:>2s}{:>2s}{:>2s}" + "  {}" + "\n"
        f.write(
            fmt.format(
                bool2str(self.rad),
                bool2str(self.itstop),
                bool2str(self.trans),
                rec,
            )
        )

        # A-6A
        if self.trans:
            rec = "A-6A CIS CIT SORP"
            if self.parent.verbose:
                print("... Writing record {}...".format(rec))
            fmt = "{:>2s}{:>2s}{:>2s}" + "  {}" + "\n"
            f.write(
                fmt.format(
                    bool2str(self.cis),
                    bool2str(self.cit),
                    bool2str(self.sorp),
                    rec,
                )
            )

        # A-7
        rec = "A-7 F11P F7P F8P F9P F6P"
        if self.parent.verbose:
            print("... Writing record {}...".format(rec))
        fmt = 5 * "{:>2s}" + "  {}" + "\n"
        f.write(
            fmt.format(
                bool2str(self.f11p),
                bool2str(self.f7p),
                bool2str(self.f8p),
                bool2str(self.f9p),
                bool2str(self.f6p),
                rec,
            )
        )

        # A-8
        rec = "A-8 THPT SPNT PPNT HPNT VPNT"
        if self.parent.verbose:
            print("... Writing record {}...".format(rec))
        fmt = 5 * "{:>2s}" + "  {}" + "\n"
        f.write(
            fmt.format(
                bool2str(self.thpt),
                bool2str(self.spnt),
                bool2str(self.ppnt),
                bool2str(self.hpnt),
                bool2str(self.vpnt),
                rec,
            )
        )

        # A-9
        write_record(
            f,
            [self.ifac, self.facx],
            "A-9 IFAC FACX",
            verbose=self.parent.verbose,
        )

        # A-10
        if self.ifac == 0:
            rec = "A-10 DXR"
            if self.parent.verbose:
                print("... Writing record {}...".format(rec))
            dxrlist = [str(val) for val in self.dxr]
            f.write("{}\n".format(" ".join(dxrlist)))
        elif self.ifac == 2:
            rec = "A-10 XMULT XMAX"
            write_record(
                f, [self.xmult, self.xmax], rec, verbose=self.parent.verbose
            )

        # A-11
        write_record(
            f,
            [self.jfac, self.facz],
            "A-11 JFAC FACZ",
            verbose=self.parent.verbose,
        )

        # A-12
        if self.jfac == 0:
            rec = "A-12 DELZ"
            if self.parent.verbose:
                print("... Writing record {}...".format(rec))
            delzlist = [str(val) for val in self.delz]
            f.write("{}\n".format(" ".join(delzlist)))
        elif self.jfac == 2:
            write_record(
                f,
                [self.zmult, self.zmax],
                "A-12 ZMULT ZMAX",
                verbose=self.parent.verbose,
            )

        # A-13
        if self.f8p:
            write_record(
                f, [self.nplt], "A-13 NPLT", verbose=self.parent.verbose
            )

        # A-14
        if self.f8p:
            rec = "A-14 PLTIM"
            if self.parent.verbose:
                print("... Writing record {}...".format(rec))
            pltimlist = [str(val) for val in self.pltim]
            f.write("{}\n".format(" ".join(pltimlist)))

        # A-15
        if self.f11p:
            write_record(
                f, [self.nobs], "A-15 NOBS", verbose=self.parent.verbose
            )

        # A-16
        if self.f11p:
            rec = "A-16 J N (obs)"
            if self.parent.verbose:
                print("... Writing record {}...".format(rec))
            for row in self.obs:
                ir = row[0] + 1
                ic = row[1] + 1
                f.write("{} {}\n".format(ir, ic))

        # A-17
        if self.f9p:
            write_record(
                f, [self.nmb9], "A-17 NMB9", verbose=self.parent.verbose
            )

        # A-18
        if self.f9p:
            rec = "A-18 MB9"
            if self.parent.verbose:
                print("... Writing record {}...".format(rec))
            mb9list = [str(val) for val in self.mb9]
            f.write("{}\n".format(" ".join(mb9list)))

        # B-1
        write_record(
            f,
            [self.eps, self.hmax, self.wus, self.eps1],
            "B-1 EPS HMAX WUS EPS1",
            verbose=self.parent.verbose,
        )

        # B-3
        write_record(
            f,
            [self.minit, self.itmax],
            "B-3 MINIT ITMAX",
            verbose=self.parent.verbose,
        )

        # B-4
        write_record(
            f, [bool2str(self.phrd)], "B-4 PHRD", verbose=self.parent.verbose
        )

        rec = "B-5 NTEX NPROP"
        vars = [self.ntex, self.nprop]
        if self.trans:
            vars.append(self.nprop1)
            rec += " NPROP1"
        write_record(f, vars, rec, verbose=self.parent.verbose)

        rec = "B-5A HFT ADT"
        vars = [self.hft, self.adt]
        write_record(f, vars, rec, verbose=self.parent.verbose)

        rec = "B-6 B-7 B-7A"
        for itex in range(self.ntex):
            write_record(
                f, [itex + 1], "B-6: itex", verbose=self.parent.verbose
            )
            write_record(
                f,
                [self.aniz[itex]] + list(self.hk[itex]),
                "B-7: aniz, hk",
                verbose=self.parent.verbose,
            )
            if self.trans:
                write_record(
                    f,
                    list(self.ht[itex]),
                    "B-7A: ht",
                    verbose=self.parent.verbose,
                )

        rec = "B-8 IROW"
        vars = [0]
        write_record(f, vars, rec, verbose=self.parent.verbose)

        rec = "B-9 JTEX"
        for k in range(self.nly):
            write_record(
                f,
                self.jtex[k],
                "{} Layer {}".format(rec, k + 1),
                verbose=self.parent.verbose,
            )

        rec = "B-11 IREAD FACTOR"
        vars = [self.iread, self.factor]
        write_record(f, vars, rec, verbose=self.parent.verbose)

        if self.iread == 2:
            rec = "B-12 DWTX HMIN"
            vars = [self.dwtx, self.hmin]
            write_record(f, vars, rec, verbose=self.parent.verbose)

        elif self.iread == 1:
            rec = "B-13 IU IFMT"
            vars = [self.iu, self.ifmt]
            write_record(f, vars, rec, verbose=self.parent.verbose)
            if self.iu == 5:
                for k in range(self.nly):
                    write_record(
                        f,
                        self.pinit[k],
                        "{} Layer {}".format(rec, k + 1),
                        verbose=self.parent.verbose,
                    )

        rec = "B-14 BCIT ETSIM"
        vars = [bool2str(self.bcit), bool2str(self.etsim)]
        write_record(f, vars, rec, verbose=self.parent.verbose)

        if self.bcit or self.etsim:
            rec = "B-15 NPV ETCYC"
            vars = [self.npv, self.etcyc]
            write_record(f, vars, rec, verbose=self.parent.verbose)

        if self.bcit:
            rec = "B-16 PEVAL"
            write_record(f, self.peval, rec, verbose=self.parent.verbose)
            rec = "B-17 RDC(1) SRES"
            write_record(f, self.rdc[0], rec, verbose=self.parent.verbose)
            rec = "B-18 RDC(2) HA"
            write_record(f, self.rdc[1], rec, verbose=self.parent.verbose)

        if self.etsim:
            rec = "B-19 PTVAL"
            write_record(f, self.ptval, rec, verbose=self.parent.verbose)
            rec = "B-20 RDC(3) ROOTING DEPTH"
            write_record(f, self.rdc[2], rec, verbose=self.parent.verbose)
            rec = "B-21 RDC(4) ROOT ACTIVITY AT BASE"
            write_record(f, self.rdc[3], rec, verbose=self.parent.verbose)
            rec = "B-22 RDC(5) ROOT ACTIVITY AT TOP"
            write_record(f, self.rdc[4], rec, verbose=self.parent.verbose)
            rec = "B-23 RDC(6) HROOT"
            write_record(f, self.rdc[5], rec, verbose=self.parent.verbose)

        if self.trans:
            rec = "B-24 IREAD_CONC FACTOR_CONC"
            vars = [self.iread_conc, self.factor_conc]
            write_record(f, vars, rec, verbose=self.parent.verbose)

            if self.iread_conc == 1:
                rec = "B-25 IU_CONC IFMT_CONC"
                vars = [self.iu_conc, self.ifmt_conc]
                write_record(f, vars, rec, verbose=self.parent.verbose)
                if self.iu_conc == 5:
                    for k in range(self.nly):
                        write_record(
                            f,
                            self.cinit[k],
                            "{} Layer {}".format(rec, k + 1),
                            verbose=self.parent.verbose,
                        )

        if self.f7p:
            numbf = len(self.boundary_faces)
            maxcells = 0
            for bf in self.boundary_faces:
                maxcells = max(maxcells, len(bf))
            rec = "B-26 NUMBF MAXCELLS"
            vars = [numbf, maxcells]
            write_record(f, vars, rec, verbose=self.parent.verbose)
            for ibf, bf in enumerate(self.boundary_faces):
                rec = "B-27 IDBF NUMCELLS"
                vars = [ibf + 1, len(bf)]
                write_record(f, vars, rec, verbose=self.parent.verbose)
                for irow, icol in bf:
                    rec = "B-28 J,N"
                    write_record(
                        f,
                        [irow + 1, icol + 1],
                        rec,
                        verbose=self.parent.verbose,
                    )

        for irech in range(abs(self.nrech)):
            rec = "C-1 TPER DELT"
            vars = [self.tper_tr, self.delt_tr]
            write_record(
                f, vars, rec, irech=irech, verbose=self.parent.verbose
            )

            rec = "C-2 TMLT DLTMX DLTMIN TRED"
            vars = [self.tmlt_tr, self.dltmx_tr, self.dltmin_tr, self.tred_tr]
            write_record(
                f, vars, rec, irech=irech, verbose=self.parent.verbose
            )

            rec = "C-3 DSMAX STERR"
            vars = [self.dsmax_tr, self.sterr_tr]
            write_record(
                f, vars, rec, irech=irech, verbose=self.parent.verbose
            )

            rec = "C-4 POND"
            vars = [self.pond_tr]
            write_record(
                f, vars, rec, irech=irech, verbose=self.parent.verbose
            )

            rec = "C-5 PRNT"
            vars = [self.prnt_tr]
            write_record(
                f, vars, rec, irech=irech, verbose=self.parent.verbose
            )

            rec = "C-6 BCIT ETSIM SEEP"
            vars = [self.bcit_tr, self.etsim_tr, self.seep_tr]
            write_record(
                f, vars, rec, irech=irech, verbose=self.parent.verbose
            )

            seep = self.seep_tr
            if isinstance(self.seep_tr, list):
                seep = seep[irech]
            if seep:
                seepage_faces = self.seepage_faces_tr[irech]
                rec = "C-7 NFCS"
                f.write("{} {}\n".format(len(seepage_faces), rec))
                for jj, jlast, seepjn in seepage_faces:
                    rec = "C-8 JJ JLAST"
                    f.write("{} {} {}\n".format(jj, jlast, rec))
                    for j, n in seepjn:
                        f.write("{} {}\n".format(j + 1, n + 1))

            rec = "C-10 IBC"
            vars = [0]
            write_record(
                f, vars, rec, irech=irech, verbose=self.parent.verbose
            )

            if self.boundary_conditions_tr is not None:
                boundary_conditions = self.boundary_conditions_tr[irech]
                for bc in boundary_conditions:
                    jj = bc.pop(0)
                    nn = bc.pop(0)
                    ntx = bc.pop(0)
                    pdfnum = bc.pop(0)
                    s = "{} {} {} {}".format(jj + 1, nn + 1, ntx, pdfnum)
                    if len(bc) > 1:
                        ntc = bc.pop(0)
                        cf = bc.pop(0)
                        s += " {} {}".format(ntc, cf)
                    s += "\n"
                    f.write(s)
            f.write("-999999 /\n")

        f.write("-999999 / C-13\n")
        f.close()
        return

    @staticmethod
    def load(f, model, ext_unit_dict=None):
        # open if f is a file name
        if not hasattr(f, "read"):
            filename = f
            f = open(filename, "r")

        # vardict
        vardict = {}

        # A-1
        titl = None
        if model.verbose:
            print("   loading A-1...")
        line = f.readline()
        titl = line.strip()
        if model.verbose:
            print("A-1: {}".format(titl))
        vardict["titl"] = titl

        # A-2
        # tmax = None
        # stim = None
        # ang = None
        # if model.verbose:
        #     print('   loading A-2...')
        # line = f.readline()
        # t = line_parse(line)
        # tmax = float(t[0])
        # stim = float(t[1])
        # ang = float(t[2])
        # if model.verbose:
        #     print('A-2: {} {} {}'.format(tmax, stim, ang))
        # vardict['tmax'] = tmax
        # vardict['stim'] = stim
        # vardict['ang'] = ang

        rec = "A-2"
        varnames = ["tmax", "stim", "ang"]
        vartypes = [float, float, float]
        vardict = load_record(
            f, rec, varnames, vartypes, vardict, model.verbose
        )

        rec = "A-3"
        zunit = None
        tunit = None
        cunx = None
        if model.verbose:
            print("   loading {}...".format(rec))
        line = f.readline()
        t = line_parse(line)
        zunit = t[0]
        tunit = t[1]
        cunx = t[2]
        if model.verbose:
            print("{}: {} {} {}".format(rec, zunit, tunit, cunx))
        vardict["zunit"] = zunit
        vardict["tunit"] = tunit
        vardict["cunx"] = cunx

        rec = "A-4"
        nxr = None
        nly = None
        if model.verbose:
            print("   loading {}...".format(rec))
        line = f.readline()
        t = line_parse(line)
        nxr = int(t[0])
        nly = int(t[1])
        if model.verbose:
            print("{}: {} {}".format(rec, nxr, nly))
        vardict["nxr"] = nxr
        vardict["nly"] = nly

        rec = "A-5"
        nrech = None
        numt = None
        if model.verbose:
            print("   loading {}...".format(rec))
        line = f.readline()
        t = line_parse(line)
        nrech = int(t[0])
        numt = int(t[1])
        if model.verbose:
            print("{}: {} {}".format(rec, nrech, numt))
        vardict["nrech"] = nrech
        vardict["numt"] = numt

        rec = "A-6"
        rad = None
        itstop = None
        trans = None
        if model.verbose:
            print("   loading {}...".format(rec))
        line = f.readline()
        t = line_parse(line)
        rad = t[0] == "T"
        itstop = t[1] == "T"
        trans = t[2] == "T"
        if model.verbose:
            print("{}: {} {} {}".format(rec, rad, itstop, trans))
        vardict["rad"] = rad
        vardict["itstop"] = itstop
        vardict["trans"] = trans

        rec = "A-6A"
        cis = None
        cit = None
        sorp = None
        if trans:
            if model.verbose:
                print("   loading {}...".format(rec))
            line = f.readline()
            t = line_parse(line)
            cis = t[0] == "T"
            cit = t[1] == "T"
            sorp = t[2] == "T"
            if model.verbose:
                print("{}: {} {} {}".format(rec, cis, cit, sorp))
        vardict["cis"] = cis
        vardict["cit"] = cit
        vardict["sorp"] = sorp

        rec = "A-7"
        variablelist = ["f11p", "f7p", "f8p", "f9p", "f6p"]
        if model.verbose:
            print("{}: {}".format(rec, " ".join(variablelist)))
        line = f.readline()
        t = line_parse(line)
        for i, variable in enumerate(variablelist):
            value = t[i] == "T"  # convert to bool
            vardict[variable] = value
        if model.verbose:
            s = ""
            for variable in variablelist:
                s += bool2str(vardict[variable]) + " "
            print("{}: {}".format(rec, s))

        rec = "A-8"
        variablelist = ["thpt", "spnt", "ppnt", "hpnt", "vpnt"]
        if model.verbose:
            print("{}: {}".format(rec, " ".join(variablelist)))
        line = f.readline()
        t = line_parse(line)
        for i, variable in enumerate(variablelist):
            value = t[i] == "T"
            vardict[variable] = value
        if model.verbose:
            s = ""
            for variable in variablelist:
                s += bool2str(vardict[variable]) + " "
            print("{}: {}".format(rec, s))

        rec = "A-9"
        ifac = None
        facx = None
        if model.verbose:
            print("   loading {}...".format(rec))
        line = f.readline()
        t = line_parse(line)
        ifac = int(t[0])
        facx = float(t[1])
        if model.verbose:
            print("{}: {} {}".format(rec, ifac, facx))
        vardict["ifac"] = ifac
        vardict["facx"] = facx

        rec = "A-10"
        dxr = None
        xmult = None
        xmax = None
        if ifac == 0:
            dxr = np.empty(nxr, dtype=np.float)
            read1d(f, dxr)
        elif ifac == 2:
            line = f.readline()
            t = line_parse(line)
            xmult = float(t[0])
            xmax = float(t[1])
        if model.verbose:
            print("{}: {} {} {}".format(rec, dxr, xmult, xmax))
        vardict["dxr"] = dxr
        vardict["xmult"] = xmult
        vardict["xmax"] = xmax

        rec = "A-11"
        jfac = None
        facz = None
        if model.verbose:
            print("   loading {}...".format(rec))
        line = f.readline()
        t = line_parse(line)
        jfac = int(t[0])
        facz = float(t[1])
        if model.verbose:
            print("{}: {} {}".format(rec, jfac, facz))
        vardict["jfac"] = jfac
        vardict["facz"] = facz

        rec = "A-12"
        delz = None
        zmult = None
        zmax = None
        if jfac == 0:
            delz = np.empty(nly, dtype=np.float)
            read1d(f, delz)
        elif jfac == 2:
            line = f.readline()
            t = line_parse(line)
            zmult = float(t[0])
            zmax = float(t[1])
        if model.verbose:
            print("{}: {} {} {}".format(rec, delz, zmult, zmax))
        vardict["delz"] = delz
        vardict["zmult"] = zmult
        vardict["zmax"] = zmax

        rec = "A-13"
        f8p = vardict.pop("f8p")
        nplt = None
        if f8p:
            line = f.readline()
            t = line_parse(line)
            nplt = int(t[0])
        if model.verbose:
            print("{}: {}".format(rec, nplt))

        rec = "A-14"
        pltim = None
        if f8p:
            pltim = np.empty(nplt, dtype=np.float)
            read1d(f, pltim)
        if model.verbose:
            print("{}: {}".format(rec, pltim))
        vardict["pltim"] = pltim

        rec = "A-15"
        f11p = vardict.pop("f11p")
        nobs = None
        if f11p:
            line = f.readline()
            t = line_parse(line)
            nobs = int(t[0])
        if model.verbose:
            print("{}: {}".format(rec, nobs))

        rec = "A-16"
        obs = None
        if f11p:
            obs = np.empty((abs(nobs), 2), dtype=np.int)
            for iobs in range(abs(nobs)):
                line = f.readline()
                t = line_parse(line)
                obs[iobs, 0] = int(t[0]) - 1
                obs[iobs, 1] = int(t[1]) - 1
        if model.verbose:
            print("{}: {}".format(rec, obs))
        vardict["obs"] = obs

        rec = "A-17"
        f9p = vardict.pop("f9p")
        nmb9 = None
        if f9p:
            line = f.readline()
            t = line_parse(line)
            nmb9 = int(t[0])
        if model.verbose:
            print("{}: {}".format(rec, nmb9))

        rec = "A-18"
        mb9 = None
        if f9p:
            mb9 = np.empty(nmb9, dtype=np.int)
            read1d(f, mb9)
        if model.verbose:
            print("{}: {}".format(rec, mb9))
        vardict["mb9"] = mb9

        # B-1
        rec = "B-1"
        eps = None
        hmax = None
        wus = None
        eps1 = None
        if model.verbose:
            print("   loading {}...".format(rec))
        line = f.readline()
        t = line_parse(line)
        eps = float(t[0])
        hmax = float(t[1])
        wus = float(t[2])
        eps1 = float(t[3])
        if model.verbose:
            print("{}: {} {} {} {}".format(rec, eps, hmax, wus, eps1))
        vardict["eps"] = eps
        vardict["hmax"] = hmax
        vardict["wus"] = wus
        vardict["eps1"] = eps1

        # B-3
        rec = "B-3"
        minit = None
        itmax = None
        if model.verbose:
            print("   loading {}...".format(rec))
        line = f.readline()
        t = line_parse(line)
        minit = int(t[0])
        itmax = int(t[1])
        if model.verbose:
            print("{}: {} {}".format(rec, minit, itmax))
        vardict["minit"] = minit
        vardict["itmax"] = itmax

        rec = "B-4"
        phrd = None
        if model.verbose:
            print("   loading {}...".format(rec))
        line = f.readline()
        t = line_parse(line)
        phrd = t[0] == "T"
        if model.verbose:
            print("{}: {}".format(rec, phrd))
        vardict["phrd"] = phrd

        rec = "B-5"
        varnames = ["ntex", "nprop"]
        vartypes = [int, int]
        if trans:
            varnames.append("nprop1")
            vartypes.append(int)
        vardict = load_record(
            f, rec, varnames, vartypes, vardict, model.verbose
        )

        rec = "B-5A"
        varnames = ["hft"]
        vartypes = [int]
        if trans:
            varnames.append("adt")
            vartypes.append(int)
        vardict = load_record(
            f, rec, varnames, vartypes, vardict, model.verbose
        )

        # B-6, B-7, B-7A
        if model.verbose:
            print("   loading {}...".format("B-6, B-7, B-7A"))
        ntex = vardict["ntex"]
        nprop = vardict["nprop"]
        aniz = np.empty(ntex, dtype=np.float)
        hk = np.empty((ntex, nprop), dtype=np.float)
        ht = None
        if trans:
            ht = np.empty((ntex, 8), dtype=np.float)
        for itex in range(ntex):
            # B-6
            line = f.readline()
            linelist = line_parse(line)
            itex_read = int(linelist[0])
            assert itex == itex_read - 1

            # B-7
            line = f.readline()
            linelist = line_parse(line)
            aniz[itex] = float(linelist.pop(0))
            for iprop in range(nprop):
                hk[itex, iprop] = float(linelist.pop(0))

            # B-7A
            if trans:
                line = f.readline()
                linelist = line_parse(line)
                for iprop in range(6):
                    ht[itex, iprop] = float(linelist.pop(0))

        vardict["aniz"] = aniz
        vardict["hk"] = hk
        vardict["ht"] = ht

        rec = "B-8 B-9 B-10"
        if model.verbose:
            print("   loading {}...".format(rec))
        line = f.readline()
        t = line_parse(line)
        irow = int(t[0])
        if model.verbose:
            print("{}: {}".format("B8", irow))
        if irow == 0:
            jtex = np.zeros((nly * nxr), dtype=np.int)
            read1d(f, jtex)
            jtex = jtex.reshape((nly, nxr))
        elif irow == 1:
            jtex = np.zeros((nly, nxr), dtype=np.int)
            ilay = 0
            while True:
                line = f.readline()
                t = line_parse(line)
                ileft = int(t.pop(0)) - 1
                iright = int(t.pop(0)) - 1
                jbt = int(t.pop(0)) - 1
                jrd = int(t.pop(0))
                jtex[ilay : jbt + 1, ileft : iright + 1] = jrd
                ilay = jbt
                if iright + 1 == nxr and jbt + 1 == nly:
                    break
        else:
            assert False, "unknown irow value {}".format(irow)
        vardict["jtex"] = jtex

        rec = "B-11"
        varnames = ["iread", "factor"]
        vartypes = [int, float]
        vardict = load_record(
            f, rec, varnames, vartypes, vardict, model.verbose
        )

        iread = vardict["iread"]
        if iread == 2:
            rec = "B-12"
            varnames = ["dwtx", "hmin"]
            vartypes = [float, float]
            vardict = load_record(
                f, rec, varnames, vartypes, vardict, model.verbose
            )
        elif iread == 1:
            rec = "B-13"
            line = f.readline()
            linelist = line.strip().split()
            iu = int(linelist.pop(0))
            vardict["iu"] = iu
            ifmt = line.strip().split("'")[1]
            vardict["ifmt"] = ifmt

        rec = "B-14"
        varnames = ["bcit", "etsim"]
        vartypes = [bool, bool]
        vardict = load_record(
            f, rec, varnames, vartypes, vardict, model.verbose
        )

        bcit = vardict["bcit"]
        etsim = vardict["etsim"]
        if bcit or etsim:
            rec = "B-15"
            varnames = ["npv", "etcyc"]
            vartypes = [int, float]
            vardict = load_record(
                f, rec, varnames, vartypes, vardict, model.verbose
            )

        rdc = None
        if bcit or etsim:
            npv = vardict["npv"]
            rdc = np.empty((6, npv), dtype=np.float)

        if bcit:
            rec = "B-16"
            npv = vardict["npv"]
            peval = np.empty(npv, dtype=np.float)
            read1d(f, peval)
            vardict["peval"] = peval

            rec = "B-17"
            read1d(f, rdc[0])
            rec = "B-18"
            read1d(f, rdc[1])

        if etsim:
            rec = "B-19"
            ptval = np.empty(npv, dtype=np.float)
            read1d(f, ptval)
            vardict["ptval"] = ptval
            rec = "B-20"
            read1d(f, rdc[2])
            rec = "B-21"
            read1d(f, rdc[3])
            rec = "B-22"
            read1d(f, rdc[4])
            rec = "B-23"
            read1d(f, rdc[5])

        vardict["rdc"] = rdc

        if trans and iread != 3:
            rec = "B-24"
            varnames = ["iread_conc", "factor_conc"]
            vartypes = [int, float]
            vardict = load_record(
                f, rec, varnames, vartypes, vardict, model.verbose
            )
            if vardict["iread_conc"] == 1:
                rec = "B-25"
                line = f.readline()
                linelist = line.strip().split()
                iu_conc = int(linelist.pop(0))
                vardict["iu_conc"] = iu_conc
                ifmt_conc = line.strip().split("'")[1]
                vardict["ifmt_conc"] = ifmt_conc
                if iu_conc == 5:
                    cinit = np.empty(nly * nxr, dtype=np.float)
                    read1d(f, cinit)
                    cinit = cinit.reshape((nly, nxr))
                    vardict["cinit"] = cinit

        f7p = vardict.pop("f7p")
        if f7p:
            rec = "B-26"
            varnames = ["numbf", "maxcells"]
            vartypes = [int, int]
            vardict = load_record(
                f, rec, varnames, vartypes, vardict, model.verbose
            )

            numbf = vardict["numbf"]
            boundary_faces = []
            for nface in range(numbf):
                line = f.readline()
                linelist = line.strip().split()
                idbf = int(linelist.pop(0))
                numcells = int(linelist.pop(0))
                facelist = []
                for icell in range(numcells):
                    line = f.readline()
                    linelist = line.strip().split()
                    irow = int(linelist.pop(0)) - 1
                    icol = int(linelist.pop(0)) - 1
                    facelist.append((irow, icol))
                boundary_faces.append(facelist)
            vardict["boundary_faces"] = boundary_faces

        for irech in range(nrech):
            rec = "C-1"
            varnames = ["tper_tr", "delt_tr"]
            vartypes = [float, float]
            vardict = load_record(
                f, rec, varnames, vartypes, vardict, model.verbose
            )

            rec = "C-2"
            varnames = ["tmlt_tr", "dltmx_tr", "dltmin_tr", "tred_tr"]
            vartypes = 4 * [float]
            vardict = load_record(
                f, rec, varnames, vartypes, vardict, model.verbose
            )

            rec = "C-3"
            varnames = ["dsmax_tr", "sterr_tr"]
            vartypes = 2 * [float]
            vardict = load_record(
                f, rec, varnames, vartypes, vardict, model.verbose
            )

            rec = "C-4"
            varnames = ["pond_tr"]
            vartypes = [float]
            vardict = load_record(
                f, rec, varnames, vartypes, vardict, model.verbose
            )

            rec = "C-5"
            varnames = ["prnt_tr"]
            vartypes = [bool]
            vardict = load_record(
                f, rec, varnames, vartypes, vardict, model.verbose
            )

            rec = "C-6"
            varnames = ["bcit_tr", "etsim_tr", "seep_tr"]
            vartypes = [bool, bool, bool]
            vardict = load_record(
                f, rec, varnames, vartypes, vardict, model.verbose
            )

            # C-7 C-8 C-9
            if irech == 0:
                vardict["seepage_faces_tr"] = []
            seepage_faces = None
            seep = vardict["seep_tr"]
            if isinstance(seep, list):
                seep = seep[-1]
            if seep:
                rec = "C-7"
                line = f.readline()
                t = line_parse(line)
                nfcs = int(t.pop(0))

                # nfcs is equal to the length of the seepage_faces list
                seepage_faces = []
                for ifcs in range(nfcs):
                    rec = "C-8"
                    line = f.readline()
                    t = line_parse(line)
                    jj = int(t.pop(0))
                    jlast = int(t.pop(0))

                    rec = "C-9"
                    seepjn = []
                    for iseep in range(jj):
                        line = f.readline()
                        t = line_parse(line)
                        j = int(t.pop(0)) - 1
                        n = int(t.pop(0)) - 1
                        seepjn.append((j, n))

                    seepage_faces.append((jj, jlast, seepjn))
            vardict["seepage_faces_tr"].append(seepage_faces)

            # C-11 C-12
            if irech == 0:
                vardict["boundary_conditions_tr"] = []
            rec = "C-10"
            line = f.readline()
            t = line_parse(line)
            ibc = int(t.pop(0))

            boundary_conditions = []
            if ibc == 0:
                while True:
                    line = f.readline()
                    if line.strip().startswith("-999999"):
                        break

                    linelist = line.strip().split()
                    jj = int(linelist.pop(0)) - 1
                    nn = int(linelist.pop(0)) - 1
                    ntx = int(linelist.pop(0))
                    pdfnum = float(linelist.pop(0))
                    bc = [jj, nn, ntx, pdfnum]

                    if trans:
                        ntc = int(linelist.pop(0))
                        cf = float(linelist.pop(0))
                        bc += [ntc, cf]

                    boundary_conditions.append(bc)

            elif ibc == 1:
                assert False, "ibc == 1 not implemented yet"

            else:
                assert False, "Invalid value for ibc {}".format(ibc)

            vardict["boundary_conditions_tr"].append(boundary_conditions)

        f.close()
        p = Vs2dt(model, **vardict)
        return p

    def get_delr(self):
        delr = np.empty(self.nxr, dtype=np.float)
        if self.ifac == 0:
            delr[:] = self.dxr * self.facx
        elif self.ifac == 1:
            delr[:] = self.facx
        elif self.ifac == 2:
            delr[:] = self.facx
            for j in range(2, self.nxr):
                delr[j] = min(delr[j - 1] * self.xmult, self.xmax)
        return delr

    def get_dell(self):
        dell = np.empty(self.nly, dtype=np.float)
        if self.jfac == 0:
            dell[:] = self.delz * self.facz
        elif self.jfac == 1:
            dell[:] = self.facz
        elif self.jfac == 2:
            dell[:] = self.facz
            for k in range(2, self.nly):
                dell[k] = min(dell[k - 1] * self.zmult, self.zmax)
        return dell
