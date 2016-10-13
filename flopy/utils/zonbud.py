from __future__ import print_function
import os
import numpy as np
from .binaryfile import CellBudgetFile


class Budget(object):
    """
    ZoneBudget Budget class. This is a wrapper around a numpy record array to allow users
    to save the record array to a formatted csv file.
    """
    def __init__(self, records, **kwargs):
        self.records = records
        self.kwargs = kwargs
        self._zonefields = [name for name in self.records.dtype.names if 'ZONE' in name]

    def get_total_inflow(self):
        idx = np.where(self.records['flow_dir'] == 'in')[0]
        ins = _numpyvoid2numeric(self.records[self._zonefields][idx])
        return ins.sum(axis=0)

    def get_total_outflow(self):
        idx = np.where(self.records['flow_dir'] == 'out')[0]
        out = _numpyvoid2numeric(self.records[self._zonefields][idx])
        return out.sum(axis=0)

    def get_percent_error(self):
        ins_minus_out = self.get_total_inflow() - self.get_total_outflow()
        ins_plus_out = self.get_total_inflow() + self.get_total_outflow()
        pcterr = 100 * ins_minus_out / (ins_plus_out / 2.)
        pcterr = np.nan_to_num(pcterr)
        return pcterr

    def to_csv(self, fname, write_format='pandas', formatter=None):
        """
        Saves the Budget object record array to a formatted csv file.

        Parameters
        ----------
        fname
        write_format
        formatter

        Returns
        -------

        """
        assert write_format.lower() in ['pandas', 'zonbud'], 'Format must be one of "pandas" or "zonbud".'

        if formatter is None:
            formatter = '{:.16e}'.format

        if write_format.lower() == 'pandas':
            with open(fname, 'w') as f:

                # Write header
                f.write(','.join(self.records.dtype.names)+'\n')

                # Write IN terms
                for rec in self.records[np.where(self.records['flow_dir'] == 'in')[0]]:
                    items = []
                    for i in rec:
                        if isinstance(i, str):
                            items.append(i)
                        else:
                            items.append(formatter(i))
                    f.write(','.join(items)+'\n')
                ins_sum = self.get_total_inflow()
                f.write(','.join([' ', 'Total IN'] + [formatter(i) for i in ins_sum])+'\n')

                # Write OUT terms
                for rec in self.records[np.where(self.records['flow_dir'] == 'out')[0]]:
                    items = []
                    for i in rec:
                        if isinstance(i, str):
                            items.append(i)
                        else:
                            items.append(formatter(i))
                    f.write(','.join(items) + '\n')
                out_sum = self.get_total_outflow()
                f.write(','.join([' ', 'Total OUT'] + [formatter(i) for i in out_sum])+'\n')

                # Write mass balance terms
                ins_minus_out = self.get_total_inflow() - self.get_total_outflow()
                pcterr = self.get_percent_error()
                f.write(','.join([' ', 'IN-OUT'] + [formatter(i) for i in ins_minus_out])+'\n')
                f.write(','.join([' ', 'Percent Error'] + [formatter(i) for i in pcterr])+'\n')

        elif write_format.lower() == 'zonbud':
            with open(fname, 'w') as f:

                # Write header
                if 'kstpkper' in self.kwargs.keys():
                    header = 'Time Step, {kstp}, Stress Period, {kper}\n'.format(kstp=self.kwargs['kstpkper'][0]+1,
                                                                                 kper=self.kwargs['kstpkper'][1]+1)
                elif 'totim' in self.kwargs.keys():
                    header = 'Sim. Time, {totim}\n'.format(totim=self.kwargs['totim'])
                else:
                    raise Exception('No stress period/time step or time specified.')

                f.write(header)
                f.write(','.join([' '] + [field for field in self.records.dtype.names[2:]])+'\n')

                # Write IN terms
                f.write(','.join([' '] + ['IN']*(len(self.records.dtype.names[1:])-1))+'\n')
                for rec in self.records[np.where(self.records['flow_dir'] == 'in')[0]]:
                    items = []
                    for i in list(rec)[1:]:
                        if isinstance(i, str):
                            items.append(i)
                        else:
                            items.append(formatter(i))
                    f.write(','.join(items)+'\n')
                ins_sum = self.get_total_inflow()
                f.write(','.join(['Total IN'] + [formatter(i) for i in ins_sum])+'\n')

                # Write OUT terms
                f.write(','.join([' '] + ['OUT']*(len(self.records.dtype.names[1:])-1))+'\n')
                for rec in self.records[np.where(self.records['flow_dir'] == 'out')[0]]:
                    items = []
                    for i in list(rec)[1:]:
                        if isinstance(i, str):
                            items.append(i)
                        else:
                            items.append(formatter(i))
                    f.write(','.join(items) + '\n')
                out_sum = self.get_total_outflow()
                f.write(','.join(['Total OUT'] + [formatter(i) for i in out_sum])+'\n')

                # Write mass balance terms
                ins_minus_out = self.get_total_inflow() - self.get_total_outflow()
                pcterr = self.get_percent_error()
                f.write(','.join(['IN-OUT'] + [formatter(i) for i in ins_minus_out])+'\n')
                f.write(','.join(['Percent Error'] + [formatter(i) for i in pcterr])+'\n')
        return


class ZoneBudget(object):
    """
    ZoneBudget class

    Example usage:

    >>>from flopy.utils import ZoneBudget
    >>>zb = ZoneBudget('zonebudtest.cbc')
    >>>bud = zb.get_budget('GWBasins.zon', kstpkper=(0, 0))
    >>>bud.to_csv('zonebudtest.csv')
    """
    def __init__(self, cbc_file):

        if isinstance(cbc_file, CellBudgetFile):
            self.cbc = cbc_file
        elif isinstance(cbc_file, str) and os.path.isfile(cbc_file):
            self.cbc = CellBudgetFile(cbc_file)
        else:
            raise Exception('Cannot load cell budget file: {}.'.format(cbc_file))

        # All record names in the cell-by-cell budget binary file
        self.record_names = [n.strip() for n in self.cbc.unique_record_names()]

        # Get imeth for each record in the CellBudgetFile record list
        self.imeth = {}
        for record in self.cbc.recordarray:
            self.imeth[record['text'].strip()] = record['imeth']

        # INTERNAL FLOW TERMS ARE USED TO CALCULATE FLOW BETWEEN ZONES.
        # CONSTANT-HEAD TERMS ARE USED TO IDENTIFY WHERE CONSTANT-HEAD CELLS ARE AND THEN USE
        # FACE FLOWS TO DETERMINE THE AMOUNT OF FLOW.
        # SWIADDTO--- terms are used by the SWI2 groundwater flow process.
        internal_flow_terms = ['CONSTANT HEAD', 'FLOW RIGHT FACE', 'FLOW FRONT FACE', 'FLOW LOWER FACE',
                               'SWIADDTOCH', 'SWIADDTOFRF', 'SWIADDTOFFF', 'SWIADDTOFLF']

        # Source/sink/storage term record names
        # These are all of the terms that are not related to constant
        # head cells or face flow terms
        self.ssst_record_names = [n for n in self.record_names
                                  if n not in internal_flow_terms]

        # Check the shape of the cbc budget file arrays
        self.cbc_shape = self.get_model_shape()
        self.nlay, self.nrow, self.ncol = self.cbc_shape

        self.float_type = np.float64
        return

    def get_model_shape(self):
        return self.cbc.get_data(idx=0, full3D=True)[0].shape

    def get_budget(self, z, **kwargs):
        """
        Creates a budget for the specified zone array. Pass keyword arguments to specify
        the time step/stress period or sim. time for which a budget is desired. Valid
        keywords are "kstpkper" and "totim". Currently, this function only supports the
        use of a single time step/stress period or time.

        :param z: A numpy.ndarray containing to zones to be used.
        :param kwargs:
        :return:
        Budget object
        """
        if 'kstpkper' in kwargs.keys():
            s = 'The specified time step/stress period' \
                ' does not exist {}'.format(kwargs['kstpkper'])
            assert kwargs['kstpkper'] in self.cbc.get_kstpkper(), print(s)
        elif 'totim' in kwargs.keys():
            s = 'The time ' \
                ' does not exist {}'.format(kwargs['totim'])
            assert kwargs['totim'] in self.cbc.get_times(), print(s)
        else:
            # FUTURE: Return a list of budgets for each time step/stress period
            # in the CellBudgetFile
            raise Exception('No stress period/time step or time specified.')
        assert isinstance(z, np.ndarray), 'Please pass zones as type {}'.format(np.ndarray)

        # Check for negative zone values
        negative_zones = [iz for iz in np.unique(z) if iz < 0]
        if len(negative_zones) > 0:
            raise Exception('Negative zone value(s) found:', negative_zones)

        # Make sure the input zone array has the same shape as the cell budget file
        if len(z.shape) == 2 and self.nlay == 1:
            izone = np.zeros(self.cbc_shape, np.int32)
            izone[0, :, :] = z[:, :]

        elif len(z.shape) == 2 and self.nlay > 1:
            raise Exception('Zone array and CellBudgetFile shapes do not match.')
        else:
            izone = z.copy()

        assert izone.shape == self.cbc_shape, \
            'Shape of input zone array {} does not' \
            ' match the cell by cell' \
            ' budget file {}'.format(izone.shape, self.cbc_shape)

        # List of unique zones numbers
        lstzon = [z for z in np.unique(izone)]

        # Initialize a constant head array
        ich = np.zeros(self.cbc_shape, np.int32)

        # Create empty array for the budget terms.
        # This array has the structure: ('flow direction', 'record name', value zone 1, value zone 2, etc.)
        self._initialize_records(lstzon)

        # Create a throwaway list of all record names
        reclist = list(self.record_names)

        if 'CONSTANT HEAD' in reclist:
            reclist.remove('CONSTANT HEAD')
            chd = self.cbc.get_data(text='CONSTANT HEAD', full3D=True, **kwargs)[0]
            ich = np.zeros(self.cbc_shape, np.int32)
            ich[chd != 0] = 1
        if 'FLOW RIGHT FACE' in reclist:
            reclist.remove('FLOW RIGHT FACE')
            self._accumulate_flow_frf('FLOW RIGHT FACE', izone, ich, **kwargs)
        if 'FLOW FRONT FACE' in reclist:
            reclist.remove('FLOW FRONT FACE')
            self._accumulate_flow_fff('FLOW FRONT FACE', izone, ich, **kwargs)
        if 'FLOW LOWER FACE' in reclist:
            reclist.remove('FLOW LOWER FACE')
            self._accumulate_flow_flf('FLOW LOWER FACE', izone, ich, **kwargs)
        if 'SWIADDTOCH' in reclist:
            reclist.remove('SWIADDTOCH')
            swichd = self.cbc.get_data(text='SWIADDTOCH', full3D=True, **kwargs)[0]
            swiich = np.zeros(self.cbc_shape, np.int32)
            swiich[swichd != 0] = 1
        if 'SWIADDTOFRF' in reclist:
            reclist.remove('SWIADDTOFRF')
            self._accumulate_flow_frf('SWIADDTOFRF', izone, swiich, **kwargs)
        if 'SWIADDTOFFF' in reclist:
            reclist.remove('SWIADDTOFFF')
            self._accumulate_flow_fff('SWIADDTOFFF', izone, swiich, **kwargs)
        if 'SWIADDTOFLF' in reclist:
            reclist.remove('SWIADDTOFLF')
            self._accumulate_flow_flf('SWIADDTOFLF', izone, swiich, **kwargs)

        # NOT AN INTERNAL FLOW TERM, SO MUST BE A SOURCE TERM OR STORAGE
        # ACCUMULATE THE FLOW BY ZONE
        # iterate over remaining items in the list
        for recname in reclist:

            imeth = self.imeth[recname]
            data = self.cbc.get_data(text=recname, **kwargs)[0]

            if imeth == 2 or imeth == 5:
                # LIST
                budin = np.ma.zeros((self.nlay * self.nrow * self.ncol), self.float_type)
                budout = np.ma.zeros((self.nlay * self.nrow * self.ncol), self.float_type)
                for [node, q] in zip(data['node'], data['q']):
                    idx = node - 1
                    if q > 0:
                        budin.data[idx] += q
                    elif q < 0:
                        budout.data[idx] += q
                budin = np.ma.reshape(budin, (self.nlay, self.nrow, self.ncol))
                budout = np.ma.reshape(budout, (self.nlay, self.nrow, self.ncol))
            elif imeth == 0 or imeth == 1:
                # FULL 3-D ARRAY
                budin = np.ma.zeros(self.cbc_shape, self.float_type)
                budout = np.ma.zeros(self.cbc_shape, self.float_type)
                budin[data > 0] = data[data > 0]
                budout[data < 0] = data[data < 0]
            elif imeth == 3:
                # 1-LAYER ARRAY WITH LAYER INDICATOR ARRAY
                rlay, rdata = data[0], data[1]
                data = np.ma.zeros(self.cbc_shape, self.float_type)
                for (r, c), l in np.ndenumerate(rlay):
                    data[l - 1, r, c] = rdata[r, c]
                budin = np.ma.zeros(self.cbc_shape, self.float_type)
                budout = np.ma.zeros(self.cbc_shape, self.float_type)
                budin[data > 0] = data[data > 0]
                budout[data < 0] = data[data < 0]
            elif imeth == 4:
                # 1-LAYER ARRAY THAT DEFINES LAYER 1
                budin = np.ma.zeros(self.cbc_shape, self.float_type)
                budout = np.ma.zeros(self.cbc_shape, self.float_type)
                r, c = np.where(data > 0)
                budin[0, r, c] = data[r, c]
                r, c = np.where(data < 0)
                budout[0, r, c] = data[r, c]
            else:
                # Should not happen
                raise Exception('Unrecognized "imeth" for {} record: {}'.format(recname, imeth))

            self._accumulate_flow_ssst(recname, budin, budout, izone, lstzon)

        return Budget(self.zonbudrecords, **kwargs)

    def _build_empty_record(self, flow_dir, recname, lstzon):
        recs = np.array(tuple([flow_dir, recname] + [0. for _ in lstzon if _ != 0]),
                        dtype=self.zonbudrecords.dtype)
        self.zonbudrecords = np.append(self.zonbudrecords, recs)

    def _initialize_records(self, lstzon):

        # Initialize the record array
        dtype_list = [('flow_dir', '|S3'), ('record', '|S20')]
        dtype_list += [('ZONE {:d}'.format(z), self.float_type) for z in lstzon if z != 0]
        dtype = np.dtype(dtype_list)
        self.zonbudrecords = np.array([], dtype=dtype)

        # Add "in" records
        if 'CONSTANT HEAD' in self.record_names:
            self._build_empty_record('in', 'CONSTANT HEAD', lstzon)
        for recname in self.ssst_record_names:
            self._build_empty_record('in', recname, lstzon)
        for z in lstzon:
            self._build_empty_record('in', 'FROM ZONE {}'.format(z), lstzon)

        # Add "out" records
        if 'CONSTANT HEAD' in self.record_names:
            self._build_empty_record('out', 'CONSTANT HEAD', lstzon)
        for recname in self.ssst_record_names:
            self._build_empty_record('out', recname, lstzon)
        for z in lstzon:
            self._build_empty_record('out', 'TO ZONE {}'.format(z), lstzon)

        return

    def _update_record(self, flow_dir, recname, colname, flux):
        if colname != 'ZONE 0':
            rowidx = np.where((self.zonbudrecords['flow_dir'] == flow_dir) &
                              (self.zonbudrecords['record'] == recname))
            self.zonbudrecords[colname][rowidx] += flux
        return

    def _accumulate_flow_frf(self, recname, izone, ich, **kwargs):
        # ACCUMULATE FLOW BETWEEN ZONES ACROSS COLUMNS. COMPUTE FLOW ONLY BETWEEN A ZONE
        # AND A HIGHER ZONE -- FLOW FROM ZONE 4 TO 3 IS THE NEGATIVE OF FLOW FROM 3 TO 4.
        # FIRST, CALCULATE FLOW BETWEEN NODE J,I,K AND J-1,I,K.
        # Accumulate flow from lower zones to higher zones from "left" to "right".
        # Flow into the higher zone will be <0 Flow Right Face from the adjacent cell to the "left".
        bud = self.cbc.get_data(text=recname, **kwargs)[0]

        nz = izone[:, :, 1:]
        nzl = izone[:, :, :-1]
        l, r, c = np.where(nz > nzl)

        # Adjust column values to account for the starting position of "nz"
        c = np.copy(c) + 1

        # Define the zone from which flow is coming
        from_zones = izone[l, r, c-1]

        # Define the zone to which flow is going
        to_zones = izone[l, r, c]

        # Get the face flow
        q = bud[l, r, c - 1]

        # Don't include CH to CH flow (can occur if CHTOCH option is used)
        q[(ich[l, r, c] == 1) & (ich[l, r, c-1] == 1)] = 0.

        # Get indices where flow face values are negative (flow into higher zone)
        idx_neg = np.where(q < 0)

        # Get indices where flow face values are positive (flow out of higher zone)
        idx_pos = np.where(q >= 0)

        # Create tuples of ("to zone", "from zone", "absolute flux")
        neg = zip(to_zones[idx_neg], from_zones[idx_neg], np.abs(q[idx_neg]))
        pos = zip(from_zones[idx_pos], to_zones[idx_pos], np.abs(q[idx_pos]))
        nzgt_l2r = neg + pos

        # CALCULATE FLOW TO CONSTANT-HEAD CELLS IN THIS DIRECTION
        l, r, c = np.where(ich == 1)

        # Can't accumulate left-to-right for cells on left edge of model (c = 0)
        l, r, c = l[c > 0], r[c > 0], c[c > 0]

        from_zones = izone[l, r, c-1]
        to_zones = izone[l, r, c]
        q = bud[l, r, c-1]
        q[(ich[l, r, c] == 1) & (ich[l, r, c-1] == 1)] = 0.
        idx_neg = np.where(q < 0)
        idx_pos = np.where(q >= 0)
        chdneg = zip(to_zones[idx_neg], from_zones[idx_neg], np.abs(q[idx_neg]))
        chdpos = zip(from_zones[idx_pos], to_zones[idx_pos], np.abs(q[idx_pos]))
        for (from_zone, to_zone, flux) in chdneg:
            self._update_record('in', 'CONSTANT HEAD', 'ZONE {}'.format(to_zone), flux)
        for (from_zone, to_zone, flux) in chdpos:
            self._update_record('out', 'CONSTANT HEAD', 'ZONE {}'.format(from_zone), flux)

        # CALCULATE FLOW BETWEEN NODE J,I,K AND J+1,I,K.
        # Accumulate flow from lower zones to higher zones from "right" to "left".
        # Flow into the higher zone will be <0 Flow Right Face from the adjacent cell to the "left".
        nz = izone[:, :, :-1]
        nzr = izone[:, :, 1:]
        l, r, c = np.where(nz > nzr)

        # Define the zone from which flow is coming
        from_zones = izone[l, r, c]

        # Define the zone to which flow is going
        to_zones = izone[l, r, c+1]

        # Get the face flow
        q = bud[l, r, c]

        # Don't include CH to CH flow (can occur if CHTOCH option is used)
        q[(ich[l, r, c] == 1) & (ich[l, r, c+1] == 1)] = 0.

        # Get indices where flow face values are negative (flow into higher zone)
        idx_neg = np.where(q < 0)

        # Get indices where flow face values are positive (flow out of higher zone)
        idx_pos = np.where(q >= 0)

        # Create tuples of ("to zone", "from zone", "absolute flux")
        neg = zip(to_zones[idx_neg], from_zones[idx_neg], np.abs(q[idx_neg]))
        pos = zip(from_zones[idx_pos], to_zones[idx_pos], np.abs(q[idx_pos]))
        nzgt_r2l = neg + pos

        # CALCULATE FLOW TO CONSTANT-HEAD CELLS IN THIS DIRECTION
        l, r, c = np.where(ich == 1)

        # Can't accumulate right-to-left for cells on right edge of model (c = ncol)
        l, r, c = l[c < self.ncol-1], r[c < self.ncol-1], c[c < self.ncol-1]

        from_zones = izone[l, r, c]
        to_zones = izone[l, r, c+1]
        q = bud[l, r, c]
        q[(ich[l, r, c] == 1) & (ich[l, r, c+1] == 1)] = 0.
        idx_neg = np.where(q < 0)
        idx_pos = np.where(q >= 0)
        chdneg = zip(to_zones[idx_neg], from_zones[idx_neg], np.abs(q[idx_neg]))
        chdpos = zip(from_zones[idx_pos], to_zones[idx_pos], np.abs(q[idx_pos]))
        for (from_zone, to_zone, flux) in chdneg:
            self._update_record('in', 'CONSTANT HEAD', 'ZONE {}'.format(to_zone), flux)
        for (from_zone, to_zone, flux) in chdpos:
            self._update_record('out', 'CONSTANT HEAD', 'ZONE {}'.format(from_zone), flux)

        # Update records
        nzgt = nzgt_l2r + nzgt_r2l
        for (from_zone, to_zone, flux) in nzgt:
            self._update_record('in', 'FROM ZONE {}'.format(from_zone), 'ZONE {}'.format(to_zone), flux)
            self._update_record('out', 'TO ZONE {}'.format(to_zone), 'ZONE {}'.format(from_zone), flux)
        return

    def _accumulate_flow_fff(self, recname, izone, ich, **kwargs):
        # ACCUMULATE FLOW BETWEEN ZONES ACROSS ROWS. COMPUTE FLOW ONLY BETWEEN A ZONE
        #  AND A HIGHER ZONE -- FLOW FROM ZONE 4 TO 3 IS THE NEGATIVE OF FLOW FROM 3 TO 4.
        # FIRST, CALCULATE FLOW BETWEEN NODE J,I,K AND J,I-1,K.
        # Accumulate flow from lower zones to higher zones from "up" to "down".
        # Returns a tuple of ("to zone", "from zone", "absolute flux")
        bud = self.cbc.get_data(text=recname, **kwargs)[0]

        nz = izone[:, 1:, :]
        nzu = izone[:, :-1, :]
        l, r, c = np.where(nz < nzu)
        # Adjust column values by +1 to account for the starting position of "nz"
        r = np.copy(r) + 1

        # Define the zone from which flow is coming
        from_zones = izone[l, r-1, c]

        # Define the zone to which flow is going
        to_zones = izone[l, r, c]

        # Get the face flow
        q = bud[l, r-1, c]

        # Don't include CH to CH flow (can occur if CHTOCH option is used)
        q[(ich[l, r, c] == 1) & (ich[l, r-1, c] == 1)] = 0.

        # Get indices where flow face values are negative (flow into higher zone)
        idx_neg = np.where(q < 0)

        # Get indices where flow face values are positive (flow out of higher zone)
        idx_pos = np.where(q >= 0)

        # Create tuples of ("to zone", "from zone", "absolute flux")
        neg = zip(to_zones[idx_neg], from_zones[idx_neg], np.abs(q[idx_neg]))
        pos = zip(from_zones[idx_pos], to_zones[idx_pos], np.abs(q[idx_pos]))
        nzgt_u2d = neg + pos

        # CALCULATE FLOW TO CONSTANT-HEAD CELLS IN THIS DIRECTION
        l, r, c = np.where(ich == 1)

        # Can't accumulate up-to-down for cells on top edge of model (r = 0)
        l, r, c = l[r > 0], r[r > 0], c[r > 0]

        from_zones = izone[l, r-1, c]
        to_zones = izone[l, r, c]
        q = bud[l, r-1, c]
        q[(ich[l, r, c] == 1) & (ich[l, r-1, c] == 1)] = 0.
        idx_neg = np.where(q < 0)
        idx_pos = np.where(q >= 0)
        chdneg = zip(to_zones[idx_neg], from_zones[idx_neg], np.abs(q[idx_neg]))
        chdpos = zip(from_zones[idx_pos], to_zones[idx_pos], np.abs(q[idx_pos]))
        for (from_zone, to_zone, flux) in chdneg:
            self._update_record('in', 'CONSTANT HEAD', 'ZONE {}'.format(to_zone), flux)
        for (from_zone, to_zone, flux) in chdpos:
            self._update_record('out', 'CONSTANT HEAD', 'ZONE {}'.format(from_zone), flux)

        # CALCULATE FLOW BETWEEN NODE J,I,K AND J,I+1,K.
        # Accumulate flow from lower zones to higher zones from "down" to "up".
        nz = izone[:, :-1, :]
        nzd = izone[:, 1:, :]
        l, r, c = np.where(nz < nzd)

        # Define the zone from which flow is coming
        from_zones = izone[l, r, c]

        # Define the zone to which flow is going
        to_zones = izone[l, r+1, c]

        # Get the face flow
        q = bud[l, r, c]

        # Don't include CH to CH flow (can occur if CHTOCH option is used)
        q[(ich[l, r, c] == 1) & (ich[l, r+1, c] == 1)] = 0.

        # Get indices where flow face values are negative (flow into higher zone)
        idx_neg = np.where(q < 0)

        # Get indices where flow face values are positive (flow out of higher zone)
        idx_pos = np.where(q >= 0)

        # Create tuples of ("to zone", "from zone", "absolute flux")
        neg = zip(to_zones[idx_neg], from_zones[idx_neg], np.abs(q[idx_neg]))
        pos = zip(from_zones[idx_pos], to_zones[idx_pos], np.abs(q[idx_pos]))
        nzgt_d2u = neg + pos

        # CALCULATE FLOW TO CONSTANT-HEAD CELLS IN THIS DIRECTION
        l, r, c = np.where(ich == 1)

        # Can't accumulate down-to-up for cells on bottom edge of model (r = nrow)
        l, r, c = l[r < self.nrow-1], r[r < self.nrow-1], c[r < self.nrow-1]

        from_zones = izone[l, r, c]
        to_zones = izone[l, r+1, c]
        q = bud[l, r, c]
        q[(ich[l, r, c] == 1) & (ich[l, r+1, c] == 1)] = 0.
        idx_neg = np.where(q < 0)
        idx_pos = np.where(q >= 0)
        chdneg = zip(to_zones[idx_neg], from_zones[idx_neg], np.abs(q[idx_neg]))
        chdpos = zip(from_zones[idx_pos], to_zones[idx_pos], np.abs(q[idx_pos]))
        for (from_zone, to_zone, flux) in chdneg:
            self._update_record('in', 'CONSTANT HEAD', 'ZONE {}'.format(to_zone), flux)
        for (from_zone, to_zone, flux) in chdpos:
            self._update_record('out', 'CONSTANT HEAD', 'ZONE {}'.format(from_zone), flux)

        # Update records
        nzgt = nzgt_u2d + nzgt_d2u
        for (from_zone, to_zone, flux) in nzgt:
            self._update_record('in', 'FROM ZONE {}'.format(from_zone), 'ZONE {}'.format(to_zone), flux)
            self._update_record('out', 'TO ZONE {}'.format(to_zone), 'ZONE {}'.format(from_zone), flux)
        return

    def _accumulate_flow_flf(self, recname, izone, ich, **kwargs):
        # ACCUMULATE FLOW BETWEEN ZONES ACROSS LAYERS. COMPUTE FLOW ONLY BETWEEN A ZONE
        #  AND A HIGHER ZONE -- FLOW FROM ZONE 4 TO 3 IS THE NEGATIVE OF FLOW FROM 3 TO 4.
        # FIRST, CALCULATE FLOW BETWEEN NODE J,I,K AND J,I,K-1.
        # Accumulate flow from lower zones to higher zones from "top" to "bottom".
        # Returns a tuple of ("to zone", "from zone", "absolute flux")
        bud = self.cbc.get_data(text=recname, **kwargs)[0]

        nz = izone[1:, :, :]
        nzt = izone[:-1, :, :]
        l, r, c = np.where(nz > nzt)
        # Adjust column values by +1 to account for the starting position of "nz"
        l = np.copy(l) + 1

        # Define the zone from which flow is coming
        from_zones = izone[l-1, r, c]

        # Define the zone to which flow is going
        to_zones = izone[l, r, c]

        # Get the face flow
        q = bud[l-1, r, c]

        # Don't include CH to CH flow (can occur if CHTOCH option is used)
        q[(ich[l, r, c] == 1) & (ich[l-1, r, c] == 1)] = 0.

        # Get indices where flow face values are negative (flow into higher zone)
        idx_neg = np.where(q < 0)

        # Get indices where flow face values are positive (flow out of higher zone)
        idx_pos = np.where(q >= 0)

        # Create tuples of ("to zone", "from zone", "absolute flux")
        neg = zip(to_zones[idx_neg], from_zones[idx_neg], np.abs(q[idx_neg]))
        pos = zip(from_zones[idx_pos], to_zones[idx_pos], np.abs(q[idx_pos]))
        nzgt_t2b = neg + pos

        # CALCULATE FLOW TO CONSTANT-HEAD CELLS IN THIS DIRECTION
        l, r, c = np.where(ich == 1)

        # Can't accumulate top-to-bottom for cells at top edge of model (l = 0)
        l, r, c = l[l > 0], r[l > 0], c[l > 0]

        from_zones = izone[l-1, r, c]
        to_zones = izone[l, r, c]
        q = bud[l-1, r, c]
        q[(ich[l, r, c] == 1) & (ich[l-1, r, c] == 1)] = 0.
        idx_neg = np.where(q < 0)
        idx_pos = np.where(q >= 0)
        chdneg = zip(to_zones[idx_neg], from_zones[idx_neg], np.abs(q[idx_neg]))
        chdpos = zip(from_zones[idx_pos], to_zones[idx_pos], np.abs(q[idx_pos]))
        for (from_zone, to_zone, flux) in chdneg:
            self._update_record('in', 'CONSTANT HEAD', 'ZONE {}'.format(to_zone), flux)
        for (from_zone, to_zone, flux) in chdpos:
            self._update_record('out', 'CONSTANT HEAD', 'ZONE {}'.format(from_zone), flux)

        # CALCULATE FLOW BETWEEN NODE J,I,K AND J+1,I,K.
        # Accumulate flow from lower zones to higher zones from "right" to "left".
        # Flow into the higher zone will be <0 Flow Right Face from the adjacent cell to the "left".
        nz = izone[:-1, :, :]
        nzb = izone[1:, :, :]
        l, r, c = np.where(nz < nzb)

        # Define the zone from which flow is coming
        from_zones = izone[l, r, c]

        # Define the zone to which flow is going
        to_zones = izone[l+1, r, c]

        # Get the face flow
        q = bud[l, r, c]

        # Don't include CH to CH flow (can occur if CHTOCH option is used)
        q[(ich[l, r, c] == 1) & (ich[l+1, r, c] == 1)] = 0.

        # Get indices where flow face values are negative (flow into higher zone)
        idx_neg = np.where(q < 0)

        # Get indices where flow face values are positive (flow out of higher zone)
        idx_pos = np.where(q >= 0)

        # Create tuples of ("to zone", "from zone", "absolute flux")
        neg = zip(to_zones[idx_neg], from_zones[idx_neg], np.abs(q[idx_neg]))
        pos = zip(from_zones[idx_pos], to_zones[idx_pos], np.abs(q[idx_pos]))
        nzgt_b2t = neg + pos

        # CALCULATE FLOW TO CONSTANT-HEAD CELLS IN THIS DIRECTION
        l, r, c = np.where(ich == 1)

        # Can't accumulate bottom-to-top for cells at bottom edge of model (l = nlay)
        l, r, c = l[l < self.nlay - 1], r[l < self.nlay - 1], c[l < self.nlay - 1]

        from_zones = izone[l, r, c]
        to_zones = izone[l+1, r, c]
        q = bud[l, r, c]
        q[(ich[l, r, c] == 1) & (ich[l+1, r, c] == 1)] = 0.
        idx_neg = np.where(q < 0)
        idx_pos = np.where(q >= 0)
        chdneg = zip(to_zones[idx_neg], from_zones[idx_neg], np.abs(q[idx_neg]))
        chdpos = zip(from_zones[idx_pos], to_zones[idx_pos], np.abs(q[idx_pos]))
        for (from_zone, to_zone, flux) in chdneg:
            self._update_record('in', 'CONSTANT HEAD', 'ZONE {}'.format(to_zone), flux)
        for (from_zone, to_zone, flux) in chdpos:
            self._update_record('out', 'CONSTANT HEAD', 'ZONE {}'.format(from_zone), flux)

        # Update records
        nzgt = nzgt_t2b + nzgt_b2t
        for (from_zone, to_zone, flux) in nzgt:
            self._update_record('in', 'FROM ZONE {}'.format(from_zone), 'ZONE {}'.format(to_zone), flux)
            self._update_record('out', 'TO ZONE {}'.format(to_zone), 'ZONE {}'.format(from_zone), flux)
        return

    def _accumulate_flow_ssst(self, recname, budin, budout, izone, lstzon):
        recin = [np.abs(budin[(izone == z)].sum()) for z in lstzon]
        recout = [np.abs(budout[(izone == z)].sum()) for z in lstzon]
        for idx, flux in enumerate(recin):
            if type(flux) == np.ma.core.MaskedConstant:
                flux = 0.
            self._update_record('in', recname, 'ZONE {}'.format(lstzon[idx]), flux)
        for idx, flux in enumerate(recout):
            if type(flux) == np.ma.core.MaskedConstant:
                flux = 0.
            self._update_record('out', recname, 'ZONE {}'.format(lstzon[idx]), flux)
        return

    def get_kstpkper(self):
        return self.cbc.get_kstpkper()

    def get_times(self):
        return self.cbc.get_times()

    def get_indices(self):
        return self.cbc.get_indices()


def _numpyvoid2numeric(a):
    return np.array([list(r) for r in a])
