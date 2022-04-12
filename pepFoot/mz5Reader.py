#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#    Copyright 2019 Jeddidiah Bellamy-Carter, University of Nottingham
#
#    This file is part of PepFoot.
#
#    PepFoot is free software: you can redistribute it and/or modify
#    it under the terms of the GNU Lesser General Public License as
#    published by the Free Software Foundation, either version 3 of
#    the License, or (at your option) any later version.
#
#    PepFoot is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU Lesser General Public License for more details.
#
#    You should have received a copy of the GNU Lesser General Public
#    License along with PepFoot.  If not, see <https://www.gnu.org/licenses/>.

import os
import numpy as np
from scipy import signal
import h5py

#  mz5 Keys:
# 'CVParam',
# 'CVReference',
# 'ChomatogramTime',
# 'ChromatogramIndex',
# 'ChromatogramIntensity',
# 'ChromatogramList',
# 'ChromatogramListBinaryData',
# 'ControlledVocabulary',
# 'DataProcessing',
# 'FileContent',
# 'FileInformation',
# 'InstrumentConfiguration',
# 'ParamGroups',
# 'RefParam',
# 'Run',
# 'Software',
# 'SourceFiles',
# 'SpectrumIndex',
# 'SpectrumIntensity',
# 'SpectrumListBinaryData',
# 'SpectrumMZ',
# 'SpectrumMetaData'


class mz5():
    """A class object for accessing data from an mz5 file"""

    def __init__(self, filename=None, in_memory=True):
        if not os.path.splitext(filename)[1] == '.mz5':
            raise ValueError('Incorrect file extension, must be .mz5')
        else:
            self.filename = filename
        self.file = h5py.File(self.filename, 'r')
        self.in_memory = in_memory
        self.num_scan = self.file['SpectrumIndex'].len()
        self.scan_lookup = np.zeros((self.num_scan),
                                    dtype=[('time', 'float32'),
                                           ('start', 'uint32'),
                                           ('end', 'uint32'),
                                           ('ms level', 'uint8'),
                                           ('precursor', 'float64'),
                                           ('min mz', 'uint32'),
                                           ('max mz', 'uint32')])
        self.time_ref = int(np.where(self.file['CVReference']['name'] ==
                                     b'scan start time')[0])
        self.msn_ref = int(np.where(self.file['CVReference']['name'] ==
                                    b'ms level')[0])
        self.min_mz_ref = int(np.where(self.file['CVReference']['name'] ==
                                       b'scan window lower limit')[0])
        self.max_mz_ref = int(np.where(self.file['CVReference']['name'] ==
                                       b'scan window upper limit')[0])

        try:
            self.precursor_ref = int(
                np.where(self.file['CVReference']['name'] == b'selected ion m/z')[0])
        except TypeError:
            print("No MSn spectra in this file.")

        self.fill_lookup()
        self.get_limits()
        if in_memory:
            self.mzs = self.get_all_mzs()

    def fill_lookup(self):
        """Extract ssential scan information from CVParam into
        self.scan_lookup"""

        self.scan_lookup['start'][1:] = self.file['SpectrumIndex'][:-1]
        self.scan_lookup['end'] = self.file['SpectrumIndex']
        meta = self.file["SpectrumMetaData"]["params"]["cvstart"]
        param_value = self.file['CVParam']['value'].astype('U')
        param_id = self.file['CVParam']['cvRefID']

        for scan in range(self.num_scan):
            start = meta[scan]
            if scan == self.num_scan-1:
                end = -1
            else:
                end = meta[scan+1]
            idx = np.where(param_id[start:end] == self.time_ref)[0]
            self.scan_lookup[scan]['time'] = param_value[start+idx][0]
            idx = np.where(param_id[start:end] == self.msn_ref)[0]
            self.scan_lookup[scan]['ms level'] = int(param_value[start+idx][0])
            if self.scan_lookup[scan]['ms level'] > 1:
                idx = np.where(param_id[start:end] == self.precursor_ref)[0]
                self.scan_lookup[scan]['precursor'] = param_value[start+idx][0]
            idx = np.where(param_id[start:end] == self.min_mz_ref)[0]
            if idx.size:
                self.scan_lookup[scan]['min mz'] = float(param_value[start+idx][0])
            else:
                print("Warning: scan #{} has no minimum m/z, it may be empty.".format(scan))
            idx = np.where(param_id[start:end] == self.max_mz_ref)[0]
            if idx.size:
                self.scan_lookup[scan]['max mz'] = float(param_value[start+idx][0])
            else:
                print("Warning: scan #{} has no maximum m/z, it may be empty.".format(scan))

    def get_limits(self):
        self.time_range = (self.scan_lookup['time'][0],
                           self.scan_lookup['time'][-1])
        scan_list = np.where(self.scan_lookup['ms level'] == 1)[0]
        self.ms1_range = (self.scan_lookup['min mz'][scan_list].min(),
                          self.scan_lookup['max mz'][scan_list].max())

    def get_mzs(self, scan):
        """Generate m/z array for a given scan"""

        if self.in_memory:
            return self.mzs[scan]
        else:
            return self.file['SpectrumMZ'][self.scan_lookup['start'][scan]:
                                       self.scan_lookup['end'][scan]].cumsum()

    def get_all_mzs(self):
        """Generate ragged array of all m/z arrays for data file.
        Only runs when in_memory == True
        """

        mzs = []
        for scan in range(self.num_scan):
            mzs.append(self.file['SpectrumMZ'][self.scan_lookup['start'][scan]:
                                       self.scan_lookup['end'][scan]].cumsum())

        return mzs

    def get_ints(self, scan, mz_idx=None):
        """Generate int array for a given scan and mz_idx if specified"""

        if mz_idx is not None:
            return self.file['SpectrumIntensity'][self.scan_lookup['start'][scan] + mz_idx[0]:
                                                  self.scan_lookup['start'][scan] + mz_idx[1]]
        else:
            return self.file['SpectrumIntensity'][self.scan_lookup['start'][scan]:
                                                  self.scan_lookup['end'][scan]]

    def get_scan_from_time(self, min_time, max_time):
        """Get scan numbers from a tuple of scan times"""

        return np.searchsorted(self.scan_lookup['time'], (min_time, max_time))

    def chromatogram(self, min_mz, max_mz, ms_level=1):
        """Generate 1xN array of ion current for m/z range by scan time.
           Returns an array of scan times and an array of intensities"""

        scan_list = np.where(self.scan_lookup['ms level'] == ms_level)[0]
        xic = np.zeros_like(scan_list, dtype='uint64')

        for i, scan in enumerate(scan_list):
            tmp_mz = self.get_mzs(scan)
            idx = np.where((tmp_mz >= min_mz) & (tmp_mz < max_mz))[0]
            if idx.any():
                xic[i] = self.get_ints(scan, (idx[0], idx[-1])).sum()

        return self.scan_lookup['time'][scan_list], xic

    def spectrum(self, min_time, max_time, ms_level=1, mz_range=None, precursor=None, tolerance=0.1):
        """Generate 1xN array of total ion current by m/z.
        ms_level allows selection for MS1 or MS2 """

        min_scan, max_scan = self.get_scan_from_time(min_time, max_time)

        if ms_level == 1:
            scan_list = np.where(self.scan_lookup['ms level'][min_scan:max_scan] == 1)[
                0] + min_scan
        elif ms_level == 2:
            if precursor is not None:
                scan_list = np.where((self.scan_lookup['ms level'][min_scan:max_scan] == 2) &
                                     (precursor - tolerance <=
                                      self.scan_lookup['precursor'][min_scan:max_scan]) &
                                     (precursor + tolerance >
                                      self.scan_lookup['precursor'][min_scan:max_scan]))[0] + min_scan
            else:
                scan_list = np.where(self.scan_lookup['ms level'][min_scan:max_scan] == 2)[0] + min_scan

        ref_scan = np.unique(self.get_mzs(scan_list[0]))
        bin_width = np.min(np.diff(ref_scan)) # Determines minimum spacing between m/z for interpolation
        ref_mz = np.arange(*self.ms1_range, bin_width)

        if mz_range:
            idx = np.where((ref_mz >= mz_range[0]) & (ref_mz < mz_range[1]))[0]
            ref_mz = ref_mz[idx]
        
        merge_int = np.zeros_like(ref_mz)

        for scan in scan_list:
            tmp_mz = self.get_mzs(scan)
            tmp_int = self.get_ints(scan)
            if mz_range:
                idx = np.where((tmp_mz >= mz_range[0]) & (tmp_mz < mz_range[1]))[0]
                tmp_mz = tmp_mz[idx]
                tmp_int = tmp_int[idx]
            if not tmp_mz.size:
                # If tmp_mz is empty after mz_range trimming, move onto next scan
                continue
            
            merge_int += np.interp(ref_mz, tmp_mz, tmp_int, left=0, right=0)

        return ref_mz, merge_int

    def get_area(self, rt_range, mz_range, ms_level=1):
        """Calculates the sum of intensities for given m/z and time ranges"""

        min_scan, max_scan = self.get_scan_from_time(*rt_range)
        min_mz, max_mz = mz_range
        scan_list = np.where(
            self.scan_lookup['ms level'][min_scan:max_scan] == ms_level)[0] + min_scan

        area = np.zeros(len(scan_list))
        for i, scan in enumerate(scan_list):
            tmp_mz = self.get_mzs(scan)
            idx = np.where((tmp_mz >= min_mz) & (tmp_mz < max_mz))[0]
            if idx.any():
                area[i] = self.get_ints(scan, (idx[0], idx[-1])).sum()

        return int(np.trapz(area, x=self.scan_lookup['time'][scan_list]*60))

    def detect_peaks(self, spectrum, threshold=500, distance=4, prominence=0.8, width=3):
        """Peak picking from a given spectrum using the relative maxima
        algorithm using a window of size order.

        Only peaks above the specified threshold are returned"""

        peaks, properties = signal.find_peaks(spectrum[1], height=threshold, prominence=prominence, width=width, distance=distance)

        return peaks
