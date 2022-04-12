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

import json
import os
import subprocess
import re
import sys
import traceback
import time
from datetime import datetime

import matplotlib.ticker as mtick
import numpy as np
from scipy import stats
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.widgets import SpanSelector
from matplotlib.backend_bases import MouseButton
from PyQt5 import QtCore as Qtc
from PyQt5 import QtGui as Qtg
from PyQt5 import QtWidgets as Qtw
from PyQt5.QtCore import Qt
from pyteomics import mass, parser

from pepFoot.gui.pepfoot import Ui_MainWindow
from pepFoot.gui.about import Ui_AboutDialog
from pepFoot.gui.preferences import Ui_Dialog

import pepFoot.mz5Reader as mz5Reader
from pepFoot.utility import *


def report_exception(exc_type, exc_value, exc_traceback):
      """Report unhandled exceptions and errors to user"""

      if issubclass(exc_type, KeyboardInterrupt):
        if Qtg.qApp:
          Qtg.qApp.quit()
        return

      Qtw.QMessageBox.warning(None, "Error", "An error has occurred. Please submit the following report to Jedd at j.s.g.bellamy-carter@bham.ac.uk\n\n{}".format("".join(traceback.format_exception(exc_type, exc_value, exc_traceback))))

      print("Closed due to an error. This is the full error report:")
      print("".join(traceback.format_exception(exc_type, exc_value, exc_traceback)))
      print("--End of report--")
      sys.exit(1)

os.environ['QTWEBENGINE_DISABLE_SANDBOX'] = '1'
sys.excepthook = report_exception

sys.argv.append('--disable-web-security')

if getattr(sys, 'frozen', False):
    BUNDLE_DIR = sys._MEIPASS
else:
    BUNDLE_DIR = os.path.dirname(os.path.abspath(__file__))

APP = Qtw.QApplication(sys.argv)
VERSION = '1.1.4b'
APP_ICON = Qtg.QIcon(os.path.join(BUNDLE_DIR, 'gui', 'icon.png'))
APP.setStyle("fusion")
APP.setAttribute(Qt.AA_EnableHighDpiScaling, True)
if hasattr(Qtw.QStyleFactory, 'AA_UseHighDpiPixmaps'):
    APP.setAttribute(Qt.AA_UseHighDpiPixmaps)

PROTON = mass.nist_mass['H+'][0][0]
NEUTRON = mass.nist_mass['C'][13][0] - mass.nist_mass['C'][12][0]

aa_mass = mass.std_aa_mass
aa_comp = mass.std_aa_comp
mods = {}
enzymes = {}

ORANGE = (0.9, 0.62, 0) #E69F00
SKYBLUE = (0.34, 0.71, 0.91) #56b4e9
GREEN = (0, 0.62, 0.45) #009E73
YELLOW = (0.94, 0.89, 0.26)
BLUE = (0, 0.45, 0.7) #0072B2
VERMILLION = (0.84, 0.37, 0) #d55e00
REDPURPLE = (0.8, 0.47, 0.65)
WHEAT = '#F5DEB3'
GREY = (0.9, 0.9, 0.9)#D3D3D3

class Peptide():
    """Represent peptide objects"""

    def __init__(self, sequence, id_, diff_mod, charges, ms1_range):
        self.sequence = sequence  # String with the unmodified peptide sequence
        self.id = id_
        self.mmass = mass.fast_mass2(self.sequence, aa_mass=aa_mass) # monoisotopic mass
        self.composition = mass.Composition(sequence=self.sequence, aa_comp=aa_comp)
        self.isotopes = get_isotopes(self.composition)
        self.mod_isotopes = get_isotopes(self.composition + aa_comp[mods[diff_mod]['ID']])
        self.mass = self.isotopes[0][self.isotopes[1].argmax()] # most abundant mass
        self.chgs = []  # List of probable charge states
        self.diff_mod = diff_mod  # Variable modification to peptide
        self.mzs = []  # m/z values for each charge state
        self.mod_mzs = []
        self.mz_ranges = [[], []]
        # Index 0 stores unmodified peptide data, index 1 stores modified peptide data
        self.rt_ranges = [[], []]
        self.chg = None
        self.areas = [0, 0]    #
        self.detected = -1  # -1, 0, 1, 2 Indicates which peptide variant is detected
        self.get_mass(charges, ms1_range)
        self.fragments = {}

    def get_mass(self, charges, ms1_range):
        """Calculates the m/z for valid ions, accounting for charge states
        where both the modified and unmodified peptide should be present"""
        self.chgs.clear()
        self.mzs.clear()
        for z in range(charges[0], charges[1]+1):
            mz = (self.mass)/z + PROTON  # Approximate for protonation
            if not self.diff_mod == 'None':
                # Approximate for protonation
                mod_mz = (self.mass + mods[self.diff_mod]['Mass'])/z + PROTON
            else:
                mod_mz = mz
            if (min([mz, mod_mz]) > ms1_range[0]) & (max([mz, mod_mz]) < ms1_range[1]):
                self.chgs.append(z)
                self.mzs.append(mz)
                self.mod_mzs.append(mod_mz)


class Main(Qtw.QMainWindow, Ui_MainWindow):
    def __init__(self,):
        super(Main, self).__init__()
        self.setupUi(self)
        self.showMaximized()
        self.settings = Qtc.QSettings(
            Qtc.QSettings.IniFormat, Qtc.QSettings.UserScope, 'pepFoot', 'pepFoot')
        self.settings.setFallbacksEnabled(False)
        self.setWindowTitle('pepFoot {}'.format(VERSION))
        self.load_settings()
        self.project = {'data files': [], 'name': ''}
        self.project_file = ''
        self.project_dir = ''
        self.peptides = []
        self.file_list = []
        self.length_range = [5, 20]
        self.charge_range = [1, 3]
        self._mz_ranges = [[],[]]
        self._rt_ranges = [[],[]]
        self.fix_mods = []
        self.diff_mod = 'None'
        self.coverExp = set()
        self.pdb_file = ''
        # Boolean to monitor states
        self.SAVED = True
        self.TREE_CHANGED = False
        self.PEPLIST_CHANGED = False
        self.LOAD_AREAS = False
        self.PDB_LOADED = False

        self.sequence = ''
        self.fill_mods()
        self.DiffMod.insertSeparator(1)
        self.Enzyme.insertSeparator(1)
        self.Enzyme.addItems(sorted(list(enzymes)))
        spacer = Qtw.QWidget(self)
        spacer.setSizePolicy(Qtw.QSizePolicy.Expanding, Qtw.QSizePolicy.Expanding)
        self.toolBar.insertWidget(self.actionHideParam, spacer)
        self.ProgressBar.reset()
        self.ProgressBar.setMaximum(1)

        # Setup triggers

        self.actionNewProject.triggered.connect(self.create_project)
        self.actionNew_from_Template.triggered.connect(self.create_from_project)
        self.actionSaveProject.triggered.connect(self.save_project)
        self.actionOpenProject.triggered.connect(self.open_project)
        self.actionExportCSV.triggered.connect(self.export_csv)
        self.actionQuit.triggered.connect(self.close)
        self.actionDocumentation.triggered.connect(self.open_documentation)
        self.actionAbout.triggered.connect(lambda: self.about_dialog(0))
        self.actionCite.triggered.connect(lambda: self.about_dialog(1))
        self.actionLicense.triggered.connect(lambda: self.about_dialog(2))
        self.actionEditMods.triggered.connect(
            lambda: self.preferences_dialog(0))
        self.actionEditCleavages.triggered.connect(
            lambda: self.preferences_dialog(1))
        self.UpdateBtn.clicked.connect(self.update_PepList)
        self.BatchBtn.clicked.connect(self.batch_process)
        self.actionBatch.triggered.connect(self.batch_process)
        self.PepList.itemClicked.connect(self.update_item)
        self.Ms1Charge.activated.connect(lambda: self.update_ms1())
        self.FileAdd.clicked.connect(self.file_add)
        self.FileRemove.clicked.connect(self.file_remove)
        self.FileConvert.clicked.connect(self.file_convert)
        self.FileConvert.setEnabled(True)
        self.RemoveAssignment.clicked.connect(self.remove_assignment)
        self.LoadPdb.clicked.connect(self.load_pdb)
        self.AnalysisUpdateBtn.clicked.connect(self.update_analysis)
        self.ZoomBtn.clicked.connect(self._ms1_zoom_active)
        self.ExtractBtn.clicked.connect(self._ms1_extract_active)
        self.IntegrateBtn.clicked.connect(self._ms1_area_active)
        self.NGLViewer.setHtml('<html><style>html {background-color: black;}</style></html>')
        self.NGL_continuous.stateChanged.connect(self.update_analysis)

        self.emblem = Qtg.QIcon()
        self.emblem.addPixmap(Qtg.QPixmap(":/icons/img/Gnome-emblem-default.svg"), Qtg.QIcon.Normal, Qtg.QIcon.Off)

        ######################
        # Initialising Plots #
        ######################

        self.ms1Fig = Figure(frameon=1)
        self.ms1Ax1 = self.ms1Fig.add_subplot(221)
        self.ms1Ax1.name = 'ms1Ax1'
        self.ms1Ax1.set_xlabel('Time')
        self.ms1Ax1.info = self.ms1Ax1.set_title('m/z: -', size='medium', loc='right')
        self.ms1Ax1.set_xlim(0, 10)
        self.ms1Ln1, = self.ms1Ax1.plot([], [], c=BLUE)
        self.ms1Lim1, = self.ms1Ax1.plot([], [], 'r', lw=1, marker=7, markersize=10, alpha=0.7)

        self.ms1Ax2 = self.ms1Fig.add_subplot(223)
        self.ms1Ax2.name = 'ms1Ax2'
        self.ms1Ax2.set_xlabel('m/z', style='italic')
        self.ms1Ax2.info = self.ms1Ax2.set_title('rt: -', size='medium', loc='right')
        self.ms1Ax2.set_xlim(100, 1000)
        self.ms1Ln2, = self.ms1Ax2.plot([], [], c=BLUE)
        self.ms1iso2, = self.ms1Ax2.plot([], [], c='hotpink', marker='o', alpha=0.7, transform=self.ms1Ax2.get_xaxis_transform())

        self.ms1Ax3 = self.ms1Fig.add_subplot(222)
        self.ms1Ax3.name = 'ms1Ax3'
        self.ms1Ax3.set_xlabel('Time')
        self.ms1Ax3.info = self.ms1Ax3.set_title('m/z: -', size='medium', loc='right')
        self.ms1Ax3.set_xlim(0, 10)
        self.ms1Ln3, = self.ms1Ax3.plot([], [], c=ORANGE)
        self.ms1Lim3, = self.ms1Ax3.plot([], [], 'r', lw=1, marker=7, markersize=10, alpha=0.7)

        self.ms1Ax4 = self.ms1Fig.add_subplot(224)
        self.ms1Ax4.name = 'ms1Ax4'
        self.ms1Ax4.set_xlabel('m/z', style='italic')
        self.ms1Ax4.info = self.ms1Ax4.set_title('rt: -', size='medium', loc='right')
        self.ms1Ax4.set_xlim(100, 1000)
        self.ms1Ln4, = self.ms1Ax4.plot([], [], c=ORANGE)
        self.ms1iso4, = self.ms1Ax4.plot([], [], c='hotpink', marker='o', alpha=0.7, transform=self.ms1Ax4.get_xaxis_transform())

        # Set common attributes
        for ax in self.ms1Fig.axes:
            ax.yaxis.set_major_formatter(
                mtick.ScalarFormatter(useMathText=True))
            ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
            ax.set_ylim(0, 10)

        self.ms1Fig.set_tight_layout(True)
        self.ms1Canvas = FigureCanvas(self.ms1Fig)
        self.Ms1PlotVL.addWidget(self.ms1Canvas)
        self.ms1Canvas.draw_idle()

        self.cid_up = self.ms1Canvas.mpl_connect(
            'button_press_event', self.zoomReset)
        self.zoom1 = SpanSelector(self.ms1Ax1, self.zoomAx1, 'horizontal', 0.01,
                                  useblit=True, rectprops=dict(alpha=0.2, facecolor='#333333'), button=1)
        self.zoom2 = SpanSelector(self.ms1Ax2, self.zoomAx2, 'horizontal', 0.01,
                                  useblit=True, rectprops=dict(alpha=0.2, facecolor='#333333'), button=1)
        self.zoom3 = SpanSelector(self.ms1Ax3, self.zoomAx3, 'horizontal', 0.01,
                                  useblit=True, rectprops=dict(alpha=0.2, facecolor='#333333'), button=1)
        self.zoom4 = SpanSelector(self.ms1Ax4, self.zoomAx4, 'horizontal',0.01,
                                  useblit=True, rectprops=dict(alpha=0.2, facecolor='#333333'), button=1)
        self.combine1 = SpanSelector(self.ms1Ax1, self.combineAx1, 'horizontal', 0.01,
                                     useblit=True, rectprops=dict(alpha=0.2, facecolor='#ef8a62'), button=1)
        self.combine3 = SpanSelector(self.ms1Ax3, self.combineAx3, 'horizontal', 0.01,
                                     useblit=True, rectprops=dict(alpha=0.2, facecolor='#ef8a62'), button=1)
        self.extract2 = SpanSelector(self.ms1Ax2, self.extractAx2, 'horizontal', 0.01,
                                     useblit=True, rectprops=dict(alpha=0.2, facecolor='#ef8a62'), button=1)
        self.extract4 = SpanSelector(self.ms1Ax4, self.extractAx4, 'horizontal', 0.01,
                                     useblit=True, rectprops=dict(alpha=0.2, facecolor='#ef8a62'), button=1)
        self.area1 = SpanSelector(self.ms1Ax1, self.areaAx1, 'horizontal',
                                       0.01, useblit=True, rectprops=dict(alpha=0.2, facecolor='#67a9cf'), button=1)
        self.area3 = SpanSelector(self.ms1Ax3, self.areaAx3, 'horizontal',
                                       0.01, useblit=True, rectprops=dict(alpha=0.2, facecolor='#67a9cf'), button=1)

        # self.Ms1PlotVL.addWidget(self.Ms1Toolbar)
        # self.toolbar.hide()
        self._ms1_zoom_active()
        self.barFig = Figure(frameon=1)
        barGS = self.barFig.add_gridspec(6,1)
        self.barAx1 = self.barFig.add_subplot(barGS[1:,0])
        self.barAx1.selected = []
        self.barAx1.tick_params(which='minor', length=0)
        self.barAx1.spines['top'].set_visible(False)
        self.barAx1.spines['right'].set_visible(False)
        self.barAx2 = self.barFig.add_subplot(barGS[0, 0])
        self.barAx2.selected = []
        self.barAx2.axis('off')
        self.barFig.set_tight_layout(True)
        self.barCanvas = FigureCanvas(self.barFig)
        self.BarPlotVL.addWidget(self.barCanvas)
        self.barCanvas.draw_idle()

        self.bar_pick = self.barCanvas.mpl_connect('button_press_event', self.on_bar_pick)

    def wait_effect(func):
        def wrapper(self, *args):
            APP.setOverrideCursor(Qtg.QCursor(Qt.BusyCursor))
            APP.processEvents()
            self.ProgressBar.setMaximum(0)
            self.ProgressBar.setValue(0)
            func(self, *args)
            self.ProgressBar.reset()
            self.ProgressBar.setMaximum(1)
            APP.restoreOverrideCursor()
        return wrapper

    def status(self, message, duration=0):
        """Update status information for user"""

        APP.processEvents()
        self.StatusBar.setText(message)
        if duration:
            APP.processEvents()
            time.sleep(duration/1000)
            self.StatusBar.setText("")

    def create_project(self):
        self.project_file, _ = Qtw.QFileDialog.getSaveFileName(self, 'Create pepFoot project', '',
                                                               'pepFoot Project (*.pfoot)')
        # Append pfoot extension if none exists
        if self.project_file:
            if not os.path.splitext(self.project_file)[1]:
                self.project_file += '.pfoot'
            self.project_dir = os.path.dirname(self.project_file)
            self.project = {'name': os.path.basename(self.project_file),
                            'creation date': str(datetime.today().strftime('%d %b %Y  %I:%M%p')),
                            'data files': [],
                            'sequence': '',
                            'charge range': [],
                            'length range': [],
                            'enzyme': '',
                            'missed cleave': 0,
                            'fixed mods': [],
                            'differential mod': '',
                            'peptides': [],
                            'charge array': [],
                            'm/z array': [[], []],
                            'rt array': [[], []],
                            'areas': [],
                            'fractional mod': [],
                            'pdb file': '',
                            'treatment': []}
            self.setWindowTitle(
                'pepFoot {} - {}'.format(VERSION, self.project['name']))
            # Reset Current Parameters
            self.FileList.clear()
            self.file_list[:] = []
            self.SeqInput.clear()
            self.PepList.clear()
            self.Tabs.setCurrentIndex(0)
            self.peptides.clear()
            self.reset_figures()
            self.ms1Canvas.draw_idle()
            self.barAx1.clear()
            self.barAx2.clear()
            self.barCanvas.draw_idle()
            self.actionSaveProject.setEnabled(True)
            self.save_project()
            self.PROJECT_OPEN = False

    def create_from_project(self):
        self.template_file, _ = Qtw.QFileDialog.getOpenFileName(
                self, 'Choose template pepFoot project', '', 'pepFoot Project (*.pfoot)')

        if self.template_file:
            self.project_file, _ = Qtw.QFileDialog.getSaveFileName(self, 'Create pepFoot project', '',
                                                           'pepFoot Project (*.pfoot)')
            if self.project_file:
                if not os.path.splitext(self.project_file)[1]:
                    self.project_file += '.pfoot'
                with open(self.template_file, 'r') as infile:
                    temp_project = json.load(infile)

                temp_project['name'] = os.path.basename(self.project_file)
                temp_project['creation date'] = str(datetime.today().strftime('%d %b %Y  %I:%M%p'))
                temp_project['data files'] = []
                temp_project['areas'] = []
                temp_project['fractional mod'] = []
                temp_project['pdb file'] = ''
                temp_project['treatment'] = []
                temp_project['pepfoot version'] = str(VERSION)

                with open(self.project_file, 'w') as outfile:
                    json.dump(temp_project, outfile, indent=1, sort_keys=True)

                self.open_project(self.project_file)

    def open_project(self, filename=None):
        if filename:
            if os.path.splitext(filename)[1] == '.pfoot':
                self.project_file = os.path.abspath(filename)
            else:
                return
        else:
            self.project_file, _ = Qtw.QFileDialog.getOpenFileName(
                self, 'Open pepFoot project', '', 'pepFoot Project (*.pfoot)')

        if self.project_file:
            self.settings.setValue('Last File', self.project_file)
            self.project_dir = os.path.dirname(self.project_file)
            self.status('Opening {} ...'.format(os.path.basename(self.project_file)), 500)

            with open(self.project_file, 'r') as infile:
                self.project = json.load(infile)

            try:
                self.actionSaveProject.setEnabled(True)
                self.setWindowTitle('pepFoot {} - {}'.format(VERSION, self.project['name']))
                self.PepList.clear()
                self.peptides.clear()
                self.reset_figures()
                self.ms1Canvas.draw_idle()
                self.barAx1.clear()
                self.barAx2.clear()
                self.barCanvas.draw_idle()
                self.PDB_LOADED = False
                self.sequence = self.project['sequence']
                self.SeqInput.setPlainText(self.sequence)
                self.Enzyme.setCurrentText(self.project['enzyme'])
                self.MissCleave.setValue(self.project['missed cleave'])
                self.length_range = self.project['length range']
                self.MinLen.setValue(self.length_range[0])
                self.MaxLen.setValue(self.length_range[1])
                self.charge_range = self.project['charge range']
                self.MinChg.setValue(self.charge_range[0])
                self.MaxChg.setValue(self.charge_range[1])
                self.DiffMod.setCurrentText(self.project['differential mod'])

                if self.project['pdb file']:
                    if not os.path.isabs(self.project['pdb file']):
                        self.pdb_file = os.path.normpath(os.path.join(self.project_dir, self.project['pdb file']))
                    else:
                        self.pdb_file = self.project['pdb file']

                for mod in self.project['fixed mods']:
                    item = self.FixModList.findItems(mod, Qt.MatchExactly)[0]
                    self.FixModList.setCurrentItem(item)
                self.FileList.clear()
                self.file_list[:] = []

                if self.project['data files']:
                    for name in self.project['data files']:
                        if name not in self.file_list:
                            self.file_list.append(name)
                            self.FileList.addItem(os.path.basename(name))

                self.ApoList.clear()
                self.HoloList.clear()

                if self.project['treatment']:
                    for i in self.project['treatment'][0]:
                        self.ApoList.addItem(os.path.basename(self.file_list[i]))
                    for i in self.project['treatment'][1]:
                        self.HoloList.addItem(os.path.basename(self.file_list[i]))
                else:
                    for file_ in self.file_list:
                        self.ApoList.addItem(os.path.basename(file_))

                if self.project['areas'][0][0]:
                    self.LOAD_AREAS = True

                    if all(len(x) == len(self.project['fractional mod'][0]) for x in self.project['fractional mod']):
                        self.update_analysis()

                        if self.pdb_file:
                            self.load_pdb(self.pdb_file)
                    else:
                        self.status('Run Batch Process!')

            except Exception as e:
                self.status(str(e), 2000)

    def save_project(self):
        """Saves project files"""

        if not self.project['name'] or self.project['name'] == '.pfoot':
            self.project_file, _ = Qtw.QFileDialog.getSaveFileName(self, 'Save pepFoot project', '',
                                                               'pepFoot Project (*.pfoot)')
            if not self.project_file:
                return

            if not os.path.splitext(self.project_file)[1]:
                self.project_file += '.pfoot'

            self.project_dir = os.path.dirname(self.project_file)
            self.project = {'name': os.path.basename(self.project_file),
                'creation date': str(datetime.today().strftime('%d %b %Y  %I:%M%p')),
                'data files': [],
                'sequence': '',
                'length range': [],
                'charge range': [],
                'enzyme': '',
                'missed cleave': 0,
                'fixed mods': [],
                'differential mod': '',
                'peptides': [],
                'charge array': [],
                'm/z array': [[], []],
                'rt array': [[], []],
                'areas': [],
                'fractional mod': [],
                'pdb file': '',
                'treatment': [],
                'pepfoot version': str(VERSION)}
        else:
            if not self.project_file:
                return

        self.status('Saving {} ...'.format(
            os.path.basename(self.project['name'])))

        if self.TREE_CHANGED:
            self.project['data files'].clear()
            self.project['areas'].clear()
            self.project['fractional mod'].clear()
            self.project['treatment'] = []
            for i, file in enumerate(self.file_list):
                if os.path.isabs(file):
                    _path = file
                else:
                    _path = os.path.normpath(os.path.join(self.project_dir, file))
                try:
                    self.project['data files'].append(os.path.relpath(_path, start=self.project_dir))
                    self.project['areas'].append([[], []])
                    self.project['fractional mod'].append([])
                except:
                    self.status('Could not resolve paths between project and {}'.format(file))
            self.TREE_CHANGED = False

        if not self.project['data files']:
            self.status('No data files present in project, could not save!')
            return

        if self.pdb_file:
            self.project['pdb file'] = os.path.relpath(self.pdb_file, start=self.project_dir)

        if self.PEPLIST_CHANGED:
            self.project['sequence'] = re.sub('[^ACDEFGHIKLMNPQRSTVWYabcdefghijklmnopqrstuvwxyz]', '', self.SeqInput.toPlainText())
            self.project['enzyme'] = str(self.Enzyme.currentText())
            self.project['missed cleave'] = self.MissCleave.value()
            self.project['differential mod'] = self.diff_mod
            self.project['fixed mods'] = self.fix_mods
            self.project['length range'] = self.length_range
            self.project['charge range'] = self.charge_range
            self.project['peptides'].clear()
            self.project['charge array'].clear()
            self.project['m/z array'][0].clear()
            self.project['m/z array'][1].clear()
            self.project['rt array'][0].clear()
            self.project['rt array'][1].clear()
            self.project['areas'].clear()
            self.project['fractional mod'].clear()
            for i in range(len(self.file_list)):
                self.project['areas'].append([[], []])
                self.project['fractional mod'].append([])
            for pep in self.peptides:
                if not pep.detected == -1:
                    self.project['peptides'].append(pep.id)
                    self.project['charge array'].append(pep.chg)
                    self.project['m/z array'][0].append(pep.mz_ranges[0])
                    self.project['rt array'][0].append(pep.rt_ranges[0])
                    self.project['m/z array'][1].append(pep.mz_ranges[1])
                    self.project['rt array'][1].append(pep.rt_ranges[1])
                    self.project['areas'][0][0].append(pep.areas[0])
                    self.project['areas'][0][1].append(pep.areas[1])
                    self.project['fractional mod'][0].append(
                        pep.areas[1]/(pep.areas[0]+pep.areas[1]))

        with open(self.project_file, 'w') as outfile:
            json.dump(self.project, outfile, indent=1, sort_keys=True)

        self.status('{} saved!'.format(os.path.basename(self.project['name'])), 500)
        self.SAVED = True

    def export_csv(self):
        """Export data from current .pfoot file into the .csv format.
        Lines 1-6 are comments and lines >7 are for data with the form:

        1. Generated by pepfoot #.#.#
        2. From file ####.pfoot
        3. Date and Time this CSV is created
        4. Data files : [###,###,.....]
        5. Treatment: [[#,...], [#,...]]
        6. Differential modification, mod mass

        id, sequence, charge, m/z, mod m/z, rt, mod rt, avg. f.mod, std. dev. f.mod,
            (avg. f.mod 2, std. dev. f.mod 2, p value), [areas by file], [f.mod by file]

        """

        file_name = self.project_file.replace('.pfoot', '.csv')
        fmod = np.array(self.project['fractional mod'])
        num_files = len(self.project['data files'])
        if self.project['treatment']:
            treatments = self.project['treatment']
            p_values = stats.ttest_ind(fmod[treatments[0]], fmod[treatments[1]], axis=0)[1]
        else:
            treatments = None

        if os.path.exists(file_name) and not os.access(file_name, os.W_OK):
            self.status('Could not write to {}, open in another application'.format(os.path.basename(file_name)),2000)
            return

        with open(file_name, 'w') as outfile:
            outfile.write('Generated by pepFoot {}\n'.format(VERSION))
            outfile.write('{}\n'.format(self.project['name']))
            outfile.write('{}\n'.format(datetime.today().strftime('%d %b %Y  %I:%M%p')))
            outfile.write(','.join(['File {}: {}'.format(i, os.path.basename(name)) for i, name in enumerate(self.project['data files'])]))
            outfile.write('\n{}'.format(treatments))
            outfile.write('\n{}\n\n'.format(self.project['differential mod']))

            header = 'Peptide ID,Sequence,Charge,m/z,Mod m/z,rt,Mod rt,'
            if treatments:
                header += 'Mean F.mod 1,S.dev. F.mod 1,Mean F.mod 2,S.dev. F.mod 2,p value, Extent Mean, Extent S.dev,'
            else:
                header += 'Mean F.mod,S.dev F.mod,'
            for i in range(num_files):
                header += 'File {0} PA,File {0} Mod PA,File {0} F.mod,'.format(i)

            outfile.write(header[:-1] + '\n')

            for i, pep in enumerate(self.project['peptides']):
                line = '"{} - {}",'.format(*pep)
                line += '{},'.format(''.join(parser.parse(self.project['sequence'])[pep[0]-1:pep[1]]))
                line += '{},'.format(self.project['charge array'][i])
                if self.project['m/z array'][0][i] and not self.project['m/z array'][1][i]:
                    line += '{:.4f} - {:.4f},'.format(*self.project['m/z array'][0][i])
                    line += ','
                    line += '{:.2f} - {:.2f},'.format(*self.project['rt array'][0][i])
                    line += ','
                elif not self.project['m/z array'][0][i] and self.project['m/z array'][1][i]:
                    line += ','
                    line += '{:.4f} - {:.4f},'.format(*self.project['m/z array'][1][i])
                    line += ','
                    line += '{:.2f} - {:.2f},'.format(*self.project['rt array'][1][i])
                else:
                    line += '{:.4f} - {:.4f},'.format(*self.project['m/z array'][0][i])
                    line += '{:.4f} - {:.4f},'.format(*self.project['m/z array'][1][i])
                    line += '{:.2f} - {:.2f},'.format(*self.project['rt array'][0][i])
                    line += '{:.2f} - {:.2f},'.format(*self.project['rt array'][1][i])

                if treatments:
                    line += '{:.4f},'.format(fmod[treatments[0]][:,i].mean())
                    line += '{:.4f},'.format(fmod[treatments[0]][:,i].std())
                    line += '{:.4f},'.format(fmod[treatments[1]][:,i].mean())
                    line += '{:.4f},'.format(fmod[treatments[1]][:,i].std())
                    line += '{:.2g},'.format(p_values[i])
                    line += '{:.3f}, {:.3f},'.format(*calc_extent_change(fmod[treatments[0]][:,i].mean(), fmod[treatments[0]][:,i].std(), fmod[treatments[1]][:,i].mean(), fmod[treatments[1]][:,i].std()))

                else:
                    line += '{:.4f},'.format(fmod[:,i].mean())
                    line += '{:.4f},'.format(fmod[:,i].std())

                for j in range(num_files):
                    line += '{},'.format(self.project['areas'][j][0][i])
                    line += '{},'.format(self.project['areas'][j][1][i])
                    line += '{:.4f},'.format(fmod[j][i])
                outfile.write(line[:-1] + '\n')

        self.status('Successfully exported to {}!'.format(os.path.basename(file_name), 500))

    def batch_process(self):
        try:
            self.save_project()
            self.PEPLIST_CHANGED = False
            self.status('Batch processing {} ...'.format(self.project['name']))
            total = len(self.project['data files'])
            count = 0
            APP.processEvents()
            self.ProgressBar.setMaximum(total)
            self.ProgressBar.setValue(0)
            num_peptides = len(self.project['peptides'])
            for i, data_file in enumerate(self.project['data files']):
                self.open_mz5(data_file)
                count += 1
                APP.processEvents()
                self.ProgressBar.setValue(count)
                self.project['areas'][i].clear()
                self.project['fractional mod'][i].clear()
                self.project['areas'][i] = [[], []]

                self.status('Batch integrating {} ...'.format(
                    os.path.basename(data_file)))
                for j in range(num_peptides):
                    unmod = 0
                    mod = 0
                    if self.project['rt array'][0][j]:
                        unmod = self.data.get_area(
                            self.project['rt array'][0][j], self.project['m/z array'][0][j])
                    self.project['areas'][i][0].append(unmod)
                    if self.project['rt array'][1][j]:
                        mod = self.data.get_area(
                            self.project['rt array'][1][j], self.project['m/z array'][1][j])
                    self.project['areas'][i][1].append(mod)
                    if unmod and mod:
                        self.project['fractional mod'][i].append(mod/(unmod+mod))
                    elif unmod and not mod:
                        self.project['fractional mod'][i].append(0)
                    elif mod and not unmod:
                        self.project['fractional mod'][i].append(1)
                    else:
                        self.project['fractional mod'][i].append(0)
            APP.processEvents()
            self.ProgressBar.reset()
            self.ProgressBar.setMaximum(1)
            self.save_project()
            self.open_project(self.project_file)
            #self.update_analysis()
            self.status('Batch processing of {} complete!'.format(
                self.project['name']), 2000)
        except AssertionError as error:
            self.status(error.args[0])

    def file_convert(self):
        """Convert RAW files to .mz5 with ProteoWizard msconvert"""
        if self.pwiz:
              if not os.path.exists(self.pwiz):
                    self.pwiz = None
        if not self.pwiz:
            pwiz_loc = Qtw.QFileDialog.getOpenFileName(self, 'Select msconvert executable', '', 'All Files(*)')[0]
            if pwiz_loc:
                self.pwiz = pwiz_loc
            else:
                return

        file_names = Qtw.QFileDialog.getOpenFileNames(self, 'Select raw file(s) to convert', '',
                                                      'Vendor MS Files (*.RAW *.raw *.d *.YEP *.FID *.BAF *.WIFF)')[0]
        if file_names:
            for name in file_names:
                try:
                    subprocess.run([self.pwiz, name, "--zlib", "--filter", "\"zeroSamples removeExtra\"", "--mz5", "-o", os.path.dirname(name)], check=True)
                    self.TREE_CHANGED = True
                    name_ = os.path.splitext(name)[0] + '.mz5'
                    if name_ not in self.file_list:
                        self.file_list.append(name_)
                        self.FileList.addItem(os.path.basename(name_))

                except subprocess.CalledProcessError:
                    self.status('Unable to convert {}'.format(name))
                    continue

    def file_add(self):
        self.TREE_CHANGED = True
        file_names = Qtw.QFileDialog.getOpenFileNames(
            self, 'Select mz5 file(s)', '', 'mz5 Files (*.mz5)')[0]
        if file_names:
            for name in file_names:
                if name not in self.file_list:
                    self.file_list.append(name)
                    self.FileList.addItem(os.path.basename(name))

    def file_remove(self):
        if not self.FileList.currentRow() == -1:
            item = self.FileList.currentRow()
            del self.file_list[item]
            self.FileList.takeItem(item)
            self.TREE_CHANGED = True

    @wait_effect
    def file_import(self):
        assert self.file_list, 'No file selected, please select a data file.'
        if self.FileList.currentRow() == -1:
            file_name = self.file_list[0]
            self.FileList.setCurrentRow(0)
        else:
            file_name = self.file_list[self.FileList.currentRow()]
        self.open_mz5(file_name)

    def open_mz5(self, file_name):
        if not os.path.isabs(file_name):
            file_name = os.path.normpath(os.path.join(self.project_dir, file_name))
        if not os.path.exists(file_name):
            raise AssertionError('Cannot find file: {}!'.format(file_name))
        self.status('Reading {}...'.format(os.path.basename(file_name)))
        try:
            self.data = mz5Reader.mz5(file_name, self.actionInMemory.isChecked())
        except ValueError:
            self.status('Incorrect file extension, must be .mz5')
        self.status('Successfully imported {}'.format(os.path.basename(file_name)))

    def fill_mods(self):
        self.FixModList.clear()
        for i, e in enumerate(sorted(list(mods))):
            self.FixModList.addItem(e)
            self.FixModList.item(i).setToolTip(
                '{}: {}'.format(mods[e]['ID'], mods[e]['Residues']))
            self.DiffMod.addItem(e)
            aa_mass[mods[e]['ID']] = mods[e]['Mass']
            aa_comp[mods[e]['ID']] = mass.Composition(formula=mods[e]['Gain']) - mass.Composition(formula=mods[e]['Loss'])

    def digest(self, sequence, enzyme, miss_cleave):
        """Generates a list of valid peptides using pyteomics::parser.cleave

        Tolerant to modX input sequences and will only produce peptides
        capable of being differentially modified according to diff_mod['Residues']
        """
        self.peptides.clear()
        _sequence = re.sub('[^ACDEFGHIKLMNPQRSTVWY]', '', sequence)
        diff_mod_res = set(mods[self.diff_mod]['Residues'])
        sequence_list = parser.parse(sequence)
        for _peptide in parser.cleave(_sequence, enzymes[enzyme], miss_cleave):
            try:
                length = len(_peptide)
                pos = _sequence.index(_peptide)
                peptide = ''.join(sequence_list[pos:pos+length])
                if ((length > self.length_range[0]) & (length < self.length_range[1])):
                    if (diff_mod_res & set(peptide)):
                        self.peptides.append(Peptide(peptide, (pos+1, pos+length), self.diff_mod, self.charge_range, self.data.ms1_range))
            except ValueError:
                pass
        self.peptides.sort(key=lambda x: (x.id[0], x.sequence))

    def update_mods(self):
        self.fix_mods = [item.text()
                         for item in self.FixModList.selectedItems()]
        self.diff_mod = self.DiffMod.currentText()
        if self.diff_mod in self.fix_mods:
            self.status(
                '{} treated as differential modification.'.format(self.diff_mod))
            self.fix_mods.remove(self.diff_mod)

    def update_PepList(self):

        # Give user warning if PepList is not empty
        if self.PepList.count():
            reply = Qtw.QMessageBox.question(self, 'Overwrite Warning',
                                             'This will overwrite the current list of peptides! Do you wish to proceed?',
                                             Qtw.QMessageBox.Yes | Qtw.QMessageBox.Cancel, Qtw.QMessageBox.Cancel)
            if reply == Qtw.QMessageBox.Cancel:
                return

        if self.DiffMod.currentText() == 'None':
            self.status('No differential modification selected, cannot proceed!')
            return

        try:
            self.PepList.clear()
            self.PEPLIST_CHANGED = True
            self.file_import()
            self.length_range = [self.MinLen.value(), self.MaxLen.value()]
            self.charge_range = [self.MinChg.value(), self.MaxChg.value()]
            self.sequence = re.sub('[^ACDEFGHIKLMNPQRSTVWYabcdefghijklmnopqrstuvwxyz]', '', self.SeqInput.toPlainText())

            covered = set()
            assert self.sequence, 'Please enter a peptide sequence'
            self.update_mods()
            mod_dict = {mods[m]['ID']: mods[m]['Residues']
                        for m in self.fix_mods}
            if self.fix_mods:
                self.sequence = next(parser.isoforms(self.sequence, fixed_mods=mod_dict))
            if self.Enzyme.currentText() == 'None':
                self.peptides.clear()
                self.peptides.append(Peptide(
                    self.sequence, 1, self.diff_mod, self.charge_range, self.data.ms1_range))
            else:
                self.status('Digesting sequence...')
                self.digest(self.sequence, str(self.Enzyme.currentText()), self.MissCleave.value())
                self.status('Sequence digested!')

            for pep in self.peptides:
                self.PepList.addItem(pep.sequence)
            for i, pep in enumerate(self.peptides):
                if not pep.mzs:
                    self.PepList.item(i).setHidden(True)
                else:
                    self.PepList.item(i).setHidden(False)
                    covered.update(range(pep.id[0], pep.id[1]+1))

            self.PepList.repaint()
            assert len(covered), 'No peptides for current parameters!'
            self.CoverPred.setText('{:.0f} %'.format(
                len(covered)*100/len(re.sub('[^ACDEFGHIKLMNPQRSTVWY]', '', self.sequence))))
            self.CoverExp.setText('N/A')
            self.coverExp = set()
            self.Tabs.setCurrentIndex(0)
            self.SAVED = False

            if self.LOAD_AREAS:
                self.fill_peptides(self.project['peptides'], self.project['charge array'],
                                   self.project['m/z array'], self.project['rt array'],
                                   self.project['areas'])
                self.LOAD_AREAS = False

        except AssertionError as error:
            self.status(error.args[0])

    def fill_peptides(self, ids, charges, mz_ranges, rt_ranges, areas):
        """Fill self.peptides with data imported from file."""

        file_ind = self.file_list.index(os.path.relpath(self.data.filename, start=self.project_dir))

        for i, pep in enumerate(self.peptides):
            if list(pep.id) in ids:
                idx = ids.index(list(pep.id))
                pep.mz_ranges = [mz_ranges[0][idx], mz_ranges[1][idx]]
                pep.rt_ranges = [rt_ranges[0][idx], rt_ranges[1][idx]]
                pep.areas = [areas[file_ind][0][idx], areas[file_ind][1][idx]]
                self.coverExp.update(range(pep.id[0], pep.id[1]+1))
                if pep.areas[0] and pep.areas[1]:
                    pep.detected = 2
                elif pep.areas[0]:
                    pep.detected = 0
                elif pep.areas[1]:
                    pep.detected = 1
                pep.chg = charges[idx]
                self.PepList.item(i).setIcon(self.emblem)

        self.CoverExp.setText('{:.0f} %'.format(len(self.coverExp)*100/len(re.sub('[^ACDEFGHIKLMNPQRSTVWY]', '', self.sequence))))

    def update_item(self):
        ind = self.PepList.currentRow()
        self.RemoveAssignment.setEnabled(True)
        self.Ms1Charge.clear()
        for each in self.peptides[ind].chgs:
            self.Ms1Charge.addItem(str(each))
        if self.peptides[ind].chg:
            idx = self.peptides[ind].chgs.index(self.peptides[ind].chg)
            if idx:
                self.Ms1Charge.setCurrentIndex(idx)
        self.update_ms1()

    @wait_effect
    def update_ms1(self):
        """Extract chromatograms for unlabelled and labelled peptide ions.
           If data are already saved to self.peptides then spectra will be
           combined.
        """
        pepind = self.PepList.currentRow()
        chgind = self.Ms1Charge.currentIndex()
        _mz = self.peptides[pepind].mzs[chgind]
        self.PepMzLbl.setText('{:.2f}'.format(_mz))
        _mod_mz = self.peptides[pepind].mod_mzs[chgind]
        self.PepModMzLbl.setText('{:.2f}'.format(_mod_mz))
        self.PepMzErr.setText('N/A')
        self.PepModMzErr.setText('N/A')

        if self.TolUnit.currentText() == 'mmu':
            _tolerance = self.TolVal.value() / 1000
            _mod_tolerance = _tolerance
        elif self.TolUnit.currentText() == 'ppm':
            _tolerance = self.TolVal.value() * _mz / 1000000
            _mod_tolerance = self.TolVal.value() * _mod_mz / 1000000
        self._mz_ranges[0] = (_mz - _tolerance, _mz + _tolerance)
        self._mz_ranges[1] = (_mod_mz - _mod_tolerance, _mod_mz + _mod_tolerance)

        if self.peptides[pepind].chg == self.peptides[pepind].chgs[chgind]:
            if self.peptides[pepind].mz_ranges[0]:
                self._mz_ranges[0] = self.peptides[pepind].mz_ranges[0]
            if self.peptides[pepind].mz_ranges[1]:
                self._mz_ranges[1] = self.peptides[pepind].mz_ranges[1]

        self.status('Extracting chromatogram...')
        self.reset_figures()
        self.ms1Ln1.set_data(self.data.chromatogram(*self._mz_ranges[0]))
        self.ms1Ln3.set_data(self.data.chromatogram(*self._mz_ranges[1]))
        self._autoscale_y(self.ms1Ax1, *self.data.time_range)
        self._autoscale_y(self.ms1Ax3, *self.data.time_range)
        self.ms1Ax1.info.set_text('m/z: {:.2f}-{:.2f}'.format(*self._mz_ranges[0]))
        self.ms1Ax2.info.set_text('rt: -')
        self.ms1Ax3.info.set_text('m/z: {:.2f}-{:.2f}'.format(*self._mz_ranges[1]))
        self.ms1Ax4.info.set_text('rt: -')
        self.ms1Canvas.draw_idle()
        self.PepIdLbl.setText('{}-{}'.format(*self.peptides[pepind].id))

        if self.peptides[pepind].chg == self.peptides[pepind].chgs[chgind]:
            if self.peptides[pepind].areas[0]:
                self.UnmodArea.setText(str(self.peptides[pepind].areas[0]))
                self.combineAx1(*self.peptides[pepind].rt_ranges[0])
            else:
                self.UnmodArea.setText('N/A')
            if self.peptides[pepind].areas[1]:
                self.ModArea.setText(str(self.peptides[pepind].areas[1]))
                self.combineAx3(*self.peptides[pepind].rt_ranges[1])
            else:
                self.ModArea.setText('N/A')
            if (self.peptides[pepind].areas[0] & self.peptides[pepind].areas[0]):
                self.FMod.setText('{:.2f}'.format(self.peptides[pepind].areas[1]/(
                    self.peptides[pepind].areas[1]+self.peptides[pepind].areas[0])))
            else:
                self.FMod.setText('N/A')
        else:
            self.UnmodArea.setText('N/A')
            self.ModArea.setText('N/A')
            self.FMod.setText('N/A')

        self.status('Chromatogram extracted!')

    def _autoscale_y(self, ax, min_, max_):

        ax.set_xlim(min_, max_)
        lines = ax.get_lines()
        top = 0

        for line in lines:
            xdata = line.get_xdata()
            ydata = line.get_ydata()
            if not (type(xdata) == list):
                if xdata.any():
                    _top = np.max(ydata[((xdata > min_) & (xdata < max_))])
                    _top *= 1.05
                    if _top > top:
                        top = _top
        if not top:
            top = 10
        ax.set_ylim(0, top)

    def reset_figures(self):
        """Resets figure data"""

        self.ms1Ln1.set_data([], [])
        self.ms1Ln2.set_data([], [])
        self.ms1Ln3.set_data([], [])
        self.ms1Ln4.set_data([], [])
        self.ms1iso2.set_data([], [])
        self.ms1iso4.set_data([], [])
        self.ms1Lim1.set_data([], [])
        self.ms1Lim3.set_data([], [])
        self.ms1Ax1.set_xlim(0, 10)
        self.ms1Ax2.set_xlim(100, 1000)
        self.ms1Ax3.set_xlim(0, 10)
        self.ms1Ax4.set_xlim(100, 1000)
        self.ms1Ax1.info.set_text('m/z: -')
        self.ms1Ax2.info.set_text('rt: -')
        self.ms1Ax3.info.set_text('m/z: -')
        self.ms1Ax4.info.set_text('rt: -')
        for txt in reversed(self.ms1Ax2.texts):
            txt.remove()
        for txt in reversed(self.ms1Ax4.texts):
            txt.remove()
        self.ms1Canvas.draw_idle()
        self._ms1_zoom_active()

    def zoomReset(self, event):
        if (event.dblclick and event.button == MouseButton.RIGHT):
            pepind = self.PepList.currentRow()
            chgind = self.Ms1Charge.currentIndex()
            chg = int(self.Ms1Charge.currentText())
            _mz = self.peptides[pepind].mzs[chgind]
            _mod_mz = self.peptides[pepind].mod_mzs[chgind]

            try:
                if not hasattr(event.inaxes, 'name'):
                    pass
                elif event.inaxes.name == 'ms1Ax1':
                    self.zoomAx1(*self.data.time_range)
                elif event.inaxes.name == 'ms1Ax2':
                    self.zoomAx2(_mz-(4/chg), _mz+(6/chg))
                elif event.inaxes.name == 'ms1Ax3':
                    self.zoomAx3(*self.data.time_range)
                elif event.inaxes.name == 'ms1Ax4':
                    self.zoomAx4(_mod_mz-(4/chg), _mod_mz+(6/chg))
            except IndexError:
                print('No data selected.')

    def label_peaks(self, axis_, line_, min_=None, max_=None, threshold=0.05, mute=False):
        data_ = line_.get_data()
        if not min_:
            min_ = data_[0][0]
        if not max_:
            max_ = data_[0][-1]
        idx = np.where((data_[0] >= min_) & (data_[0] <= max_))[0]
        peaks_ = self.data.detect_peaks((data_[0][idx], data_[1][idx]), threshold=data_[1][idx].max()*threshold)
        if not mute:
            for txt in reversed(axis_.texts):
                txt.remove()
            for peak_ in peaks_:
                axis_.text(data_[0][idx][peak_], data_[1][idx][peak_],
                                 "%.2f" % data_[0][idx][peak_], ha='center', va='bottom', clip_on=True)
        return data_[0][idx][peaks_], data_[1][idx][peaks_]

    def zoomAx1(self, min_, max_):
        self._autoscale_y(self.ms1Ax1, min_, max_)
        self.ms1Canvas.draw_idle()

    def zoomAx2(self, min_, max_):
        peaks = self.label_peaks(self.ms1Ax2, self.ms1Ln2, min_, max_)
        self._autoscale_y(self.ms1Ax2, min_, max_)
        self.ms1Canvas.draw_idle()

    def zoomAx3(self, min_, max_):
        self._autoscale_y(self.ms1Ax3, min_, max_)
        self.ms1Canvas.draw_idle()

    def zoomAx4(self, min_, max_):
        peaks = self.label_peaks(self.ms1Ax4, self.ms1Ln4, min_, max_)
        self._autoscale_y(self.ms1Ax4, min_, max_)
        self.ms1Canvas.draw_idle()

    @wait_effect
    def combineAx1(self, min_, max_):
        self.status('Combining spectra...')
        chg = int(self.Ms1Charge.currentText())
        _mz = self.peptides[self.PepList.currentRow()].mzs[self.Ms1Charge.currentIndex()]
        spectrum = self.data.spectrum(min_, max_, mz_range=(_mz-(4/chg), _mz+(6/chg)))
        self.ms1Ln2.set_data(spectrum[0], spectrum[1])
        peaks = self.label_peaks(self.ms1Ax2, self.ms1Ln2)
        ma_peaks = self.label_peaks(self.ms1Ax2, self.ms1Ln2, _mz-0.1, _mz+0.1, 0.01, mute=True)
        if ma_peaks[0].size > 0:
            _nearest = np.abs(_mz - ma_peaks[0]).argmin() # nearest peak to mass
            if self.TolUnit.currentText() == 'mmu':
                  self.PepMzErr.setText('{:.2f} mmu'.format((_mz - ma_peaks[0][_nearest])*1e3))
            else:
                  self.PepMzErr.setText('{:.2f} ppm'.format(((_mz - ma_peaks[0][_nearest])/_mz)*1e6))
        else:
            self.PepMzErr.setText('N/A')
        self.ms1iso2.set_data([], [])
        self._autoscale_y(self.ms1Ax2, _mz-(4/chg), _mz+(6/chg))
        self.ms1Ax2.info.set_text('rt: {:.1f}-{:.1f}'.format(min_, max_))
        self.ms1iso2.set_data(self.peptides[self.PepList.currentRow()].isotopes[0]/chg + PROTON, self.peptides[self.PepList.currentRow()].isotopes[1]/1.05)
        self.ms1Lim1.set_data([min_, max_], 0)
        self.ms1Canvas.draw_idle()
        self._ms1_zoom_active()
        self.status('Spectra combined!')

    @wait_effect
    def combineAx3(self, min_, max_):
        self.status('Combining spectra...')
        chg = int(self.Ms1Charge.currentText())
        _mz = self.peptides[self.PepList.currentRow()].mod_mzs[self.Ms1Charge.currentIndex()]
        spectrum = self.data.spectrum(min_, max_, mz_range=(_mz-(4/chg), _mz+(6/chg)))
        self.ms1Ln4.set_data(spectrum[0], spectrum[1])
        peaks = self.label_peaks(self.ms1Ax4, self.ms1Ln4)
        ma_peaks = self.label_peaks(self.ms1Ax4, self.ms1Ln4, _mz-0.1, _mz+0.1, 0.01, mute=True)
        if ma_peaks[0].size > 0:
            _nearest = np.abs(_mz - ma_peaks[0]).argmin() # nearest peak to mass
            if self.TolUnit.currentText() == 'mmu':
                  self.PepModMzErr.setText('{:.2f} mmu'.format((_mz - ma_peaks[0][_nearest])*1e3))
            else:
                  self.PepModMzErr.setText('{:.2f} ppm'.format(((_mz - ma_peaks[0][_nearest])/_mz)*1e6))
        else:
            self.PepModMzErr.setText('N/A')
        self.ms1iso4.set_data([], [])
        self._autoscale_y(self.ms1Ax4, _mz-(4/chg), _mz+(6/chg))
        self.ms1Ax4.info.set_text('rt: {:.1f}-{:.1f}'.format(min_, max_))
        self.ms1iso4.set_data(self.peptides[self.PepList.currentRow()].mod_isotopes[0]/chg + PROTON, self.peptides[self.PepList.currentRow()].mod_isotopes[1]/1.05)
        self.ms1Lim3.set_data([min_, max_], 0)
        self.ms1Canvas.draw_idle()
        self._ms1_zoom_active()
        self.status('Spectra combined!')

    @wait_effect
    def extractAx2(self, min_, max_):
        self.status('Extracting chromatogram...')
        self._mz_ranges[0] = (min_, max_)
        self.ms1Ln1.set_data(self.data.chromatogram(min_, max_))
        self._autoscale_y(self.ms1Ax1, *self.data.time_range)
        self.ms1Ax1.info.set_text('m/z: {:.2f}-{:.2f}'.format(min_, max_))
        self.ms1Canvas.draw_idle()
        self._ms1_zoom_active()
        self.status('Chromatogram extracted!', 500)

    @wait_effect
    def extractAx4(self, min_, max_):
        self.status('Extracting chromatogram...')
        self._mz_ranges[1] = (min_, max_)
        self.ms1Ln3.set_data(self.data.chromatogram(min_, max_))
        self._autoscale_y(self.ms1Ax3, *self.data.time_range)
        self.ms1Ax3.info.set_text('m/z: {:.2f}-{:.2f}'.format(min_, max_))
        self.ms1Canvas.draw_idle()
        self._ms1_zoom_active()
        self.status('Chromatogram extracted!', 500)

    @wait_effect
    def areaAx1(self, min_, max_):
        pepind = self.PepList.currentRow()
        chg = int(self.Ms1Charge.currentText())
        self.status('Calculating area...')
        self.peptides[pepind].rt_ranges[0] = (min_, max_)
        self.peptides[pepind].mz_ranges[0] = self._mz_ranges[0]

        if self.peptides[pepind].chg != chg:
            self.peptides[pepind].areas[1] = 0
            self.peptides[pepind].rt_ranges[1] = []
            self.peptides[pepind].mz_ranges[1] = []
            self.peptides[pepind].detected = -1
        self.peptides[pepind].chg = chg
        self.peptides[pepind].areas[0] = self.data.get_area((min_, max_), self._mz_ranges[0])
        self.UnmodArea.setText(str(self.peptides[pepind].areas[0]))
        self.PepList.currentItem().setIcon(self.emblem)

        if self.peptides[pepind].detected == -1:
            self.peptides[pepind].detected = 0
        else:
            self.peptides[pepind].detected = 2
            self.FMod.setText('{:.2f}'.format(self.peptides[pepind].areas[1]/(self.peptides[pepind].areas[1]+self.peptides[pepind].areas[0])))
        self.status('Area calculated!', 500)
        self.coverExp.update(range(*self.peptides[pepind].id))
        self.CoverExp.setText('{:.0f} %'.format(
            len(self.coverExp)*100/len(re.sub('[^ACDEFGHIKLMNPQRSTVWY]', '', self.sequence))))
        self.ms1Lim1.set_data([min_, max_], 0)
        self.ms1Canvas.draw_idle()
        self._ms1_zoom_active()
        self.SAVED = False

    @wait_effect
    def areaAx3(self, min_, max_):
        pepind = self.PepList.currentRow()
        chg = int(self.Ms1Charge.currentText())
        self.status('Calculating area...')
        self.peptides[pepind].rt_ranges[1] = (min_, max_)
        self.peptides[pepind].mz_ranges[1] = self._mz_ranges[1]

        if self.peptides[pepind].chg != chg:
            self.peptides[pepind].areas[0] = 0
            self.peptides[pepind].rt_ranges[0] = []
            self.peptides[pepind].mz_ranges[0] = []
            self.peptides[pepind].detected = -1
        self.peptides[pepind].chg = chg
        self.peptides[pepind].areas[1] = self.data.get_area((min_, max_), self._mz_ranges[1])
        self.ModArea.setText(str(self.peptides[pepind].areas[1]))
        self.PepList.currentItem().setIcon(self.emblem)

        if self.peptides[pepind].detected == -1:
            self.peptides[pepind].detected = 1
        else:
            self.peptides[pepind].detected = 2
            self.FMod.setText('{:.2f}'.format(self.peptides[pepind].areas[1]/(self.peptides[pepind].areas[1]+self.peptides[pepind].areas[0])))
        self.status('Area calculated!', 500)
        self.coverExp.update(range(*self.peptides[pepind].id))
        self.CoverExp.setText('{:.0f} %'.format(
            len(self.coverExp)*100/len(re.sub('[^ACDEFGHIKLMNPQRSTVWY]', '', self.sequence))))
        self.ms1Lim3.set_data([min_, max_], 0)
        self.ms1Canvas.draw_idle()
        self._ms1_zoom_active()
        self.SAVED = False

    def remove_assignment(self, assignment):
        """Remove assigned from PepList and peptides"""

        pepind = self.PepList.currentRow()
        self.PepList.currentItem().setIcon(Qtg.QIcon())
        self.peptides[pepind].chg = None
        self.peptides[pepind].detected = -1
        self.peptides[pepind].mz_ranges = [[], []]
        self.peptides[pepind].rt_ranges = [[], []]
        self.peptides[pepind].areas = [0, 0]
        self.RemoveAssignment.setEnabled(False)
        self.SAVED = False

    def update_analysis(self):
        """Fill bar graph with fractional modification data."""

        treatment = []
        if self.HoloList.count():
            treatment.extend([[], []])
            file_list_basename = [os.path.basename(x) for x in self.file_list]
            for i in range(self.ApoList.count()):
                item = self.ApoList.item(i).text()
                if item in file_list_basename:
                    treatment[0].append(file_list_basename.index(item))

            for i in range(self.HoloList.count()):
                item = self.HoloList.item(i).text()
                if item in file_list_basename:
                    treatment[1].append(file_list_basename.index(item))

            treatment[0].sort()
            treatment[1].sort()

        self.project['treatment'] = treatment

        self.barAx1.clear()
        self.barAx2.clear()
        threshold = self.ThreshVal.value()
        significance = self.SigVal.value()

        try:
            fmod = np.array(self.project['fractional mod'])
            pep_ids = np.array(self.project['peptides'])
            num_peptides = len(fmod[0])
            idx = np.arange(num_peptides)

            if treatment:
                treatment = self.project['treatment']
                mean1 = fmod[treatment[0]].mean(axis=0)
                mean2 = fmod[treatment[1]].mean(axis=0)
                mask1 = mean1 >= threshold
                mask2 = mean2 >= threshold
                stdev1 = fmod[treatment[0]].std(axis=0)
                stdev2 = fmod[treatment[1]].std(axis=0)
                tvalue, pvalue = stats.ttest_ind(fmod[treatment[0]], fmod[treatment[1]], axis=0)
                sig_mask = (pvalue <= significance) & (mask1 | mask2)
                sign_mask = np.sign(mean1[sig_mask]-mean2[sig_mask])
                
                extent_mean, extent_stdev = calc_extent_change(mean1, stdev1, mean2, stdev2)
                
                mask3 = extent_mean > 0
                mean = -extent_mean

                width = 0.3
                offset = width/2

                if not self.DiffPlotChk.isChecked():

                    self.barAx1.bar(idx[mask1]-offset, mean1[mask1], width, yerr=stdev1[mask1], color=ORANGE, error_kw={'elinewidth':1}, capsize=2)
                    self.barAx1.bar(idx[~mask1]-offset, mean1[~mask1], width, yerr=stdev1[~mask1], color='lightgray', error_kw={'elinewidth':1}, capsize=2)
                    self.barAx1.bar(idx[mask2]+offset, mean2[mask2], width, yerr=stdev2[mask2], color=BLUE, error_kw={'elinewidth':1}, capsize=2)
                    self.barAx1.bar(idx[~mask2]+offset, mean2[~mask2], width, yerr=stdev2[~mask2], color='lightgray', error_kw={'elinewidth':1}, capsize=2)
                    self.barAx1.plot(idx[sig_mask], np.maximum(mean1, mean2)[sig_mask]+0.1, 'ko', markersize=5)
                    self.barAx1.axhline(threshold, ls=':', lw=1, c='k')

                    widths = pep_ids[:,1]-pep_ids[:,0]+1

                    if not self.UnityChk.isChecked():
                        heights1 = np.sqrt(mean1/widths)
                        bottoms1 = np.zeros_like(mean1)
                        heights2 = -np.sqrt(mean2/widths)
                        bottoms2 = bottoms1
                    else:
                        heights1 = np.ones_like(mean1)*0.9
                        bottoms1 = get_figure_rows(pep_ids)
                        heights2 = -heights1
                        bottoms2 = -bottoms1

                    self.barAx2.bar(pep_ids[:,0], heights1, widths, bottoms1, align='edge', color=get_colours(BLUE, mean1, threshold))
                    self.barAx2.bar(pep_ids[:,0], heights2, widths, bottoms2, align='edge', color=get_colours(ORANGE, mean2, threshold))
                    self.barAx2.axhline(0, lw=1, c='0.5')
                    self.barAx1.set_ylabel('Fractional Modification')

                else:

                    self.barAx1.bar(idx[mask3], extent_mean[mask3], 0.4, yerr=extent_stdev[mask3], color=GREEN, error_kw={'elinewidth':1}, capsize=2)
                    self.barAx1.bar(idx[~mask3], extent_mean[~mask3], 0.4, yerr=extent_stdev[~mask3], color=REDPURPLE, error_kw={'elinewidth':1}, capsize=2)
                    self.barAx1.axhline(0, lw=1, c='0.5')
                    self.barAx1.yaxis.set_major_formatter(mtick.PercentFormatter(xmax=1))
                    self.barAx1.set_ylabel('Modification Difference %')

                    widths = pep_ids[:,1]-pep_ids[:,0]+1
                    extent_max = np.abs(extent_mean).max()

                    if not self.UnityChk.isChecked():
                        heights = np.sign(extent_mean)*np.sqrt(np.abs(extent_mean)/widths)
                        bottoms = np.zeros_like(extent_mean)
                    else:
                        heights = np.ones_like(extent_mean)*0.9
                        bottoms = get_figure_rows(pep_ids)
                        self.barAx2.set_ylim(-0.05, bottoms.max()+1.05)

                    self.barAx2.bar(pep_ids[:,0][mask3], heights[mask3], widths[mask3], bottoms[mask3], align='edge', color=get_colours(GREEN, extent_mean[mask3]/extent_max, threshold))
                    self.barAx2.bar(pep_ids[:,0][~mask3], heights[~mask3], widths[~mask3], bottoms[~mask3], align='edge', color=get_colours(REDPURPLE, -extent_mean[~mask3]/extent_max, threshold))
                    self.barAx2.axhline(0, lw=1, c='0.5')

            else:
                mean = fmod.mean(axis=0)
                sig_mask = mean >= threshold
                sign_mask = np.ones_like(mean)[sig_mask]
                stdev = fmod.std(axis=0)
                width = 0.4

                self.barAx1.bar(idx[sig_mask], mean[sig_mask], width, yerr=stdev[sig_mask], color=ORANGE, error_kw={'elinewidth':1}, capsize=2)
                self.barAx1.bar(idx[~sig_mask], mean[~sig_mask], width, yerr=stdev[~sig_mask], color=BLUE, error_kw={'elinewidth':1}, capsize=2)
                self.barAx1.axhline(threshold, ls=':', lw=1, c='k')
                self.barAx1.set_ylabel('Fractional Modification')

                widths = pep_ids[:,1]-pep_ids[:,0]+1

                if not self.UnityChk.isChecked():
                    heights = np.sqrt(mean/widths)
                    bottoms = -heights/2
                else:
                    heights = np.ones_like(mean)*0.9
                    bottoms = get_figure_rows(pep_ids)
                    self.barAx2.set_ylim(-0.05, bottoms.max()+1.05)

                self.barAx2.bar(pep_ids[:,0][sig_mask], heights[sig_mask], widths[sig_mask], bottoms[sig_mask], align='edge', color=get_colours(ORANGE, mean[sig_mask], 0.01))
                self.barAx2.bar(pep_ids[:,0][~sig_mask], heights[~sig_mask], widths[~sig_mask], bottoms[~sig_mask], align='edge', color=get_colours(BLUE, mean[~sig_mask], 0.01))
                self.barAx2.axhline(0, lw=1, c='0.5', zorder=0)

            labels = ['{}-{}'.format(*pep) for pep in self.project['peptides']]
            self.barAx1.set_xlim(-0.5, idx[-1]+0.5)
            self.barAx1.set_xticks(idx-0.5)
            self.barAx1.set_xticks(idx, minor=True)
            self.barAx1.set_xticklabels([])
            self.barAx1.set_xticklabels(labels, rotation='vertical', minor=True)
            self.barAx2.set_xlim(1, len(re.sub('[^ACDEFGHIKLMNPQRSTVWY]', '', self.sequence))+1.05)
            self.barAx2.axis('off')
            self.barCanvas.draw_idle()
            #self.Tabs.setCurrentIndex(1)

            self.actionExportCSV.setEnabled(True)

            if self.PDB_LOADED:
                self.update_bfactor(np.array(self.project['peptides']), sig_mask, sign_mask, mean, continuous=self.NGL_continuous.isChecked())

        except KeyError:
            self.status('No fractional modification data detected.', 2000)

    def on_bar_pick(self, event):
        """Activate show_bar_info on pick event for bar in bar chart"""

        if event.inaxes == self.barAx1:
            idx = int(np.floor(event.xdata+0.5))
            if self.barAx1.selected:
                self.barAx1.selected.remove()
                self.barAx2.selected.remove()
            self.barAx1.selected = self.barAx1.axvspan(idx-0.5, idx+0.5, -0.1, 1.1, color='0.95', ec='k', zorder=0)
            self.barAx2.selected = self.barAx2.axvspan(self.project['peptides'][idx][0], self.project['peptides'][idx][1]+1, -0.1, 1.1, color=(0.8, 0.8, 0.8, 0.2), ec='k')
            self.barCanvas.draw_idle()
            self.show_bar_info(idx)
        elif event.inaxes == self.barAx2:
            idx = int(np.floor(event.xdata))
            pep_ids = np.array(self.project['peptides'])
            bar = list(np.where((pep_ids[:,0] <= idx) & (pep_ids[:,1] > idx))[0])
            if bar:
                if self.barAx1.selected:
                    self.barAx1.selected.remove()
                    self.barAx2.selected.remove()
                self.barAx1.selected = self.barAx1.axvspan(bar[0]-0.5, bar[0]+0.5, -0.1, 1.1, color='0.95', ec='k', zorder=0)
                self.barAx2.selected = self.barAx2.axvspan(self.project['peptides'][bar[0]][0], self.project['peptides'][bar[0]][1]+1, -0.1, 1.1, color=(0.8, 0.8, 0.8, 0.2), ec='k')
                self.barCanvas.draw_idle()
                self.show_bar_info(bar[0])
        else:
            if self.barAx1.selected:
                self.barAx1.selected.remove()
                self.barAx2.selected.remove()
            self.barAx1.selected = []
            self.barAx2.selected = []
            self.barCanvas.draw_idle()

    def show_bar_info(self, idx):
        """On hover/pick show fractional modification data for given peptide

        Peptide sequence
        (Treatment #)
        fmod 1, fmod 2, .. fmod n, mean, std
        """

        fmod = np.array(self.project['fractional mod'])

        pep_id = self.project['peptides'][idx]

        self.PeptideLbl.setText('{}-{}-{}'.format(pep_id[0], ''.join(parser.parse(self.project['sequence'])[pep_id[0]-1:pep_id[1]]), pep_id[1]))

        _text = ''
        if self.project['treatment']:
            treatments = self.project['treatment']
            _text += 'Apo\n'
            for j in treatments[0]:
                _text += '{:.2f}    '.format(fmod[j][idx])
            _mean1 = fmod[treatments[0]][:,idx].mean()
            _stdev1 = fmod[treatments[0]][:,idx].std()
            _text += '    Average: {:.2f} ± {:.2f}\n\n'.format(_mean1, _stdev1)
            _text += 'Holo\n'
            for j in treatments[1]:
                _text += '{:.2f}    '.format(fmod[j][idx])
            _mean2 = fmod[treatments[1]][:,idx].mean()
            _stdev2 = fmod[treatments[1]][:,idx].std()
            _text += '    Average: {:.2f} ± {:.2f}\n\n'.format(_mean2, _stdev2)
            _text += 'Extent of change:'
            _text += ' {:.1%} ± {:.1%}'.format(*calc_extent_change(_mean1, _stdev1, _mean2, _stdev2))
            _text += '    p = {:.2g}\n'.format(stats.ttest_ind(fmod[treatments[0]][:,idx], fmod[treatments[1]][:,idx])[1])

        else:
            for i in range(len(fmod)):
                _text += '{:.2f}    '.format(fmod[i][idx])
            _text += '    Average: {:.2f} ± {:.2f}\n'.format(fmod[:,idx].mean(), fmod[:,idx].std())

        self.BarInfo.setText(_text)

    def load_pdb(self, pdb_file=None):
        """Open .pdb and process with PDB object"""
        if not pdb_file:
            file_name = Qtw.QFileDialog.getOpenFileName(
                self, 'Select PDB file', '', 'Protein Data Bank Files (*.pdb)')[0]
            if file_name:
                pass
        else:
            file_name = pdb_file

        try:
            if not os.path.exists(file_name):
                raise AssertionError('PDB file {} not found!'.format(file_name))
            self.pdb = PDB(file_name)
            self.PdbName.setText(os.path.split(file_name)[1])
            self.PDB_LOADED = True
            self.Tabs.setCurrentIndex(2)
            self.update_ngl(file_name)
        except Exception as e:
            self.status(str(e), 2000)
            print(e)

    def update_bfactor(self, ids, sig_mask, sign_mask, mean, continuous=False):
        """Calculate B-factors from fmod data and apply to self.pdb

        Not detected    -2  Wheat
        Insignificant    0  Grey
        Significant      1  Red
        Significant -   -1  Blue

        """

        # Only handle normal amino acid code, removes lowercase modX
        strip_sequence = re.sub('[^ACDEFGHIKLMNPQRSTVWY]', '', self.sequence)

        bfactors = np.zeros(len(strip_sequence))

        bfactors[:] = -2

        if not continuous:

            for peptide in ids[~sig_mask]:
                bfactors[peptide[0]:peptide[1]] = 0

            for i, peptide in enumerate(ids[sig_mask]):
                if sign_mask[i] == 1:
                    bfactors[peptide[0]-1:peptide[1]] = 1
                else:
                    bfactors[peptide[0]-1:peptide[1]] = -1

        else:

            for i, peptide in enumerate(ids):
                bfactors[peptide[0]-1:peptide[1]] = mean[i]
                #bfactors[peptide[0]-1:peptide[1]] = -mean[i]

        try:
            self.pdb.bfactor_by_residue(strip_sequence, bfactors)
            self.save_pdb(self.pdb.filename)
            self.update_ngl(self.pdb.filename)
        except ValueError as error:
            self.status(error.args[0], 2000)


    def save_pdb(self, filename=None):
        if filename == None:
            self.pdb_file = Qtw.QFileDialog.getSaveFileName(
                self, 'Save PDB file', '', 'Protein Data Bank Files (*.pdb)')[0]
        else:
            self.pdb_file = filename

        if not os.path.splitext(self.pdb_file)[1]:
            self.pdb_file += '.pdb'
        self.pdb.write(self.pdb_file)

    def update_ngl(self, filename):
        """Update NGLViewer with modified PDB file."""

        ngl_path = Qtc.QUrl.fromLocalFile(os.path.join(BUNDLE_DIR, 'ngl.js')).toString()
        pdb_path = Qtc.QUrl.fromLocalFile(filename).toString()

        HTML = '<html><body>'
        HTML += '<script src="{}"></script>'.format(ngl_path)
        HTML += '''
  <script>
    document.addEventListener("DOMContentLoaded", function () {
      var stage = new NGL.Stage("viewport", {tooltip: false, cameraType:'orthographic'});
      var myBfactor = NGL.ColormakerRegistry.addScheme(function (params) {
          this.parameters.scale = "RdYlBu"
          this.parameters.reverse = true
          this.parameters.domain = [-1.1, 1.1]
          this.Bscale = this.getScale()
          this.atomColor = function (atom) {
            if (atom.bfactor == -2) {
              return 0xFFFFFF // Grey
            } else {
              return this.Bscale(atom.bfactor);
            }
          };
        });'''
        HTML += 'stage.loadFile("{}").then('.format(pdb_path)
        HTML += '''function (protein) {
        var cartoon = protein.addRepresentation("cartoon", {colorScheme: myBfactor, aspectRatio: 7, smoothedSheet:true});
        var hetatom = protein.addRepresentation("licorice", {sele: "hetero", multipleBond: true, radius: 0.3});
        var surface = protein.addRepresentation("surface", {colorScheme: myBfactor, surfaceType: "av", visible:false});
        protein.autoView();
        function toggleCartoon () {
            surface.setVisibility(!surface.visible)
        };
        function toggleHetatom () {
            hetatom.setVisibility(!hetatom.visible)
        };
        stage.keyControls.add("c", toggleCartoon);
        stage.keyControls.add("h", toggleHetatom);

        });

      var tooltip = document.createElement("div");
      Object.assign(tooltip.style, {
        display: "none",
        position: "absolute",
        zIndex: 10,
        pointerEvents: "none",
        backgroundColor: "rgba(0, 0, 0, 0.6)",
        color: "lightgrey",
        padding: "0.5em",
        fontFamily: "sans-serif"
      });
      stage.viewer.container.appendChild(tooltip);

      stage.signals.hovered.add(function (pickingProxy) {
       if (pickingProxy && (pickingProxy.atom || pickingProxy.bond)){
         var atom = pickingProxy.atom || pickingProxy.closestBondAtom;
         var cp = pickingProxy.canvasPosition;
         tooltip.innerText = atom.resname + " " + atom.resno;
         tooltip.style.bottom = cp.y + 0 + "px";
         tooltip.style.left = cp.x + 0 + "px";
         tooltip.style.display = "block";
       }else{
         tooltip.style.display = "none";
       }
      });

    });
  </script>
  '''
        HTML += '<div id="viewport" style="width:100%; height:100%;"></div></body></html>'

        self.NGLViewer.setHtml(HTML)


    def _ms1_zoom_active(self):
        self.ExtractBtn.setChecked(False)
        self.ZoomBtn.setChecked(True)
        self.IntegrateBtn.setChecked(False)
        self.zoom1.active = True
        self.zoom2.active = True
        self.zoom3.active = True
        self.zoom4.active = True
        self.area1.active = False
        self.area3.active = False
        self.combine1.active = False
        self.combine3.active = False
        self.extract2.active = False
        self.extract4.active = False

    def _ms1_extract_active(self):
        self.ExtractBtn.setChecked(True)
        self.ZoomBtn.setChecked(False)
        self.IntegrateBtn.setChecked(False)
        self.zoom1.active = False
        self.zoom2.active = False
        self.zoom3.active = False
        self.zoom4.active = False
        self.area1.active = False
        self.area3.active = False
        self.combine1.active = True
        self.combine3.active = True
        self.extract2.active = True
        self.extract4.active = True

    def _ms1_area_active(self):
        self.ExtractBtn.setChecked(False)
        self.ZoomBtn.setChecked(False)
        self.IntegrateBtn.setChecked(True)
        self.zoom1.active = False
        self.zoom2.active = False
        self.zoom3.active = False
        self.zoom4.active = False
        self.area1.active = True
        self.area3.active = True
        self.combine1.active = False
        self.combine3.active = False
        self.extract2.active = False
        self.extract4.active = False

    def keyPressEvent(self, event):
        """Reimplement the keyPressEvent() event handler to include modifying
        SpanSelector activity"""

        if self.Tabs.currentIndex() == 0:
            if event.key() == Qt.Key_E:
                self._ms1_extract_active()

            elif event.key() == Qt.Key_I:
                self._ms1_area_active()

            elif event.key() == Qt.Key_Z:
                self._ms1_zoom_active()

    def load_settings(self):
        """Parse settings stored in self.settings to relevant objects/widgets"""

        if 'pwiz' in self.settings.childKeys():
            self.pwiz = self.settings.value('pwiz')
        else:
            self.pwiz = None

        # Get enzymes
        enzymes.clear()
        if 'Enzymes' in self.settings.childGroups():
            self.settings.beginGroup('Enzymes')
            enz_keys = self.settings.childKeys()
            for key in enz_keys:
                enzymes[key] = self.settings.value(key)
            self.settings.endGroup()
        else:
            for enzyme in default_enzymes:
                enzymes[enzyme] = default_enzymes[enzyme]

        # Get Mods
        mods.clear()
        if 'Modifications' in self.settings.childGroups():
            self.settings.beginGroup('Modifications')
            mod_count = self.settings.beginReadArray('mod')
            for mod in range(mod_count):
                self.settings.setArrayIndex(mod)
                name = self.settings.value('name')
                mods[name] = {}
                mods[name]['ID'] = self.settings.value('id')
                mods[name]['Gain'] = self.settings.value('gain')
                mods[name]['Loss'] = self.settings.value('loss')
                mods[name]['Mass'] = float(self.settings.value('mass'))
                mods[name]['Residues'] = self.settings.value('residue')
            self.settings.endArray()
            self.settings.endGroup()
        else:
            for mod in default_mods:
                mods[mod] = default_mods[mod]

    def save_settings(self):
        """Store current settings to pepFoot.ini"""

        if self.pwiz:
            self.settings.setValue('pwiz', self.pwiz)

        # Save Enzymes

        self.settings.beginGroup('Enzymes')
        self.settings.remove("")
        for enzyme in enzymes:
            self.settings.setValue(enzyme, enzymes[enzyme])
        self.settings.endGroup()

        # Save Mods

        self.settings.beginGroup('Modifications')
        self.settings.remove("")
        mod_keys = [*mods]
        mod_count = len(mod_keys)
        self.settings.beginWriteArray('mod')
        for mod in range(mod_count):
            self.settings.setArrayIndex(mod)
            name = mod_keys[mod]
            self.settings.setValue('name', name)
            self.settings.setValue('id', mods[name]['ID'])
            self.settings.setValue('gain', mods[name]['Gain'])
            self.settings.setValue('loss', mods[name]['Loss'])
            self.settings.setValue('mass', mods[name]['Mass'])
            self.settings.setValue('residue', mods[name]['Residues'])
        self.settings.endArray()
        self.settings.endGroup()

    def about_dialog(self, tab):
        about = About(tab)
        about.exec_()
        about.show()

    def preferences_dialog(self, tab):
        prefs = Preferences(tab, self)
        prefs.exec_()
        prefs.show()
        
    def open_documentation(self):
        url = Qtc.QUrl("https://github.com/jbellamycarter/pepfoot/wiki")
        Qtg.QDesktopServices.openUrl(url)

    def closeEvent(self, event):
        """Reimplement the closeEvent() event handler to include a 'Question'
        dialog with options on how to proceed - Save, Close, Cancel buttons
        """
        self.save_settings()
        if not self.SAVED:
            reply = Qtw.QMessageBox.question(
                self, "Quit PepFoot",
                "Are you sure you want to quit? Any unsaved work will be lost.",
                Qtw.QMessageBox.Close | Qtw.QMessageBox.Cancel, Qtw.QMessageBox.Cancel)

            if reply == Qtw.QMessageBox.Close:
                event.accept()
            else:
                event.ignore()
        else:
            event.accept()


class About(Qtw.QDialog, Ui_AboutDialog):
    def __init__(self, tab):
        super(About, self).__init__()
        self.setupUi(self)
        self.Tabs.setCurrentIndex(tab)


class Preferences(Qtw.QDialog, Ui_Dialog):
    def __init__(self, tab, main_):
        super(Preferences, self).__init__()
        self.main_ = main_
        self.setupUi(self)
        self.Tabs.setCurrentIndex(tab)
        self.fill_mods()
        self.fill_enzymes()
        self.ModTable.itemSelectionChanged.connect(self.update_modbox)
        self.ModAdd.clicked.connect(self.update_ModTable)
        self.ModRemove.clicked.connect(self.remove_mod)
        self.ModID.textChanged.connect(self.check_id)
        self.EnzymeTable.itemSelectionChanged.connect(self.update_enzymebox)
        self.EnzymeAdd.clicked.connect(self.update_EnzymeTable)
        self.EnzymeRemove.clicked.connect(self.remove_enzyme)
        self.EnzymeName.textChanged.connect(self.check_name)
        self.accepted.connect(self.handleOK)

    def fill_enzymes(self):
        self.EnzymeTable.setSortingEnabled(False)
        for i, e in enumerate(enzymes):
            clv = enzymes[e].find('(?=')
            if clv == -1:
                cterm = ''
                nterm = enzymes[e][1:-1]
            elif clv == 0:
                cterm = enzymes[e][4:-2]
                nterm = ''
            else:
                nterm = enzymes[e][1:clv-1]
                cterm = enzymes[e][clv+4:-2]
            self.EnzymeTable.insertRow(i)
            self.EnzymeTable.setItem(i, 0, Qtw.QTableWidgetItem(e))
            self.EnzymeTable.setItem(i, 1, Qtw.QTableWidgetItem(nterm))
            self.EnzymeTable.setItem(i, 2, Qtw.QTableWidgetItem(cterm))
        self.EnzymeTable.setSortingEnabled(True)
        self.EnzymeTable.sortItems(0)
        self.EnzymeTable.horizontalHeader().setSectionResizeMode(Qtw.QHeaderView.ResizeToContents)

    def fill_mods(self):
        self.ModTable.setSortingEnabled(False)
        for i, e in enumerate(mods):
            self.ModTable.insertRow(i)
            self.ModTable.setItem(i, 0, Qtw.QTableWidgetItem(mods[e]['ID']))
            self.ModTable.setItem(i, 1, Qtw.QTableWidgetItem(mods[e]['Loss']))
            self.ModTable.setItem(i, 2, Qtw.QTableWidgetItem(mods[e]['Gain']))
            self.ModTable.setItem(i, 3, Qtw.QTableWidgetItem(
                "{:.4f}".format(mods[e]['Mass'])))
            self.ModTable.setItem(
                i, 4, Qtw.QTableWidgetItem(mods[e]['Residues']))
            self.ModTable.setItem(i, 5, Qtw.QTableWidgetItem(e))
        self.ModTable.setSortingEnabled(True)
        self.ModTable.sortItems(0)
        self.ModTable.horizontalHeader().setSectionResizeMode(Qtw.QHeaderView.ResizeToContents)

    def check_id(self):
        if self.ModID.text():
            self.ModAdd.setEnabled(True)
        else:
            self.ModAdd.setEnabled(False)

    def check_name(self):
        if self.EnzymeName.text():
            self.EnzymeAdd.setEnabled(True)
        else:
            self.EnzymeAdd.setEnabled(False)

    def remove_mod(self):
        self.ModTable.removeRow(self.ModTable.currentRow())
        self.ModTable.clearSelection()
        self.ModRemove.setEnabled(False)

    def update_modbox(self):
        row = self.ModTable.currentRow()
        if not row == -1:
            self.ModID.setText(self.ModTable.item(row, 0).text())
            self.ModLoss.setText(self.ModTable.item(row, 1).text())
            self.ModGain.setText(self.ModTable.item(row, 2).text())
            self.ModResi.setText(self.ModTable.item(row, 4).text())
            self.ModDesc.setText(self.ModTable.item(row, 5).text())
            self.ModAdd.setEnabled(True)
            self.ModRemove.setEnabled(True)

    def update_ModTable(self):
        existing = self.ModTable.findItems(self.ModID.text(), Qt.MatchFlags(8))
        shift = mass.calculate_mass(formula=self.ModGain.text(
        )) - mass.calculate_mass(formula=self.ModLoss.text())
        self.ModTable.setSortingEnabled(False)
        if existing:
            row = existing[0].row()
            self.ModTable.item(row, 1).setText(self.ModLoss.text())
            self.ModTable.item(row, 2).setText(self.ModGain.text())
            self.ModTable.item(row, 3).setText('{:.4f}'.format(shift))
            self.ModTable.item(row, 4).setText(self.ModResi.text())
            self.ModTable.item(row, 5).setText(self.ModDesc.text())
        else:
            row = self.ModTable.rowCount()
            self.ModTable.insertRow(row)
            self.ModTable.setItem(
                row, 0, Qtw.QTableWidgetItem(self.ModID.text()))
            self.ModTable.setItem(
                row, 1, Qtw.QTableWidgetItem(self.ModLoss.text()))
            self.ModTable.setItem(
                row, 2, Qtw.QTableWidgetItem(self.ModGain.text()))
            self.ModTable.setItem(
                row, 3, Qtw.QTableWidgetItem('{:.4f}'.format(shift)))
            self.ModTable.setItem(
                row, 4, Qtw.QTableWidgetItem(self.ModResi.text()))
            self.ModTable.setItem(
                row, 5, Qtw.QTableWidgetItem(self.ModDesc.text()))
        self.ModTable.setSortingEnabled(True)

    def remove_enzyme(self):
        self.EnzymeTable.removeRow(self.EnzymeTable.currentRow())
        self.EnzymeTable.clearSelection()
        self.EnzymeRemove.setEnabled(False)

    def update_enzymebox(self):
        row = self.EnzymeTable.currentRow()
        if not row == -1:
            self.EnzymeName.setText(self.EnzymeTable.item(row, 0).text())
            self.EnzymeNterm.setText(self.EnzymeTable.item(row, 1).text())
            self.EnzymeCterm.setText(self.EnzymeTable.item(row, 2).text())
            self.EnzymeAdd.setEnabled(True)
            self.EnzymeRemove.setEnabled(True)

    def update_EnzymeTable(self):
        existing = self.EnzymeTable.findItems(self.EnzymeName.text(), Qt.MatchFlags(8))
        self.EnzymeTable.setSortingEnabled(False)
        if existing:
            row = existing[0].row()
            self.EnzymeTable.item(row, 1).setText(self.EnzymeNterm.text())
            self.EnzymeTable.item(row, 2).setText(self.EnzymeCterm.text())

        else:
            row = self.EnzymeTable.rowCount()
            self.EnzymeTable.insertRow(row)
            self.EnzymeTable.setItem(row, 0, Qtw.QTableWidgetItem(self.EnzymeName.text()))
            self.EnzymeTable.setItem(row, 1, Qtw.QTableWidgetItem(self.EnzymeNterm.text()))
            self.EnzymeTable.setItem(row, 2, Qtw.QTableWidgetItem(self.EnzymeCterm.text()))
        self.EnzymeTable.setSortingEnabled(True)

    def handleOK(self):
        mods.clear()
        for row in range(self.ModTable.rowCount()):
            mods[self.ModTable.item(row, 5).text()] = {
                'Loss': self.ModTable.item(row, 1).text(),
                'Gain': self.ModTable.item(row, 2).text(),
                'Mass': float(self.ModTable.item(row, 3).text()),
                'Residues': self.ModTable.item(row, 4).text(),
                'ID': self.ModTable.item(row, 0).text()}
        self.main_.DiffMod.clear()
        self.main_.DiffMod.addItem('None')
        self.main_.DiffMod.insertSeparator(1)
        self.main_.fill_mods()

        enzymes.clear()
        for row in range(self.EnzymeTable.rowCount()):
            _name = self.EnzymeTable.item(row, 0).text()
            _nterm = self.EnzymeTable.item(row, 1).text()
            _cterm = self.EnzymeTable.item(row, 2).text()
            if _nterm and _cterm:
                enzymes[_name] = '[{}](?=[{}])'.format(_nterm, _cterm)
            elif _nterm:
                enzymes[_name] = '[{}]'.format(_nterm)
            else:
                enzymes[_name] = '(?=[{}])'.format(_cterm)
        self.main_.Enzyme.clear()
        self.main_.Enzyme.addItem('None')
        self.main_.Enzyme.insertSeparator(1)
        self.main_.Enzyme.addItems(sorted(list(enzymes)))

def run_app():
    start = time.time()
    splash_pix = Qtg.QPixmap(os.path.join(BUNDLE_DIR, 'gui', 'splash.png'))
    SPLASH = Qtw.QSplashScreen(splash_pix, Qt.WindowStaysOnTopHint)

    SPLASH.show()
    APP.setWindowIcon(APP_ICON)

    while time.time() - start < 0.5:
        time.sleep(0.001)
        APP.processEvents()

    main = Main()
    main.show()

    if len(sys.argv) > 2:
        main.open_project(sys.argv[1])

    SPLASH.finish(main)

    sys.exit(APP.exec_())

if __name__ == '__main__':
    run_app()
