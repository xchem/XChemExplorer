import getpass
import glob
import math
import os
import pickle
import sys
from datetime import datetime

from PyQt4 import QtCore, QtGui

from xce.gui_scripts import layout, layout_functions, stylesheet
from xce.lib import (
    XChemDB,
    XChemDeposit,
    XChemLog,
    XChemMain,
    XChemPANDDA,
    XChemPlots,
    XChemThread,
    XChemToolTips,
)
from xce.lib.XChemUtils import parse
from xce.lib.cluster.slurm import get_token
from xce.web import XChemWeb


class XChemExplorer(QtGui.QApplication):
    def __init__(self, args):
        # init a QApplication object to hold XCE
        QtGui.QApplication.__init__(self, args)

        # start GUI
        self.start_GUI()

        # set stylesheet - how the gui looks
        stylesheet.set_stylesheet(self)

        self.exec_()

    def start_GUI(self):
        # check http://doc.qt.io/qt-4.8/stylesheet-customizing.html#the-box-model
        # This needs moving somewhere more appropriate...
        self.headlineLabelfont = QtGui.QFont("Arial", 20, QtGui.QFont.Bold)

        layout.setup().settings(self)
        layout.setup().preferences(self)
        layout.setup().tables(self)

        # GUI setup
        self.window = QtGui.QWidget()
        self.window.setWindowTitle("XChemExplorer")
        self.screen = QtGui.QDesktopWidget().screenGeometry()

        layout.LayoutObjects().workflow(self)
        layout.LayoutObjects().main_layout(self)
        layout_functions.add_widgets_layouts(self)

        self.checkLabXChemDir()

        if os.path.isfile(os.path.join(self.database_directory, self.data_source_file)):
            self.backup_soakDB()

    def backup_soakDB(self):
        XChemMain.backup_soakDB(
            os.path.join(self.database_directory, self.data_source_file),
            self.xce_logfile,
        )

    def checkLabXChemDir(self):
        dirCheck = QtGui.QMessageBox()
        dirCheckLayout = dirCheck.layout()
        vbox = QtGui.QVBoxLayout()
        try:
            warning = (
                "Are you sure you want to launch XCE here:\n\n"
                + self.labxchem_directory_current
                + "\n\n"
                "If this is not where you should be running XCE, please close!\n"
            )
        except AttributeError:
            return
        vbox.addWidget(QtGui.QLabel(warning))
        dirCheckLayout.addLayout(vbox, 0, 0)
        dirCheck.exec_()

    # function to update datasource

    def datasource_menu_reload_samples(self):
        self.update_log.insert(
            "reading samples from data source: "
            + os.path.join(self.database_directory, self.data_source_file)
        )
        self.update_status_bar(
            "reading samples from data source: "
            + os.path.join(self.database_directory, self.data_source_file)
        )
        self.update_header_and_data_from_datasource()
        self.update_all_tables()
        self.overview_datasource_table.resizeColumnsToContents()

    # function to create new datasource
    def create_new_data_source(self):
        file_name = str(
            QtGui.QFileDialog.getSaveFileName(
                self.window, "Save file", self.database_directory
            )
        )
        # make sure that the file always has .sqlite extension
        if file_name.rfind(".") != -1:
            file_name = file_name[: file_name.rfind(".")] + ".sqlite"
        else:
            file_name = file_name + ".sqlite"
        self.db = XChemDB.data_source(file_name)
        print("==> XCE: creating new data source")
        self.db.create_empty_data_source_file()
        self.db.create_missing_columns()
        self.database_directory = file_name[: file_name.rfind("/")]
        self.data_source_file = file_name[file_name.rfind("/") + 1 :]
        self.data_source_file_label.setText(
            os.path.join(self.database_directory, self.data_source_file)
        )
        self.settings["database_directory"] = self.database_directory
        self.settings["data_source"] = self.data_source_file
        self.data_source_set = True
        self.datasource_menu_reload_samples()

    ####################################################################################
    #                                                                                  #
    #                                   DATASETS TAB                                   #
    #                                                                                  #
    ####################################################################################

    def continously_check_for_new_data_collection(self, state):
        self.timer_to_check_for_new_data_collection.timeout.connect(
            lambda: self.check_for_new_autoprocessing_or_rescore(False)
        )
        if state == QtCore.Qt.Checked:
            print("==> XCE: checking automatically every 120s for new data collection")
            self.timer_to_check_for_new_data_collection.start(120000)
        else:
            print("==> XCE: stopped checking for new data collections")
            self.timer_to_check_for_new_data_collection.stop()

    def target_selection_combobox_activated(self, text):
        self.target = str(text)

    def check_for_new_autoprocessing_or_rescore(self, rescore_only):
        self.update_log.insert("checking for new data collection")
        start_thread = False
        if rescore_only:
            # first pop up a warning message as this will overwrite all user selections
            msgBox = QtGui.QMessageBox()
            msgBox.setText(
                "*** WARNING ***\n"
                "This will overwrite all your manual selections!\n"
                "Do you want to continue?"
            )
            msgBox.addButton(QtGui.QPushButton("Yes"), QtGui.QMessageBox.YesRole)
            msgBox.addButton(QtGui.QPushButton("No"), QtGui.QMessageBox.RejectRole)
            reply = msgBox.exec_()
            if reply == 0:
                start_thread = True
            else:
                start_thread = False
        else:
            start_thread = True

        if start_thread:
            if self.target == "=== SELECT TARGET ===":
                msgBox = QtGui.QMessageBox()
                warning = (
                    "*** WARNING ***\n"
                    "Please select a target or\n"
                    'select "=== project directory ==="'
                    " if you want to read reprocessed results\n"
                    "In case target list is empty,"
                    " make sure that you have selected the actual\n"
                    "data collection visit (e.g. /dls/i04-1/data/2018/lb18145-70)"
                )
                msgBox.setText(warning)
                start_thread = False

        if start_thread:
            self.work_thread = XChemThread.read_autoprocessing_results_from_disc(
                self.visit_list,
                self.target,
                self.reference_file_list,
                self.database_directory,
                self.data_collection_dict,
                self.preferences,
                self.datasets_summary_file,
                self.initial_model_directory,
                rescore_only,
                self.acceptable_low_resolution_limit_for_data,
                os.path.join(self.database_directory, self.data_source_file),
                self.xce_logfile,
            )
            self.explorer_active = 1
            self.connect(
                self.work_thread,
                QtCore.SIGNAL("update_progress_bar"),
                self.update_progress_bar,
            )
            self.connect(
                self.work_thread,
                QtCore.SIGNAL("update_status_bar(QString)"),
                self.update_status_bar,
            )
            self.connect(
                self.work_thread, QtCore.SIGNAL("finished()"), self.thread_finished
            )
            self.connect(
                self.work_thread,
                QtCore.SIGNAL("create_widgets_for_autoprocessing_results_only"),
                self.create_widgets_for_autoprocessing_results_only,
            )
            self.work_thread.start()

    ####################################################################################
    #
    #
    #
    # => for new module from hell
    # > start

    def update_gdaLog_parsing_instructions_and_score(self, gdaLogInstructions):
        self.gdaLogInstructions = gdaLogInstructions
        self.select_best_autoprocessing_result()

    def read_pinIDs_from_gda_logs(self):
        self.update_log.insert("reading pinIDs from gda logfiles...")
        visit, beamline = XChemMain.getVisitAndBeamline(self.beamline_directory)
        self.work_thread = XChemThread.read_pinIDs_from_gda_logs(
            beamline,
            visit,
            os.path.join(self.database_directory, self.data_source_file),
            self.gdaLogInstructions,
            self.xce_logfile,
        )

        self.explorer_active = 1
        self.connect(
            self.work_thread,
            QtCore.SIGNAL("update_progress_bar"),
            self.update_progress_bar,
        )
        self.connect(
            self.work_thread,
            QtCore.SIGNAL("update_status_bar(QString)"),
            self.update_status_bar,
        )
        self.connect(
            self.work_thread, QtCore.SIGNAL("finished()"), self.thread_finished
        )
        self.connect(
            self.work_thread,
            QtCore.SIGNAL("update_gdaLog_parsing_instructions_and_score"),
            self.update_gdaLog_parsing_instructions_and_score,
        )
        self.work_thread.start()

    def check_for_new_autoprocessing_results(self):
        self.update_log.insert("checking for new data collection")
        if self.target == "=== SELECT TARGET ===":
            self.update_log.error(
                "NO TARGET SELECTED, PLEASE SELECT A TARGET AND TRY AGAIN!"
            )
            start_thread = False
        elif self.target == "=== project directory ===":
            processedDir = self.initial_model_directory
            start_thread = True
        else:
            processedDir = os.path.join(
                self.beamline_directory, "processed", self.target
            )
            start_thread = True

        if start_thread:
            self.work_thread = (
                XChemThread.read_write_autoprocessing_results_from_to_disc(
                    processedDir,
                    os.path.join(self.database_directory, self.data_source_file),
                    self.initial_model_directory,
                    self.xce_logfile,
                    self.target,
                    self.read_agamemnon.isChecked(),
                )
            )

            self.explorer_active = 1
            self.connect(
                self.work_thread,
                QtCore.SIGNAL("update_progress_bar"),
                self.update_progress_bar,
            )
            self.connect(
                self.work_thread,
                QtCore.SIGNAL("update_status_bar(QString)"),
                self.update_status_bar,
            )
            self.connect(
                self.work_thread, QtCore.SIGNAL("finished()"), self.thread_finished
            )
            self.connect(
                self.work_thread,
                QtCore.SIGNAL("read_pinIDs_from_gda_logs"),
                self.read_pinIDs_from_gda_logs,
            )
            self.work_thread.start()

    def select_best_autoprocessing_result(self):
        if self.rescore:
            # first pop up a warning message as this will overwrite all user selections
            msgBox = QtGui.QMessageBox()
            msgBox.setText(
                "*** WARNING ***\n"
                "This will overwrite all your manual selections!\n"
                "Do you want to continue?"
            )
            msgBox.addButton(QtGui.QPushButton("Yes"), QtGui.QMessageBox.YesRole)
            msgBox.addButton(QtGui.QPushButton("No"), QtGui.QMessageBox.RejectRole)
            reply = msgBox.exec_()
            if reply != 0:
                start_thread = False
            else:
                start_thread = True
        else:
            start_thread = True

        if start_thread:
            self.update_log.insert("selecting best autoprocessing result")
            self.update_log.insert(
                "samples where user made manual changes will be ignored!"
            )

            if self.target == "=== project directory ===":
                processedDir = self.initial_model_directory
            else:
                processedDir = os.path.join(
                    self.beamline_directory, "processed", self.target
                )

            visit, beamline = XChemMain.getVisitAndBeamline(processedDir)

            if self.read_agamemnon.isChecked():
                visit = []
                for v in glob.glob(
                    os.path.join(
                        self.beamline_directory[
                            : self.beamline_directory.rfind("-") + 1
                        ]
                        + "*"
                    )
                ):
                    visit.append(v[v.rfind("/") + 1 :])

            self.work_thread = XChemThread.choose_autoprocessing_outcome(
                os.path.join(self.database_directory, self.data_source_file),
                visit,
                self.reference_file_list,
                self.preferences,
                self.initial_model_directory,
                self.rescore,
                self.xce_logfile,
                self.read_agamemnon.isChecked(),
            )

            self.explorer_active = 1
            self.connect(
                self.work_thread,
                QtCore.SIGNAL("update_progress_bar"),
                self.update_progress_bar,
            )
            self.connect(
                self.work_thread,
                QtCore.SIGNAL("update_status_bar(QString)"),
                self.update_status_bar,
            )
            self.connect(
                self.work_thread, QtCore.SIGNAL("finished()"), self.thread_finished
            )
            self.connect(
                self.work_thread,
                QtCore.SIGNAL("populate_datasets_summary_table_NEW"),
                self.populate_datasets_summary_table_NEW,
            )
            self.work_thread.start()

    # < end
    ####################################################################################

    ####################################################################################
    #                                                                                  #
    #                                                 MAPS TAB                         #
    #                                                                                  #
    ####################################################################################

    def set_new_reference_if_applicable(self):
        print("hallo")
        reference_root = str(self.reference_file_selection_combobox.currentText())
        pg_ref = ""
        ucVol_ref = 0.0
        for reference in self.reference_file_list:
            print((reference[0], reference_root))
            if reference[0] == reference_root:
                pg_ref = reference[5]
                ucVol_ref = reference[4]
                break
        if ucVol_ref == 0.0:
            self.update_log.insert(
                "cannot set reference file"
                " since unit cell volume of reference pdb is 0!"
            )
            return

        for xtal in self.initial_model_dimple_dict:
            reference_file_selection_combobox = self.initial_model_dimple_dict[xtal][1]
            self.populate_reference_combobox(reference_file_selection_combobox)
            db_dict = self.xtal_db_dict[xtal]
            pg_xtal = db_dict["DataProcessingPointGroup"]
            ucVol_xtal = db_dict["DataProcessingUnitCellVolume"]

            try:
                difference = math.fabs(1 - (float(ucVol_xtal) / float(ucVol_ref))) * 100
            except ValueError:
                self.update_log.insert(
                    xtal + " -> cannot calculate unit cell volume difference"
                )
                continue

            if (
                pg_xtal == pg_ref
                and difference < self.allowed_unitcell_difference_percent
            ):
                print((xtal, pg_xtal, ucVol_xtal))
                index = reference_file_selection_combobox.findText(
                    reference_root, QtCore.Qt.MatchFixedString
                )
                reference_file_selection_combobox.setCurrentIndex(index)
                self.update_log.insert(
                    xtal
                    + " -> setting "
                    + reference_root
                    + " as input PDB file for DIMPLE"
                )

    def refresh_reference_file_list(self):
        self.reference_file_list = self.get_reference_file_list(" ")
        self.populate_reference_combobox(self.reference_file_selection_combobox)

    def on_context_menu_initial_model(self, point):
        # show context menu
        self.popMenu_for_maps_table.exec_(self.sender().mapToGlobal(point))

    ####################################################################################
    #                                                                                  #
    #                                    PANDDA TAB                                    #
    #                                                                                  #
    ####################################################################################
    def select_pandda_input_template(self):
        mtzin = ""
        filepath_temp = QtGui.QFileDialog.getOpenFileNameAndFilter(
            self.window,
            "Select Example PDB or MTZ File",
            self.initial_model_directory,
            "*.pdb;;*.mtz",
        )
        filepath = str(tuple(filepath_temp)[0])
        pdbin = filepath.split("/")[-1]
        if filepath.endswith(".pdb"):
            pdbin = filepath.split("/")[-1]
            mtzin_temp = pdbin.replace(".pdb", ".mtz")
            if os.path.isfile(filepath.replace(pdbin, mtzin_temp)):
                mtzin = mtzin_temp
            else:
                mtzin = ""
        if filepath.endswith(".mtz"):
            mtzin = filepath.split("/")[-1]
            pdbin_temp = pdbin.replace(".mtz", ".pdb")
            if os.path.isfile(filepath.replace(mtzin, pdbin_temp)):
                pdbin = pdbin_temp
            else:
                pdbin = ""

        try:
            self.pandda_input_data_dir_entry.setText(
                "/"
                + os.path.join(*filepath.split("/")[0 : len(filepath.split("/")) - 2])
            )
        except TypeError:
            self.update_log.error("directory selection invalid")
        self.pandda_pdb_style_entry.setText(pdbin)
        self.pandda_mtz_style_entry.setText(mtzin)

    def change_pandda_spg_label(self):
        combo_text = str(self.pandda_reference_file_selection_combobox.currentText())
        for file in self.reference_file_list:
            if file[0] == combo_text:
                self.pandda_reference_file_spg_label.setText(file[1])
                break

    def on_context_menu_pandda(self, point):
        # show context menu
        self.popMenu_for_pandda_table.exec_(self.sender().mapToGlobal(point))

    ####################################################################################
    #                                                                                  #
    #                                     DEPO TAB                                     #
    #                                                                                  #
    ####################################################################################
    def export_to_html(self):
        XChemWeb.export_to_html(
            self.html_export_directory,
            self.initial_model_directory,
            os.path.join(self.database_directory, self.data_source_file),
            self.xce_logfile,
        ).prepare("0")

    def export_to_html_CompChem(self):
        XChemWeb.export_to_html(
            self.html_export_directory,
            self.initial_model_directory,
            os.path.join(self.database_directory, self.data_source_file),
            self.xce_logfile,
        ).prepare("4")

    def export_to_html_deposition_ready(self):
        XChemWeb.export_to_html(
            self.html_export_directory,
            self.initial_model_directory,
            os.path.join(self.database_directory, self.data_source_file),
            self.xce_logfile,
        ).prepare("5")

    def add_ground_state_db(self):
        pdb, mtz = self.auto_select_ground_state_reference_PDB()
        if pdb is not None:
            db_dict = {
                "DimplePANDDApath": self.panddas_directory,
                "PDB_file": pdb,
                "MTZ_file": mtz,
            }
            self.db.create_or_remove_missing_records_in_depositTable(
                self.xce_logfile, "ground_state", "ground_state", db_dict
            )
        else:
            self.update_log.error(
                "could not find a suitable reference file; see messages above!"
            )

    def auto_select_ground_state_reference_PDB(self):
        pdb = None
        mtz = None
        xtalList = []
        for dirs in glob.glob(
            os.path.join(self.panddas_directory, "processed_datasets", "*")
        ):
            xtal = dirs[dirs.rfind("/") + 1 :]
            if os.path.isfile(os.path.join(dirs, xtal + "-pandda-input.pdb")):
                pdbHeader = parse().PDBheader(
                    os.path.join(dirs, xtal + "-pandda-input.pdb")
                )
                try:
                    xtalList.append(
                        [
                            xtal,
                            float(pdbHeader["Rfree"]),
                            float(pdbHeader["ResolutionHigh"]),
                        ]
                    )
                except ValueError:
                    self.update_log.error(
                        "%s: cannot read Rfree or Resolution from PDB header;"
                        " skipping..."
                    )
                    pass
        self.update_log.insert(
            "found %s PDB files in %s"
            % (
                str(len(xtalList)),
                os.path.join(self.panddas_directory, "processed_datasets"),
            )
        )
        if len(xtalList) >= 10:
            self.update_log.insert(
                "sorting PDBs by Rfree and selecting the 10 with lowest value"
            )
            rfree = sorted(xtalList, key=lambda x: x[1])[:10]
            self.update_log.insert("top 10 PDB files with lowest Rfree:")
            for item in rfree:
                self.update_log.insert(
                    "%s: Rfree = %s | Resolution = %s"
                    % (item[0], str(round(item[1], 3)), str(round(item[2], 2)))
                )
            self.update_log.insert("selecting PDB with highest resolution")
            reso = sorted(rfree, key=lambda x: x[2])[:1]
            self.update_log.insert(
                "selected the following PDB file: %s: Rfree = %s | Resolution = %s"
                % (reso[0][0], str(round(reso[0][1], 3)), str(round(reso[0][2], 2)))
            )
            pdb = os.path.join(
                self.panddas_directory,
                "processed_datasets",
                reso[0][0],
                reso[0][0] + "-pandda-input.pdb",
            )
            mtz = os.path.join(
                self.panddas_directory,
                "processed_datasets",
                reso[0][0],
                reso[0][0] + "-pandda-input.mtz",
            )
        else:
            self.update_log.error(
                "found less than 10 valid PDB files in %s"
                % os.path.join(self.panddas_directory, "processed_datasets")
            )
        return pdb, mtz

    def prepare_ground_state_mmcif(self):
        self.update_log.insert("preparing mmcif file for apo structure deposition")
        self.prepare_models_for_deposition_ligand_bound("ground_state")

    ####################################################################################
    #                                                                                  #
    #                                   SETTINGS TAB                                   #
    #                                                                                  #
    ####################################################################################
    def settings_button_clicked(self):
        if self.sender().text() == "Select Project Directory":
            self.initial_model_directory = str(
                QtGui.QFileDialog.getExistingDirectory(self.window, "Select Directory")
            )
            self.initial_model_directory_label.setText(self.initial_model_directory)
            self.pandda_input_data_dir_entry.setText(self.initial_model_directory)
            self.settings["initial_model_directory"] = self.initial_model_directory
        if self.sender().text() == "Select Reference Structure Directory":
            reference_directory_temp = str(
                QtGui.QFileDialog.getExistingDirectory(self.window, "Select Directory")
            )
            if reference_directory_temp != self.reference_directory:
                self.reference_directory = reference_directory_temp
                self.update_reference_files(" ")
            self.reference_directory_label.setText(self.reference_directory)
            self.settings["reference_directory"] = self.reference_directory
        if self.sender().text() == "Select Data Source File":
            filepath_temp = QtGui.QFileDialog.getOpenFileNameAndFilter(
                self.window, "Select File", self.database_directory, "*.sqlite"
            )
            filepath = str(tuple(filepath_temp)[0])
            self.data_source_file = filepath.split("/")[-1]
            self.database_directory = filepath[: filepath.rfind("/")]
            self.settings["database_directory"] = self.database_directory
            self.settings["data_source"] = os.path.join(
                self.database_directory, self.data_source_file
            )
            write_enabled = self.check_write_permissions_of_data_source()
            if not write_enabled:
                self.data_source_set = False
            else:
                self.data_source_set = True
                self.data_source_file_label.setText(
                    os.path.join(self.database_directory, self.data_source_file)
                )
                self.db = XChemDB.data_source(
                    os.path.join(self.database_directory, self.data_source_file)
                )
                self.db.create_missing_columns()
                self.datasource_menu_reload_samples()
        if self.sender().text() == "Select Data Collection Directory":
            dir_name = str(
                QtGui.QFileDialog.getExistingDirectory(self.window, "Select Directory")
            )
            if dir_name != self.beamline_directory:
                self.beamline_directory = dir_name
                self.target_list, self.visit_list = XChemMain.get_target_and_visit_list(
                    self.beamline_directory, self.read_agamemnon.isChecked()
                )
                self.populate_target_selection_combobox(self.target_selection_combobox)
            self.beamline_directory_label.setText(self.beamline_directory)
            self.settings["beamline_directory"] = self.beamline_directory

        if self.sender().text() == "Select Existing\nCollection Summary File":
            if self.datasets_summary_file != "":
                filepath_temp = QtGui.QFileDialog.getOpenFileNameAndFilter(
                    self.window,
                    "Select File",
                    self.datasets_summary_file[: self.datasets_summary_file.rfind("/")],
                    "*.pkl",
                )
            else:
                filepath_temp = QtGui.QFileDialog.getOpenFileNameAndFilter(
                    self.window, "Select File", os.getcwd(), "*.pkl"
                )
            filepath = str(tuple(filepath_temp)[0])
            self.datasets_summary_file = filepath
            self.datasets_summary_file_label.setText(self.datasets_summary_file)
            self.settings["datasets_summary"] = self.datasets_summary_file

        if self.sender().text() == "Assign New\nCollection Summary File":
            if self.datasets_summary_file != "":
                file_name = str(
                    QtGui.QFileDialog.getSaveFileName(
                        self.window,
                        "New file",
                        self.datasets_summary_file[
                            : self.datasets_summary_file.rfind("/")
                        ],
                    )
                )
            else:
                file_name = str(
                    QtGui.QFileDialog.getSaveFileName(
                        self.window, "New file", self.current_directory
                    )
                )
            # make sure that the file always has .pkl extension
            if str(file_name).rfind(".") != -1:
                file_name = file_name[: file_name.rfind(".")] + ".pkl"
            else:
                file_name = file_name + ".pkl"
            self.datasets_summary_file = file_name
            self.datasets_summary_file_label.setText(self.datasets_summary_file)
            self.settings["datasets_summary"] = self.datasets_summary_file

        if self.sender().text() == "Select CCP4_SCR Directory":
            self.ccp4_scratch_directory = str(
                QtGui.QFileDialog.getExistingDirectory(self.window, "Select Directory")
            )
            self.ccp4_scratch_directory_label.setText(self.ccp4_scratch_directory)
            self.settings["ccp4_scratch"] = self.ccp4_scratch_directory
        if self.sender().text() == "Select PanDDA Directory":
            self.panddas_directory = str(
                QtGui.QFileDialog.getExistingDirectory(self.window, "Select Directory")
            )
            self.panddas_directory_label.setText(self.panddas_directory)
            self.pandda_output_data_dir_entry.setText(self.panddas_directory)
            self.ground_state_pandda_directory_label.setText(self.panddas_directory)
            print(("PANDDA", self.panddas_directory))
            self.settings["panddas_directory"] = self.panddas_directory

            layout_functions.pandda_html(self)

        if self.sender().text() == "Select HTML Export Directory":
            self.html_export_directory = str(
                QtGui.QFileDialog.getExistingDirectory(self.window, "Select Directory")
            )
            self.html_export_directory_label.setText(self.html_export_directory)
            self.settings["html_export_directory"] = self.html_export_directory

        if self.sender().text() == "Select Group deposition Directory":
            self.group_deposit_directory = str(
                QtGui.QFileDialog.getExistingDirectory(self.window, "Select Directory")
            )
            self.group_deposition_directory_label.setText(self.group_deposit_directory)
            self.settings["group_deposit_directory"] = self.group_deposit_directory

    # ############################## sort stuff below here #############################

    def select_sample_for_dimple(self):
        indexes = self.maps_table.selectionModel().selectedRows()
        for index in sorted(indexes):
            xtal = str(self.maps_table.item(index.row(), 0).text())
            self.update_log.insert("{0!s} is marked for DIMPLE".format(index.row()))
            self.initial_model_dimple_dict[xtal][0].setChecked(True)

    def update_summary_plot(self):
        if self.data_source_set:
            XChemPlots.summary_plot(
                os.path.join(self.database_directory, self.data_source_file),
                self.overview_axes,
            ).update_overview()
            self.overview_canvas.draw()

    def show_preferences(self):
        preferences = QtGui.QMessageBox()
        preferencesLayout = preferences.layout()

        vbox = QtGui.QVBoxLayout()
        settings_hbox_filename_root = QtGui.QHBoxLayout()
        filename_root_label = QtGui.QLabel("filename root:")
        settings_hbox_filename_root.addWidget(filename_root_label)
        filename_root_input = QtGui.QLineEdit()
        filename_root_input.setFixedWidth(400)
        filename_root_input.setText(str(self.filename_root))
        filename_root_input.textChanged[str].connect(self.change_filename_root)
        settings_hbox_filename_root.addWidget(filename_root_input)
        vbox.addLayout(settings_hbox_filename_root)

        settings_hbox_adjust_allowed_unit_cell_difference = QtGui.QHBoxLayout()
        adjust_allowed_unit_cell_difference_label = QtGui.QLabel(
            "Max. Allowed Unit Cell Difference between Reference and Target (%):"
        )
        settings_hbox_adjust_allowed_unit_cell_difference.addWidget(
            adjust_allowed_unit_cell_difference_label
        )
        adjust_allowed_unit_cell_difference = QtGui.QLineEdit()
        adjust_allowed_unit_cell_difference.setFixedWidth(200)
        adjust_allowed_unit_cell_difference.setText(
            str(self.allowed_unitcell_difference_percent)
        )
        adjust_allowed_unit_cell_difference.textChanged[str].connect(
            self.change_allowed_unitcell_difference_percent
        )
        settings_hbox_adjust_allowed_unit_cell_difference.addWidget(
            adjust_allowed_unit_cell_difference
        )
        vbox.addLayout(settings_hbox_adjust_allowed_unit_cell_difference)

        settings_hbox_acceptable_low_resolution_limit = QtGui.QHBoxLayout()
        adjust_acceptable_low_resolution_limit_label = QtGui.QLabel(
            "Acceptable low resolution limit for datasets (in Angstrom):"
        )
        settings_hbox_acceptable_low_resolution_limit.addWidget(
            adjust_acceptable_low_resolution_limit_label
        )
        adjust_acceptable_low_resolution_limit = QtGui.QLineEdit()
        adjust_acceptable_low_resolution_limit.setFixedWidth(200)
        adjust_acceptable_low_resolution_limit.setText(
            str(self.acceptable_low_resolution_limit_for_data)
        )
        adjust_acceptable_low_resolution_limit.textChanged[str].connect(
            self.change_acceptable_low_resolution_limit
        )
        settings_hbox_acceptable_low_resolution_limit.addWidget(
            adjust_acceptable_low_resolution_limit
        )
        vbox.addLayout(settings_hbox_acceptable_low_resolution_limit)

        vbox_data = QtGui.QVBoxLayout()
        vbox_data.addWidget(
            QtGui.QLabel(
                "Select amount of processed data you wish to copy to initial_model"
                " directory:"
            )
        )
        self.preferences_data_to_copy_combobox = QtGui.QComboBox()
        for item in self.preferences_data_to_copy:
            self.preferences_data_to_copy_combobox.addItem(item[0])
        self.preferences_data_to_copy_combobox.currentIndexChanged.connect(
            self.preferences_data_to_copy_combobox_changed
        )
        vbox_data.addWidget(self.preferences_data_to_copy_combobox)
        vbox.addLayout(vbox_data)

        vbox_select = QtGui.QVBoxLayout()
        vbox_select.addWidget(QtGui.QLabel("Dataset Selection Mechanism:"))
        self.preferences_selection_mechanism_combobox = QtGui.QComboBox()
        for item in self.preferences_selection_mechanism:
            self.preferences_selection_mechanism_combobox.addItem(item)
        self.preferences_selection_mechanism_combobox.currentIndexChanged.connect(
            self.preferences_selection_mechanism_combobox_changed
        )
        index = self.preferences_selection_mechanism_combobox.findText(
            self.preferences["dataset_selection_mechanism"], QtCore.Qt.MatchFixedString
        )
        self.preferences_selection_mechanism_combobox.setCurrentIndex(index)
        vbox_select.addWidget(self.preferences_selection_mechanism_combobox)
        vbox.addLayout(vbox_select)

        vbox_restraints = QtGui.QVBoxLayout()
        vbox_restraints.addWidget(QtGui.QLabel("Restraints generation program:"))
        self.preferences_restraints_generation_combobox = QtGui.QComboBox()
        program_list = []

        if self.external_software["acedrg"]:
            program_list.append("acedrg")
        if self.external_software["phenix.elbow"]:
            program_list.append("phenix.elbow")
        if self.external_software["grade"]:
            program_list.append("grade")
        for item in program_list:
            self.preferences_restraints_generation_combobox.addItem(item)
        self.preferences_restraints_generation_combobox.currentIndexChanged.connect(
            self.preferences_restraints_generation_combobox_changed
        )
        index = self.preferences_restraints_generation_combobox.findText(
            self.restraints_program, QtCore.Qt.MatchFixedString
        )
        self.preferences_restraints_generation_combobox.setCurrentIndex(index)
        vbox_restraints.addWidget(self.preferences_restraints_generation_combobox)
        vbox.addLayout(vbox_restraints)

        hbox = QtGui.QHBoxLayout()
        hbox.addWidget(QtGui.QLabel("XCE logfile:"))
        self.xce_logfile_label = QtGui.QLabel(self.xce_logfile)
        hbox.addWidget(self.xce_logfile_label)
        button = QtGui.QPushButton("Change")
        button.clicked.connect(self.set_xce_logfile)
        hbox.addWidget(button)
        vbox.addLayout(hbox)

        settings_hbox_max_queue_jobs = QtGui.QHBoxLayout()
        adjust_max_queue_jobs_label = QtGui.QLabel(
            "Max. number of jobs running at once on DLS cluster:"
        )
        settings_hbox_max_queue_jobs.addWidget(adjust_max_queue_jobs_label)
        adjust_max_queue_jobs = QtGui.QLineEdit()
        adjust_max_queue_jobs.setFixedWidth(200)
        adjust_max_queue_jobs.setText(str(self.max_queue_jobs))
        adjust_max_queue_jobs.textChanged[str].connect(self.change_max_queue_jobs)
        settings_hbox_max_queue_jobs.addWidget(adjust_max_queue_jobs)
        vbox.addLayout(settings_hbox_max_queue_jobs)

        settings_hbox_dimple_twin_mode = QtGui.QHBoxLayout()
        self.dimple_twin_mode_label_checkbox = QtGui.QCheckBox(
            "run DIMPLE in TWIN mode"
        )
        if self.preferences["dimple_twin_mode"]:
            self.dimple_twin_mode_label_checkbox.setChecked(True)
        self.dimple_twin_mode_label_checkbox.toggled.connect(
            self.dimple_change_twin_mode
        )
        settings_hbox_dimple_twin_mode.addWidget(self.dimple_twin_mode_label_checkbox)
        vbox.addLayout(settings_hbox_dimple_twin_mode)

        hbox = QtGui.QHBoxLayout()
        hbox.addWidget(QtGui.QLabel("Additional CIF file for non-standard ligand:"))
        self.second_cif_file_label = QtGui.QLabel(self.second_cif_file)
        hbox.addWidget(self.second_cif_file_label)
        button = QtGui.QPushButton("Select")
        button.clicked.connect(self.set_second_cif_file)
        hbox.addWidget(button)
        vbox.addLayout(hbox)

        preferencesLayout.addLayout(vbox, 0, 0)

        preferences.exec_()

    def dimple_change_twin_mode(self):
        if self.preferences["dimple_twin_mode"]:
            self.update_log.insert(
                "changing preferences: turning off DIMPLE in TWIN mode"
            )
            self.preferences["dimple_twin_mode"] = False
        else:
            self.update_log.insert("changing preferences: changing DIMPLE to TWIN mode")
            self.preferences["dimple_twin_mode"] = True

    def enter_pdb_codes(self):
        pdbID_entry = QtGui.QMessageBox()
        pdbID_entryLayout = pdbID_entry.layout()

        vbox = QtGui.QVBoxLayout()

        frame = QtGui.QFrame()
        frame.setFrameShape(QtGui.QFrame.StyledPanel)

        grid = QtGui.QGridLayout()

        grid.addWidget(QtGui.QLabel("Text from PDB email"), 0, 0)
        self.pdb_code_entry = QtGui.QTextEdit()
        self.pdb_code_entry.setText("")
        self.pdb_code_entry.setFixedWidth(500)
        grid.addWidget(self.pdb_code_entry, 1, 0, 20, 1)

        frame.setLayout(grid)
        vbox.addWidget(frame)

        hbox = QtGui.QHBoxLayout()
        button = QtGui.QPushButton("Update Database")
        button.clicked.connect(self.update_database_with_pdb_codes)
        hbox.addWidget(button)

        vbox.addLayout(hbox)
        pdbID_entryLayout.addLayout(vbox, 0, 0)
        pdbID_entry.exec_()

    def add_label_information(self):
        label_entry = QtGui.QMessageBox()
        label_entryLayout = label_entry.layout()

        try:
            labelInfo = self.db.get_label_info_from_db()
        except AttributeError:
            self.update_log.warning("please specify DB file first")
            return None

        vbox = QtGui.QVBoxLayout()

        frame = QtGui.QFrame()
        frame.setFrameShape(QtGui.QFrame.StyledPanel)

        grid = QtGui.QGridLayout()
        grid.addWidget(QtGui.QLabel("label"), 0, 0)
        grid.addWidget(QtGui.QLabel("description"), 0, 1)

        self.labelList = []
        for i in range(5):
            labelEdit = QtGui.QLineEdit()
            descriptionEdit = QtGui.QLineEdit()
            grid.addWidget(labelEdit, i + 1, 0)
            grid.addWidget(descriptionEdit, i + 1, 1)
            try:
                labelEdit.setText(labelInfo[i][0])
                descriptionEdit.setText(labelInfo[i][1])
            except IndexError:
                labelEdit.setText("")
                descriptionEdit.setText("")
            labelEdit.setFixedWidth(100)
            descriptionEdit.setFixedWidth(500)
            self.labelList.append([labelEdit, descriptionEdit])
        frame.setLayout(grid)
        vbox.addWidget(frame)

        hbox = QtGui.QHBoxLayout()
        button = QtGui.QPushButton("Update Database")
        button.clicked.connect(self.update_database_with_labelInfo)
        hbox.addWidget(button)

        vbox.addLayout(hbox)
        label_entryLayout.addLayout(vbox, 0, 0)
        label_entry.exec_()

    def prepare_models_for_deposition_ligand_bound(self, structureType):
        start_thread = True
        self.update_log.insert("preparing mmcif files for PDB group deposition...")
        ignore_event_map = False
        if structureType == "ground_state":
            try:
                self.update_log.insert("ground-state deposition")
                data_template_dict = self.db.get_deposit_dict_for_sample("ground_state")
                pdb = data_template_dict["PDB_file"]
                self.update_log.insert("looking for ground-state PDB: " + pdb)
                if not os.path.isfile(pdb):
                    self.update_log.error(
                        "ground-state PDB does not exist; stopping..."
                    )
                    start_thread = False
                mtz = data_template_dict["MTZ_file"]
                self.update_log.insert("looking for ground-state MTZ: " + mtz)
                if not os.path.isfile(mtz):
                    self.update_log.error(
                        "ground-state MTZ does not exist; stopping..."
                    )
                    start_thread = False
                ground_state = [pdb, mtz, self.panddas_directory]
            except KeyError:
                self.update_log.error(
                    "seems like there is no entry for ground-state in database"
                )
                start_thread = False
        else:
            ground_state = []
            if self.deposition_bounnd_state_preparation_ignore_event_map.isChecked():
                ignore_event_map = True

        if start_thread:
            if ground_state != []:
                self.update_log.insert("apo PDB: " + ground_state[0])
                self.update_log.insert("apo MTZ: " + ground_state[1])
                self.update_log.insert("pandda directory: " + ground_state[2])
            overwrite_existing_mmcif = True
            self.work_thread = XChemDeposit.prepare_mmcif_files_for_deposition(
                os.path.join(self.database_directory, self.data_source_file),
                self.xce_logfile,
                overwrite_existing_mmcif,
                self.initial_model_directory,
                ground_state,
                ignore_event_map,
            )
            self.explorer_active = 1
            self.connect(
                self.work_thread,
                QtCore.SIGNAL("update_progress_bar"),
                self.update_progress_bar,
            )
            self.connect(
                self.work_thread,
                QtCore.SIGNAL("update_status_bar(QString)"),
                self.update_status_bar,
            )
            self.connect(
                self.work_thread, QtCore.SIGNAL("finished()"), self.thread_finished
            )
            self.work_thread.start()

    def prepare_for_group_deposition_upload_ligand_bound(self):
        self.work_thread = XChemDeposit.prepare_for_group_deposition_upload(
            os.path.join(self.database_directory, self.data_source_file),
            self.xce_logfile,
            self.group_deposit_directory,
            self.initial_model_directory,
            "ligand_bound",
        )
        self.explorer_active = 1
        self.connect(
            self.work_thread,
            QtCore.SIGNAL("update_progress_bar"),
            self.update_progress_bar,
        )
        self.connect(
            self.work_thread,
            QtCore.SIGNAL("update_status_bar(QString)"),
            self.update_status_bar,
        )
        self.connect(
            self.work_thread, QtCore.SIGNAL("finished()"), self.thread_finished
        )
        self.work_thread.start()

    def prepare_for_group_deposition_upload_ground_state(self):
        self.work_thread = XChemDeposit.prepare_for_group_deposition_upload(
            os.path.join(self.database_directory, self.data_source_file),
            self.xce_logfile,
            self.group_deposit_directory,
            self.initial_model_directory,
            "ground_state",
        )
        self.explorer_active = 1
        self.connect(
            self.work_thread,
            QtCore.SIGNAL("update_progress_bar"),
            self.update_progress_bar,
        )
        self.connect(
            self.work_thread,
            QtCore.SIGNAL("update_status_bar(QString)"),
            self.update_status_bar,
        )
        self.connect(
            self.work_thread, QtCore.SIGNAL("finished()"), self.thread_finished
        )
        self.work_thread.start()

    def check_smiles_in_db_and_pdb(self):
        self.work_thread = XChemDeposit.compare_smiles_in_db_with_ligand_in_pdb(
            self.initial_model_directory,
            os.path.join(self.database_directory, self.data_source_file),
            self.xce_logfile,
        )
        self.explorer_active = 1
        self.connect(
            self.work_thread,
            QtCore.SIGNAL("update_progress_bar"),
            self.update_progress_bar,
        )
        self.connect(
            self.work_thread,
            QtCore.SIGNAL("update_status_bar(QString)"),
            self.update_status_bar,
        )
        self.connect(
            self.work_thread, QtCore.SIGNAL("finished()"), self.thread_finished
        )
        self.connect(
            self.work_thread, QtCore.SIGNAL("show_error_dict"), self.show_error_dict
        )
        self.work_thread.start()

    def deposition_data(self):
        depositData = QtGui.QMessageBox()
        depositDataLayout = depositData.layout()

        vbox = QtGui.QVBoxLayout()

        deposit_tab_widget = QtGui.QTabWidget()
        deposit_tab_list = [
            "Contact",
            "General",
            "Authors",
            "Citation",
            "Molecule",
            "Misc",
            "Methods",
            "Software",
            "Funding",
        ]

        deposit_tab_dict = {}
        for page in deposit_tab_list:
            tab = QtGui.QWidget()
            vb = QtGui.QVBoxLayout(tab)
            deposit_tab_widget.addTab(tab, page)
            deposit_tab_dict[page] = [tab, vb]

        # PI and scientist info
        vb = QtGui.QVBoxLayout()
        hbox = QtGui.QHBoxLayout()

        frame = QtGui.QFrame()
        frame.setFrameShape(QtGui.QFrame.StyledPanel)

        grid = QtGui.QGridLayout()
        grid.addWidget(QtGui.QLabel("Principal Investigator"), 0, 0)

        grid.addWidget(QtGui.QLabel("Salutation"), 1, 0)
        self.contact_author_PI_salutation = QtGui.QLineEdit()
        self.contact_author_PI_salutation.setText("Dr.")
        self.contact_author_PI_salutation.setFixedWidth(200)
        grid.addWidget(self.contact_author_PI_salutation, 1, 1)

        grid.addWidget(QtGui.QLabel("First name"), 2, 0)
        self.contact_author_PI_first_name = QtGui.QLineEdit()
        self.contact_author_PI_first_name.setText("")
        self.contact_author_PI_first_name.setFixedWidth(200)
        grid.addWidget(self.contact_author_PI_first_name, 2, 1)

        grid.addWidget(QtGui.QLabel("Last name"), 3, 0)
        self.contact_author_PI_last_name = QtGui.QLineEdit()
        self.contact_author_PI_last_name.setText("")
        self.contact_author_PI_last_name.setFixedWidth(200)
        grid.addWidget(self.contact_author_PI_last_name, 3, 1)

        grid.addWidget(QtGui.QLabel("Middle name"), 4, 0)
        self.contact_author_PI_middle_name = QtGui.QLineEdit()
        self.contact_author_PI_middle_name.setText("")
        self.contact_author_PI_middle_name.setFixedWidth(200)
        self.contact_author_PI_middle_name.setStyleSheet(
            "background-color: rgb(192, 192, 192);"
        )
        grid.addWidget(self.contact_author_PI_middle_name, 4, 1)

        grid.addWidget(QtGui.QLabel("PI role"), 5, 0)
        self.contact_author_PI_role = QtGui.QComboBox()
        PIroles = ["principal investigator/group leader"]
        for item in PIroles:
            self.contact_author_PI_role.addItem(item)
        grid.addWidget(self.contact_author_PI_role, 5, 1)

        grid.addWidget(QtGui.QLabel("Organization type"), 6, 0)
        self.contact_author_PI_organization_type = QtGui.QComboBox()
        Organizations = ["academic", "commercial", "government"]
        for item in Organizations:
            self.contact_author_PI_organization_type.addItem(item)
        grid.addWidget(self.contact_author_PI_organization_type, 6, 1)

        grid.addWidget(QtGui.QLabel("Organization Name"), 7, 0)
        self.contact_author_PI_organization_name = QtGui.QLineEdit()
        self.contact_author_PI_organization_name.setText("")
        self.contact_author_PI_organization_name.setFixedWidth(200)
        grid.addWidget(self.contact_author_PI_organization_name, 7, 1)

        grid.addWidget(QtGui.QLabel("Email"), 8, 0)
        self.contact_author_PI_email = QtGui.QLineEdit()
        self.contact_author_PI_email.setText("")
        self.contact_author_PI_email.setFixedWidth(200)
        grid.addWidget(self.contact_author_PI_email, 8, 1)

        grid.addWidget(QtGui.QLabel("Street"), 9, 0)
        self.contact_author_PI_address = QtGui.QLineEdit()
        self.contact_author_PI_address.setText("")
        self.contact_author_PI_address.setFixedWidth(200)
        grid.addWidget(self.contact_author_PI_address, 9, 1)

        grid.addWidget(QtGui.QLabel("City"), 10, 0)
        self.contact_author_PI_city = QtGui.QLineEdit()
        self.contact_author_PI_city.setText("")
        self.contact_author_PI_city.setFixedWidth(200)
        grid.addWidget(self.contact_author_PI_city, 10, 1)

        grid.addWidget(QtGui.QLabel("State"), 11, 0)
        self.contact_author_PI_State_or_Province = QtGui.QLineEdit()
        self.contact_author_PI_State_or_Province.setText("")
        self.contact_author_PI_State_or_Province.setFixedWidth(200)
        self.contact_author_PI_State_or_Province.setStyleSheet(
            "background-color: rgb(192, 192, 192);"
        )
        grid.addWidget(self.contact_author_PI_State_or_Province, 11, 1)

        grid.addWidget(QtGui.QLabel("ZIP code"), 12, 0)
        self.contact_author_PI_Zip_Code = QtGui.QLineEdit()
        self.contact_author_PI_Zip_Code.setText("")
        self.contact_author_PI_Zip_Code.setFixedWidth(200)
        grid.addWidget(self.contact_author_PI_Zip_Code, 12, 1)

        grid.addWidget(QtGui.QLabel("Country"), 13, 0)
        self.contact_author_PI_Country = QtGui.QLineEdit()
        self.contact_author_PI_Country.setText("")
        self.contact_author_PI_Country.setFixedWidth(200)
        grid.addWidget(self.contact_author_PI_Country, 13, 1)

        grid.addWidget(QtGui.QLabel("Phone"), 14, 0)
        self.contact_author_PI_phone_number = QtGui.QLineEdit()
        self.contact_author_PI_phone_number.setText("")
        self.contact_author_PI_phone_number.setFixedWidth(200)
        grid.addWidget(self.contact_author_PI_phone_number, 14, 1)

        grid.addWidget(QtGui.QLabel("ORCID"), 15, 0)
        self.contact_author_PI_ORCID = QtGui.QLineEdit()
        self.contact_author_PI_ORCID.setText("")
        self.contact_author_PI_ORCID.setFixedWidth(200)
        grid.addWidget(self.contact_author_PI_ORCID, 15, 1)

        frame.setLayout(grid)
        hbox.addWidget(frame)

        frame = QtGui.QFrame()
        frame.setFrameShape(QtGui.QFrame.StyledPanel)
        grid = QtGui.QGridLayout()
        grid.addWidget(QtGui.QLabel("Responsible Scientist"), 0, 0)

        grid.addWidget(QtGui.QLabel("Salutation"), 1, 0)
        self.contact_author_salutation = QtGui.QLineEdit()
        self.contact_author_salutation.setText("Dr.")
        self.contact_author_salutation.setFixedWidth(200)
        grid.addWidget(self.contact_author_salutation, 1, 1)

        grid.addWidget(QtGui.QLabel("First name"), 2, 0)
        self.contact_author_first_name = QtGui.QLineEdit()
        self.contact_author_first_name.setText("")
        self.contact_author_first_name.setFixedWidth(200)
        grid.addWidget(self.contact_author_first_name, 2, 1)

        grid.addWidget(QtGui.QLabel("Last name"), 3, 0)
        self.contact_author_last_name = QtGui.QLineEdit()
        self.contact_author_last_name.setText("")
        self.contact_author_last_name.setFixedWidth(200)
        grid.addWidget(self.contact_author_last_name, 3, 1)

        grid.addWidget(QtGui.QLabel("Middle name"), 4, 0)
        self.contact_author_middle_name = QtGui.QLineEdit()
        self.contact_author_middle_name.setText("")
        self.contact_author_middle_name.setFixedWidth(200)
        self.contact_author_middle_name.setStyleSheet(
            "background-color: rgb(192, 192, 192);"
        )
        grid.addWidget(self.contact_author_middle_name, 4, 1)

        grid.addWidget(QtGui.QLabel("Role"), 5, 0)

        self.contact_author_role = QtGui.QComboBox()
        ScientistRoles = ["responsible scientist", "investigator"]
        for item in ScientistRoles:
            self.contact_author_role.addItem(item)
        grid.addWidget(self.contact_author_role, 5, 1)

        grid.addWidget(QtGui.QLabel("Organization type"), 6, 0)

        self.contact_author_organization_type = QtGui.QComboBox()
        for item in Organizations:
            self.contact_author_organization_type.addItem(item)
        grid.addWidget(self.contact_author_organization_type, 6, 1)

        grid.addWidget(QtGui.QLabel("Organization Name"), 7, 0)
        self.contact_author_organization_name = QtGui.QLineEdit()
        self.contact_author_organization_name.setText("")
        self.contact_author_organization_name.setFixedWidth(200)
        grid.addWidget(self.contact_author_organization_name, 7, 1)

        grid.addWidget(QtGui.QLabel("Email"), 8, 0)
        self.contact_author_email = QtGui.QLineEdit()
        self.contact_author_email.setText("")
        self.contact_author_email.setFixedWidth(200)
        grid.addWidget(self.contact_author_email, 8, 1)

        grid.addWidget(QtGui.QLabel("Street"), 9, 0)
        self.contact_author_address = QtGui.QLineEdit()
        self.contact_author_address.setText("")
        self.contact_author_address.setFixedWidth(200)
        grid.addWidget(self.contact_author_address, 9, 1)

        grid.addWidget(QtGui.QLabel("City"), 10, 0)
        self.contact_author_city = QtGui.QLineEdit()
        self.contact_author_city.setText("")
        self.contact_author_city.setFixedWidth(200)
        grid.addWidget(self.contact_author_city, 10, 1)

        grid.addWidget(QtGui.QLabel("State"), 11, 0)
        self.contact_author_State_or_Province = QtGui.QLineEdit()
        self.contact_author_State_or_Province.setText("")
        self.contact_author_State_or_Province.setFixedWidth(200)
        self.contact_author_State_or_Province.setStyleSheet(
            "background-color: rgb(192, 192, 192);"
        )
        grid.addWidget(self.contact_author_State_or_Province, 11, 1)

        grid.addWidget(QtGui.QLabel("ZIP code"), 12, 0)
        self.contact_author_Zip_Code = QtGui.QLineEdit()
        self.contact_author_Zip_Code.setText("")
        self.contact_author_Zip_Code.setFixedWidth(200)
        grid.addWidget(self.contact_author_Zip_Code, 12, 1)

        grid.addWidget(QtGui.QLabel("Country"), 13, 0)
        self.contact_author_Country = QtGui.QLineEdit()
        self.contact_author_Country.setText("")
        self.contact_author_Country.setFixedWidth(200)
        grid.addWidget(self.contact_author_Country, 13, 1)

        grid.addWidget(QtGui.QLabel("Phone"), 14, 0)
        self.contact_author_phone_number = QtGui.QLineEdit()
        self.contact_author_phone_number.setText("")
        self.contact_author_phone_number.setFixedWidth(200)
        grid.addWidget(self.contact_author_phone_number, 14, 1)

        grid.addWidget(QtGui.QLabel("ORCID"), 15, 0)
        self.contact_author_ORCID = QtGui.QLineEdit()
        self.contact_author_ORCID.setText("")
        self.contact_author_ORCID.setFixedWidth(200)
        grid.addWidget(self.contact_author_ORCID, 15, 1)

        frame.setLayout(grid)
        hbox.addWidget(frame)

        vb.addLayout(hbox)
        vb.addWidget(QtGui.QLabel(XChemToolTips.deposition_interface_note()))
        vb.addStretch(1)

        deposit_tab_dict["Contact"][1].addLayout(vb)

        # release status
        vb = QtGui.QVBoxLayout()

        frame = QtGui.QFrame()
        frame.setFrameShape(QtGui.QFrame.StyledPanel)

        grid = QtGui.QGridLayout()
        grid.addWidget(QtGui.QLabel("Release status"), 0, 0)

        grid.addWidget(QtGui.QLabel("Release Status for sequence"), 4, 0)

        self.Release_status_for_sequence = QtGui.QComboBox()
        codeStatus = ["RELEASE NOW", "HOLD FOR RELEASE"]
        for item in codeStatus:
            self.Release_status_for_sequence.addItem(item)
        grid.addWidget(self.Release_status_for_sequence, 4, 1)

        grid.addWidget(QtGui.QLabel("Release Status for coordinates/ SF"), 8, 0)
        self.Release_status_for_coordinates = QtGui.QComboBox()
        coordStatus = [
            "RELEASE NOW",
            "HOLD FOR PUBLICATION",
            "HOLD FOR 4 WEEKS",
            "HOLD FOR 6 MONTHS",
            "HOLD FOR 1 YEAR",
        ]
        for item in coordStatus:
            self.Release_status_for_coordinates.addItem(item)
        grid.addWidget(self.Release_status_for_coordinates, 8, 1)

        frame.setLayout(grid)
        vb.addWidget(frame)

        frame = QtGui.QFrame()
        frame.setFrameShape(QtGui.QFrame.StyledPanel)

        grid = QtGui.QGridLayout()
        grid.addWidget(QtGui.QLabel("Title & Details"), 0, 0)
        note = (
            "Note: supported wildcards: $ProteinName,$CompoundName;"
            ' e.g. "Crystal Structure of human JMJD2D in complex with N2317a"'
        )
        grid.addWidget(QtGui.QLabel(note), 1, 0)

        grid.addWidget(QtGui.QLabel("Group deposition title"), 2, 0)
        self.group_deposition_title = QtGui.QLineEdit()
        self.group_deposition_title.setText("PanDDA analysis group deposition")
        self.group_deposition_title.setFixedWidth(600)

        grid.addWidget(self.group_deposition_title, 2, 1)

        grid.addWidget(QtGui.QLabel("Description"), 3, 0)
        self.group_description = QtGui.QLineEdit()
        self.group_description.setText(
            "XDomainX of XOrganismX $ProteinName screened against the"
            " XXX Fragment Library by X-ray Crystallography at the XChem facility of"
            " Diamond Light Source beamline I04-1"
        )
        self.group_description.setFixedWidth(600)
        grid.addWidget(self.group_description, 3, 1)

        grid.addWidget(QtGui.QLabel("Structure Title (ligand bound)"), 4, 0)
        self.structure_title = QtGui.QLineEdit()
        self.structure_title.setText(
            "Crystal Structure of $ProteinName in complex with $CompoundName"
        )
        self.structure_title.setFixedWidth(600)
        grid.addWidget(self.structure_title, 4, 1)

        note = "\n\nApo Structure:\nonly use if you want to deposit PanDDA models!"
        grid.addWidget(QtGui.QLabel(note), 6, 0)

        grid.addWidget(QtGui.QLabel("Structure Title (apo)"), 7, 0)
        self.structure_title_apo = QtGui.QLineEdit()
        self.structure_title_apo.setText(
            "PanDDA analysis group deposition of ground-state model of $ProteinName"
        )
        self.structure_title_apo.setFixedWidth(600)
        grid.addWidget(self.structure_title_apo, 7, 1)

        frame.setLayout(grid)
        vb.addWidget(frame)

        vb.addStretch(1)

        deposit_tab_dict["General"][1].addLayout(vb)

        # authors
        vb = QtGui.QVBoxLayout()

        frame = QtGui.QFrame()
        frame.setFrameShape(QtGui.QFrame.StyledPanel)

        grid = QtGui.QGridLayout()
        grid.addWidget(QtGui.QLabel("Deposition authors (e.g. Surname, F.M.)"), 0, 0)

        self.structure_author_name_List = []

        for column in range(0, 2):
            for row in range(1, 15):
                structure_author_name = QtGui.QLineEdit()
                structure_author_name.setText("")
                structure_author_name.setFixedWidth(300)
                grid.addWidget(structure_author_name, row, column)
                self.structure_author_name_List.append(structure_author_name)

        frame.setLayout(grid)
        vb.addWidget(frame)

        vb.addStretch(1)

        deposit_tab_dict["Authors"][1].addLayout(vb)

        # primary citation
        vb = QtGui.QVBoxLayout()

        frame = QtGui.QFrame()
        frame.setFrameShape(QtGui.QFrame.StyledPanel)

        grid = QtGui.QGridLayout()
        grid.addWidget(QtGui.QLabel("Primary Citation"), 0, 0)

        grid.addWidget(QtGui.QLabel("ID"), 1, 0)
        self.primary_citation_id = QtGui.QLineEdit()
        self.primary_citation_id.setText("primary")
        self.primary_citation_id.setFixedWidth(500)
        grid.addWidget(self.primary_citation_id, 1, 1)

        grid.addWidget(QtGui.QLabel("Journal"), 2, 0)
        self.primary_citation_journal_abbrev = QtGui.QLineEdit()
        self.primary_citation_journal_abbrev.setText("To be published")
        self.primary_citation_journal_abbrev.setFixedWidth(500)
        grid.addWidget(self.primary_citation_journal_abbrev, 2, 1)

        grid.addWidget(QtGui.QLabel("Title"), 3, 0)
        self.primary_citation_title = QtGui.QLineEdit()
        self.primary_citation_title.setText("")
        self.primary_citation_title.setFixedWidth(500)
        self.primary_citation_title.setStyleSheet(
            "background-color: rgb(192, 192, 192);"
        )
        grid.addWidget(self.primary_citation_title, 3, 1)

        grid.addWidget(QtGui.QLabel("Year"), 4, 0)
        self.primary_citation_year = QtGui.QLineEdit()
        self.primary_citation_year.setText("")
        self.primary_citation_year.setFixedWidth(500)
        self.primary_citation_year.setStyleSheet(
            "background-color: rgb(192, 192, 192);"
        )
        grid.addWidget(self.primary_citation_year, 4, 1)

        grid.addWidget(QtGui.QLabel("Volume"), 5, 0)
        self.primary_citation_journal_volume = QtGui.QLineEdit()
        self.primary_citation_journal_volume.setText("")
        self.primary_citation_journal_volume.setFixedWidth(500)
        self.primary_citation_journal_volume.setStyleSheet(
            "background-color: rgb(192, 192, 192);"
        )
        grid.addWidget(self.primary_citation_journal_volume, 5, 1)

        grid.addWidget(QtGui.QLabel("Page, first"), 6, 0)
        self.primary_citation_page_first = QtGui.QLineEdit()
        self.primary_citation_page_first.setText("")
        self.primary_citation_page_first.setFixedWidth(500)
        self.primary_citation_page_first.setStyleSheet(
            "background-color: rgb(192, 192, 192);"
        )
        grid.addWidget(self.primary_citation_page_first, 6, 1)

        grid.addWidget(QtGui.QLabel("Page, last"), 7, 0)
        self.primary_citation_page_last = QtGui.QLineEdit()
        self.primary_citation_page_last.setText("")
        self.primary_citation_page_last.setFixedWidth(500)
        self.primary_citation_page_last.setStyleSheet(
            "background-color: rgb(192, 192, 192);"
        )
        grid.addWidget(self.primary_citation_page_last, 7, 1)

        frame.setLayout(grid)
        vb.addWidget(frame)

        # citation authors
        frame = QtGui.QFrame()
        frame.setFrameShape(QtGui.QFrame.StyledPanel)

        grid = QtGui.QGridLayout()
        self.set_primary_citation_authors = QtGui.QCheckBox(
            "same as deposition authors"
        )
        layout_functions.add_checkbox(
            self,
            self.set_primary_citation_authors,
            "xce_object.set_primary_citation_as_structure_authors",
        )
        grid.addWidget(self.set_primary_citation_authors, 0, 0)

        self.primary_citation_author_name_List = []

        for column in range(0, 2):
            for row in range(1, 15):
                primary_citation_author_name = QtGui.QLineEdit()
                primary_citation_author_name.setText("")
                primary_citation_author_name.setFixedWidth(300)
                grid.addWidget(primary_citation_author_name, row, column)
                self.primary_citation_author_name_List.append(
                    primary_citation_author_name
                )

        frame.setLayout(grid)
        vb.addWidget(frame)

        vb.addStretch(1)

        deposit_tab_dict["Citation"][1].addLayout(vb)

        # molecule info
        vb = QtGui.QVBoxLayout()

        frame = QtGui.QFrame()
        frame.setFrameShape(QtGui.QFrame.StyledPanel)

        grid = QtGui.QGridLayout()

        grid.addWidget(QtGui.QLabel("Entity 1"), 1, 0)

        grid.addWidget(QtGui.QLabel("Molecule Name"), 2, 0)
        self.molecule_name = QtGui.QLineEdit()
        self.molecule_name.setText("")
        self.molecule_name.setFixedWidth(300)

        grid.addWidget(self.molecule_name, 2, 1)
        grid.addWidget(QtGui.QLabel("(e.g. RNA Hammerhead Ribozyme)"), 2, 2)

        grid.addWidget(QtGui.QLabel("Fragment Name"), 3, 0)
        self.fragment_name_one = QtGui.QLineEdit()
        self.fragment_name_one.setText("")
        self.fragment_name_one.setFixedWidth(300)
        self.fragment_name_one.setStyleSheet("background-color: rgb(192, 192, 192);")
        grid.addWidget(self.fragment_name_one, 3, 1)
        grid.addWidget(QtGui.QLabel("(e.g. ligand binding domain, hairpin)"), 3, 2)

        grid.addWidget(QtGui.QLabel("Specific Mutation"), 4, 0)
        self.fragment_name_one_specific_mutation = QtGui.QLineEdit()
        self.fragment_name_one_specific_mutation.setText("")
        self.fragment_name_one_specific_mutation.setFixedWidth(300)
        self.fragment_name_one_specific_mutation.setStyleSheet(
            "background-color: rgb(192, 192, 192);"
        )
        grid.addWidget(self.fragment_name_one_specific_mutation, 4, 1)
        grid.addWidget(QtGui.QLabel("(e.g. C280S)"), 4, 2)

        grid.addWidget(QtGui.QLabel("Enzyme Comission Number"), 5, 0)
        self.fragment_name_one_enzyme_comission_number = QtGui.QLineEdit()
        self.fragment_name_one_enzyme_comission_number.setText("")
        self.fragment_name_one_enzyme_comission_number.setFixedWidth(300)
        self.fragment_name_one_enzyme_comission_number.setStyleSheet(
            "background-color: rgb(192, 192, 192);"
        )
        grid.addWidget(self.fragment_name_one_enzyme_comission_number, 5, 1)
        grid.addWidget(QtGui.QLabel("(if known: e.g. 2.7.7.7)"), 5, 2)

        grid.addWidget(QtGui.QLabel("Genetically Manipulated Source"), 6, 0)

        grid.addWidget(QtGui.QLabel("Source organism scientific name"), 7, 0)

        self.Source_organism_scientific_name = QtGui.QComboBox()
        taxonomy_dict = XChemMain.NCBI_taxonomy_ID()
        for item in taxonomy_dict:
            self.Source_organism_scientific_name.addItem(taxonomy_dict[item])
        grid.addWidget(self.Source_organism_scientific_name, 7, 1)

        grid.addWidget(QtGui.QLabel("Source organism gene"), 8, 0)
        self.Source_organism_gene = QtGui.QLineEdit()
        self.Source_organism_gene.setText("")
        self.Source_organism_gene.setFixedWidth(300)
        grid.addWidget(self.Source_organism_gene, 8, 1)
        grid.addWidget(QtGui.QLabel("(e.g. RPOD, ALKA...)"), 8, 2)

        grid.addWidget(QtGui.QLabel("Source organism strain"), 9, 0)
        self.Source_organism_strain = QtGui.QLineEdit()
        self.Source_organism_strain.setText("")
        self.Source_organism_strain.setFixedWidth(300)
        self.Source_organism_strain.setStyleSheet(
            "background-color: rgb(192, 192, 192);"
        )
        grid.addWidget(self.Source_organism_strain, 9, 1)
        grid.addWidget(QtGui.QLabel("(e.g. BH10 ISOLATE, K-12...)"), 9, 2)

        grid.addWidget(QtGui.QLabel("Expression system scientific name"), 10, 0)

        self.Expression_system_scientific_name = QtGui.QComboBox()
        for item in taxonomy_dict:
            self.Expression_system_scientific_name.addItem(taxonomy_dict[item])
        grid.addWidget(self.Expression_system_scientific_name, 10, 1)

        grid.addWidget(QtGui.QLabel("Expression system strain"), 11, 0)
        self.Expression_system_strain = QtGui.QLineEdit()
        self.Expression_system_strain.setText("")
        self.Expression_system_strain.setFixedWidth(300)
        self.Expression_system_strain.setStyleSheet(
            "background-color: rgb(192, 192, 192);"
        )
        grid.addWidget(self.Expression_system_strain, 11, 1)
        grid.addWidget(QtGui.QLabel("(e.g. BL21(DE3))"), 11, 2)

        grid.addWidget(QtGui.QLabel("Expression system vector type"), 12, 0)
        self.Expression_system_vector_type = QtGui.QLineEdit()
        self.Expression_system_vector_type.setText("")
        self.Expression_system_vector_type.setFixedWidth(300)
        self.Expression_system_vector_type.setStyleSheet(
            "background-color: rgb(192, 192, 192);"
        )
        grid.addWidget(self.Expression_system_vector_type, 12, 1)
        grid.addWidget(QtGui.QLabel("(e.g. plasmid)"), 12, 2)

        grid.addWidget(QtGui.QLabel("Expression_system_plasmid_name"), 13, 0)
        self.Expression_system_plasmid_name = QtGui.QLineEdit()
        self.Expression_system_plasmid_name.setText("")
        self.Expression_system_plasmid_name.setFixedWidth(300)
        self.Expression_system_plasmid_name.setStyleSheet(
            "background-color: rgb(192, 192, 192);"
        )
        grid.addWidget(self.Expression_system_plasmid_name, 13, 1)
        grid.addWidget(QtGui.QLabel("(e.g. pET26)"), 13, 2)

        grid.addWidget(QtGui.QLabel("Manipulated_source_details"), 14, 0)
        self.Manipulated_source_details = QtGui.QLineEdit()
        self.Manipulated_source_details.setText("")
        self.Manipulated_source_details.setFixedWidth(300)
        self.Manipulated_source_details.setStyleSheet(
            "background-color: rgb(192, 192, 192);"
        )
        grid.addWidget(self.Manipulated_source_details, 14, 1)
        grid.addWidget(QtGui.QLabel("(any other relevant information)"), 14, 2)

        grid.addWidget(QtGui.QLabel("Chains"), 15, 0)
        self.molecule_chain_one = QtGui.QLineEdit()
        self.molecule_chain_one.setText("")
        self.molecule_chain_one.setFixedWidth(300)
        grid.addWidget(self.molecule_chain_one, 15, 1)
        grid.addWidget(QtGui.QLabel("(e.g. A or A,B)"), 15, 2)

        frame.setLayout(grid)
        vb.addWidget(frame)

        # entity 2

        frame = QtGui.QFrame()
        frame.setFrameShape(QtGui.QFrame.StyledPanel)

        grid = QtGui.QGridLayout()

        grid.addWidget(
            QtGui.QLabel(
                "Entity 2 (IMPORTANT: only fill in if you are working with a"
                " protein-protein complex!)"
            ),
            1,
            0,
        )

        grid.addWidget(QtGui.QLabel("Molecule Name"), 2, 0)
        self.molecule_name_two = QtGui.QLineEdit()
        self.molecule_name_two.setText("")
        self.molecule_name_two.setFixedWidth(300)

        grid.addWidget(self.molecule_name_two, 2, 1)
        grid.addWidget(QtGui.QLabel("(e.g. RNA Hammerhead Ribozyme)"), 2, 2)

        grid.addWidget(QtGui.QLabel("Fragment Name"), 3, 0)
        self.fragment_name_two = QtGui.QLineEdit()
        self.fragment_name_two.setText("")
        self.fragment_name_two.setFixedWidth(300)
        self.fragment_name_two.setStyleSheet("background-color: rgb(192, 192, 192);")
        grid.addWidget(self.fragment_name_two, 3, 1)
        grid.addWidget(QtGui.QLabel("(e.g. ligand binding domain, hairpin)"), 3, 2)

        grid.addWidget(QtGui.QLabel("Specific Mutation"), 4, 0)
        self.fragment_name_two_specific_mutation = QtGui.QLineEdit()
        self.fragment_name_two_specific_mutation.setText("")
        self.fragment_name_two_specific_mutation.setFixedWidth(300)
        self.fragment_name_two_specific_mutation.setStyleSheet(
            "background-color: rgb(192, 192, 192);"
        )
        grid.addWidget(self.fragment_name_two_specific_mutation, 4, 1)
        grid.addWidget(QtGui.QLabel("(e.g. C280S)"), 4, 2)

        grid.addWidget(QtGui.QLabel("Enzyme Comission Number"), 5, 0)
        self.fragment_name_two_enzyme_comission_number = QtGui.QLineEdit()
        self.fragment_name_two_enzyme_comission_number.setText("")
        self.fragment_name_two_enzyme_comission_number.setFixedWidth(300)
        self.fragment_name_two_enzyme_comission_number.setStyleSheet(
            "background-color: rgb(192, 192, 192);"
        )
        grid.addWidget(self.fragment_name_two_enzyme_comission_number, 5, 1)
        grid.addWidget(QtGui.QLabel("(if known: e.g. 2.7.7.7)"), 5, 2)

        grid.addWidget(QtGui.QLabel("Genetically Manipulated Source"), 6, 0)

        grid.addWidget(QtGui.QLabel("Source organism scientific name"), 7, 0)

        self.Source_organism_scientific_name_two = QtGui.QComboBox()
        taxonomy_dict = XChemMain.NCBI_taxonomy_ID()
        for item in taxonomy_dict:
            self.Source_organism_scientific_name_two.addItem(taxonomy_dict[item])
        grid.addWidget(self.Source_organism_scientific_name_two, 7, 1)

        grid.addWidget(QtGui.QLabel("Source organism gene"), 8, 0)
        self.Source_organism_gene_two = QtGui.QLineEdit()
        self.Source_organism_gene_two.setText("")
        self.Source_organism_gene_two.setFixedWidth(300)
        grid.addWidget(self.Source_organism_gene_two, 8, 1)
        grid.addWidget(QtGui.QLabel("(e.g. RPOD, ALKA...)"), 8, 2)

        grid.addWidget(QtGui.QLabel("Source organism strain"), 9, 0)
        self.Source_organism_strain_two = QtGui.QLineEdit()
        self.Source_organism_strain_two.setText("")
        self.Source_organism_strain_two.setFixedWidth(300)
        self.Source_organism_strain_two.setStyleSheet(
            "background-color: rgb(192, 192, 192);"
        )
        grid.addWidget(self.Source_organism_strain_two, 9, 1)
        grid.addWidget(QtGui.QLabel("(e.g. BH10 ISOLATE, K-12...)"), 9, 2)

        grid.addWidget(QtGui.QLabel("Expression system scientific name"), 10, 0)

        self.Expression_system_scientific_name_two = QtGui.QComboBox()
        for item in taxonomy_dict:
            self.Expression_system_scientific_name_two.addItem(taxonomy_dict[item])
        grid.addWidget(self.Expression_system_scientific_name_two, 10, 1)

        grid.addWidget(QtGui.QLabel("Expression system strain"), 11, 0)
        self.Expression_system_strain_two = QtGui.QLineEdit()
        self.Expression_system_strain_two.setText("")
        self.Expression_system_strain_two.setFixedWidth(300)
        self.Expression_system_strain_two.setStyleSheet(
            "background-color: rgb(192, 192, 192);"
        )
        grid.addWidget(self.Expression_system_strain_two, 11, 1)
        grid.addWidget(QtGui.QLabel("(e.g. BL21(DE3))"), 11, 2)

        grid.addWidget(QtGui.QLabel("Expression system vector type"), 12, 0)
        self.Expression_system_vector_type_two = QtGui.QLineEdit()
        self.Expression_system_vector_type_two.setText("")
        self.Expression_system_vector_type_two.setFixedWidth(300)
        self.Expression_system_vector_type_two.setStyleSheet(
            "background-color: rgb(192, 192, 192);"
        )
        grid.addWidget(self.Expression_system_vector_type_two, 12, 1)
        grid.addWidget(QtGui.QLabel("(e.g. plasmid)"), 12, 2)

        grid.addWidget(QtGui.QLabel("Expression_system_plasmid_name"), 13, 0)
        self.Expression_system_plasmid_name_two = QtGui.QLineEdit()
        self.Expression_system_plasmid_name_two.setText("")
        self.Expression_system_plasmid_name_two.setFixedWidth(300)
        self.Expression_system_plasmid_name_two.setStyleSheet(
            "background-color: rgb(192, 192, 192);"
        )
        grid.addWidget(self.Expression_system_plasmid_name_two, 13, 1)
        grid.addWidget(QtGui.QLabel("(e.g. pET26)"), 13, 2)

        grid.addWidget(QtGui.QLabel("Manipulated_source_details"), 14, 0)
        self.Manipulated_source_details_two = QtGui.QLineEdit()
        self.Manipulated_source_details_two.setText("")
        self.Manipulated_source_details_two.setFixedWidth(300)
        self.Manipulated_source_details_two.setStyleSheet(
            "background-color: rgb(192, 192, 192);"
        )
        grid.addWidget(self.Manipulated_source_details_two, 14, 1)
        grid.addWidget(QtGui.QLabel("(any other relevant information)"), 14, 2)

        grid.addWidget(QtGui.QLabel("Chains"), 15, 0)
        self.molecule_chain_two = QtGui.QLineEdit()
        self.molecule_chain_two.setText("")
        self.molecule_chain_two.setFixedWidth(300)
        grid.addWidget(self.molecule_chain_two, 15, 1)
        grid.addWidget(QtGui.QLabel("(e.g. A or A,B)"), 15, 2)

        frame.setLayout(grid)

        vb.addWidget(frame)

        # entity 2 --- END

        vb.addStretch(1)

        deposit_tab_dict["Molecule"][1].addLayout(vb)

        # misc
        vb = QtGui.QVBoxLayout()

        frame = QtGui.QFrame()
        frame.setFrameShape(QtGui.QFrame.StyledPanel)

        grid = QtGui.QGridLayout()

        grid.addWidget(QtGui.QLabel("Keywords"), 1, 0)
        self.structure_keywords = QtGui.QLineEdit()
        self.structure_keywords.setText(
            "SGC - Diamond I04-1 fragment screening, PanDDA, XChemExplorer"
        )
        self.structure_keywords.setFixedWidth(300)
        grid.addWidget(self.structure_keywords, 1, 1)
        grid.addWidget(QtGui.QLabel("(e.g. beta barrel, protein-DNA complex)"), 1, 2)

        grid.addWidget(QtGui.QLabel("Type"), 2, 0)
        self.structure_keywords_type = QtGui.QComboBox()
        self.structure_keywords_type.setStyleSheet(
            "background-color: rgb(192, 192, 192);"
        )
        for item in XChemMain.pdbx_keywords():
            self.structure_keywords_type.addItem(item)
        grid.addWidget(self.structure_keywords_type, 2, 1)

        grid.addWidget(QtGui.QLabel("Biological Assembly"), 3, 0)
        self.biological_assembly_chain_number = QtGui.QLineEdit()
        self.biological_assembly_chain_number.setText("")
        self.biological_assembly_chain_number.setFixedWidth(300)
        grid.addWidget(self.biological_assembly_chain_number, 3, 1)
        grid.addWidget(QtGui.QLabel("(e.g.  1 for monomer, 2 for dimer ..)"), 3, 2)

        grid.addWidget(QtGui.QLabel("Sequence UNIPROT ID"), 4, 0)
        self.molecule_one_letter_sequence_uniprot_id = QtGui.QLineEdit()
        self.molecule_one_letter_sequence_uniprot_id.setText("")
        self.molecule_one_letter_sequence_uniprot_id.setFixedWidth(300)
        grid.addWidget(self.molecule_one_letter_sequence_uniprot_id, 4, 1)
        grid.addWidget(QtGui.QLabel("(e.g.  Q6B0I6)"), 4, 2)

        grid.addWidget(QtGui.QLabel("Sequence"), 5, 0)
        self.molecule_one_letter_sequence = QtGui.QTextEdit()
        self.molecule_one_letter_sequence.setStyleSheet(
            "background-color: rgb(255, 255, 255);"
        )

        self.molecule_one_letter_sequence.setText("")
        self.molecule_one_letter_sequence.setFixedWidth(300)
        grid.addWidget(self.molecule_one_letter_sequence, 5, 1, 8, 2)

        grid.addWidget(QtGui.QLabel("Sequence UNIPROT ID (Entity 2) - optional"), 13, 0)
        self.molecule_one_letter_sequence_uniprot_id_two = QtGui.QLineEdit()
        self.molecule_one_letter_sequence_uniprot_id_two.setText("")
        self.molecule_one_letter_sequence_uniprot_id_two.setStyleSheet(
            "background-color: rgb(192, 192, 192);"
        )
        self.molecule_one_letter_sequence_uniprot_id_two.setFixedWidth(300)
        grid.addWidget(self.molecule_one_letter_sequence_uniprot_id_two, 13, 1)
        grid.addWidget(QtGui.QLabel("(e.g.  Q6B0I6)"), 13, 2)

        grid.addWidget(QtGui.QLabel("Sequence (Entity 2) - optional"), 14, 0)
        self.molecule_one_letter_sequence_two = QtGui.QTextEdit()
        self.molecule_one_letter_sequence_two.setText("")
        self.molecule_one_letter_sequence_two.setFixedWidth(300)
        grid.addWidget(self.molecule_one_letter_sequence_two, 14, 1, 19, 2)

        grid.addWidget(QtGui.QLabel("Structural Genomic (optional)"), 21, 0)

        grid.addWidget(QtGui.QLabel("Project Name"), 22, 0)
        self.SG_project_name = QtGui.QLineEdit()
        self.SG_project_name.setText("")
        self.SG_project_name.setStyleSheet("background-color: rgb(192, 192, 192);")
        self.SG_project_name.setFixedWidth(300)
        grid.addWidget(self.SG_project_name, 22, 1)
        grid.addWidget(
            QtGui.QLabel("(e.g. SGC, Structural Genomics Consortium)"), 22, 2
        )

        grid.addWidget(QtGui.QLabel("Full Name"), 23, 0)
        self.full_name_of_SG_center = QtGui.QLineEdit()
        self.full_name_of_SG_center.setText("")
        self.full_name_of_SG_center.setStyleSheet(
            "background-color: rgb(192, 192, 192);"
        )
        self.full_name_of_SG_center.setFixedWidth(300)
        grid.addWidget(self.full_name_of_SG_center, 23, 1)
        grid.addWidget(QtGui.QLabel("(e.g. Structural Genomics Consortium)"), 23, 2)

        frame.setLayout(grid)
        vb.addWidget(frame)

        vb.addStretch(1)

        deposit_tab_dict["Misc"][1].addLayout(vb)

        # methods
        vb = QtGui.QVBoxLayout()

        frame = QtGui.QFrame()
        frame.setFrameShape(QtGui.QFrame.StyledPanel)

        grid = QtGui.QGridLayout()

        grid.addWidget(QtGui.QLabel("Crystallization"), 1, 0)

        grid.addWidget(QtGui.QLabel("Method"), 2, 0)

        self.crystallization_method = QtGui.QComboBox()
        for item in XChemMain.crystal_growth_methods():
            self.crystallization_method.addItem(item)
        grid.addWidget(self.crystallization_method, 2, 1)

        grid.addWidget(QtGui.QLabel("pH"), 3, 0)
        self.crystallization_pH = QtGui.QLineEdit()
        self.crystallization_pH.setText("")
        self.crystallization_pH.setFixedWidth(300)
        grid.addWidget(self.crystallization_pH, 3, 1)
        grid.addWidget(QtGui.QLabel("(e.g. 7.5 ...)"), 3, 2)

        grid.addWidget(QtGui.QLabel("Temperature"), 4, 0)
        self.crystallization_temperature = QtGui.QLineEdit()
        self.crystallization_temperature.setText("")
        self.crystallization_temperature.setFixedWidth(300)
        grid.addWidget(self.crystallization_temperature, 4, 1)
        grid.addWidget(QtGui.QLabel("(e.g. 298) (in Kelvin)"), 4, 2)

        grid.addWidget(QtGui.QLabel("Condition"), 5, 0)
        self.crystallization_details = QtGui.QLineEdit()
        self.crystallization_details.setText("")
        self.crystallization_details.setFixedWidth(300)
        grid.addWidget(self.crystallization_details, 5, 1)
        grid.addWidget(QtGui.QLabel("(e.g. PEG 4000, NaCl etc.)"), 5, 2)

        grid.addWidget(QtGui.QLabel("Diffraction Experiment"), 6, 0)
        note = (
            "Note: this information will only be used if it is\n"
            "not already available in the mainTable!\n"
            "Ignore if data were collected at DLS"
        )
        grid.addWidget(QtGui.QLabel(note), 7, 0)

        grid.addWidget(QtGui.QLabel("Source"), 8, 0)

        self.radiation_source = QtGui.QComboBox()
        for item in XChemMain.radiationSource():
            self.radiation_source.addItem(item)
        grid.addWidget(self.radiation_source, 8, 1)

        grid.addWidget(QtGui.QLabel("Source Type"), 9, 0)

        self.radiation_source_type = QtGui.QComboBox()
        for item in XChemMain.wwBeamlines():
            self.radiation_source_type.addItem(item)
        grid.addWidget(self.radiation_source_type, 9, 1)

        grid.addWidget(QtGui.QLabel("Wavelength"), 10, 0)
        self.radiation_wavelengths = QtGui.QLineEdit()
        self.radiation_wavelengths.setText("")
        self.radiation_wavelengths.setFixedWidth(300)
        grid.addWidget(self.radiation_wavelengths, 10, 1)
        grid.addWidget(QtGui.QLabel("(e.g. 1.502)"), 10, 2)

        grid.addWidget(QtGui.QLabel("Detector"), 11, 0)

        self.radiation_detector = QtGui.QComboBox()
        for item in XChemMain.detector():
            self.radiation_detector.addItem(item)
        grid.addWidget(self.radiation_detector, 11, 1)

        grid.addWidget(QtGui.QLabel("Detector Type"), 12, 0)

        self.radiation_detector_type = QtGui.QComboBox()
        for item in XChemMain.detectorType():
            self.radiation_detector_type.addItem(item)
        grid.addWidget(self.radiation_detector_type, 12, 1)

        grid.addWidget(QtGui.QLabel("Date"), 13, 0)
        self.data_collection_date = QtGui.QLineEdit()
        self.data_collection_date.setText("")
        self.data_collection_date.setFixedWidth(300)
        grid.addWidget(self.data_collection_date, 13, 1)
        grid.addWidget(QtGui.QLabel("(e.g. 2004-01-07)"), 13, 2)

        grid.addWidget(QtGui.QLabel("Temperature"), 14, 0)
        self.data_collection_temperature = QtGui.QLineEdit()
        self.data_collection_temperature.setText("")
        self.data_collection_temperature.setFixedWidth(300)
        grid.addWidget(self.data_collection_temperature, 14, 1)
        grid.addWidget(QtGui.QLabel("(e.g. 100) (in Kelvin)"), 14, 2)

        grid.addWidget(QtGui.QLabel("Protocol"), 15, 0)
        self.data_collection_protocol = QtGui.QLineEdit()
        self.data_collection_protocol.setText("SINGLE WAVELENGTH")
        self.data_collection_protocol.setFixedWidth(300)
        grid.addWidget(self.data_collection_protocol, 15, 1)
        grid.addWidget(QtGui.QLabel("(e.g. SINGLE WAVELENGTH, MAD, ...)"), 15, 2)

        frame.setLayout(grid)
        vb.addWidget(frame)

        vb.addStretch(1)

        deposit_tab_dict["Methods"][1].addLayout(vb)

        # software
        vb = QtGui.QVBoxLayout()

        frame = QtGui.QFrame()
        frame.setFrameShape(QtGui.QFrame.StyledPanel)

        grid = QtGui.QGridLayout()

        grid.addWidget(QtGui.QLabel("PDB starting model"), 1, 0)
        self.pdbx_starting_model = QtGui.QLineEdit()
        self.pdbx_starting_model.setText("")
        self.pdbx_starting_model.setFixedWidth(300)
        grid.addWidget(self.pdbx_starting_model, 1, 1)
        grid.addWidget(QtGui.QLabel("(e.g. 7.5 ...)"), 1, 2)

        grid.addWidget(QtGui.QLabel("Data reduction"), 2, 0)
        self.data_integration_software = QtGui.QComboBox()
        for item in XChemMain.data_integration_software():
            self.data_integration_software.addItem(item)
        grid.addWidget(self.data_integration_software, 2, 1)

        grid.addWidget(QtGui.QLabel("Phasing"), 3, 0)
        self.phasing_software = QtGui.QComboBox()
        for item in XChemMain.phasing_software():
            self.phasing_software.addItem(item)
        grid.addWidget(self.phasing_software, 3, 1)

        frame.setLayout(grid)
        vb.addWidget(frame)
        vb.addStretch(1)

        deposit_tab_dict["Software"][1].addLayout(vb)

        # Funding

        vb = QtGui.QVBoxLayout()

        frame = QtGui.QFrame()
        frame.setFrameShape(QtGui.QFrame.StyledPanel)

        grid = QtGui.QGridLayout()

        grid.addWidget(QtGui.QLabel("Funding Organization"), 1, 0)
        self.pdbx_funding_organization_one = QtGui.QLineEdit()
        self.pdbx_funding_organization_one.setText("")
        self.pdbx_funding_organization_one.setFixedWidth(700)
        grid.addWidget(self.pdbx_funding_organization_one, 1, 1)

        grid.addWidget(QtGui.QLabel("Grant Number"), 2, 0)
        self.pdbx_grant_number_one = QtGui.QLineEdit()
        self.pdbx_grant_number_one.setText("")
        self.pdbx_grant_number_one.setFixedWidth(700)
        grid.addWidget(self.pdbx_grant_number_one, 2, 1)

        grid.addWidget(QtGui.QLabel("Country"), 3, 0)
        self.pdbx_grant_country_one = QtGui.QComboBox()
        for item in XChemMain.pdbx_country():
            self.pdbx_grant_country_one.addItem(item)
        grid.addWidget(self.pdbx_grant_country_one, 3, 1)

        frame.setLayout(grid)
        vb.addWidget(frame)

        frame = QtGui.QFrame()
        frame.setFrameShape(QtGui.QFrame.StyledPanel)

        grid = QtGui.QGridLayout()

        grid.addWidget(QtGui.QLabel("Funding Organization"), 1, 0)
        self.pdbx_funding_organization_two = QtGui.QLineEdit()
        self.pdbx_funding_organization_two.setText("")
        self.pdbx_funding_organization_two.setFixedWidth(700)
        grid.addWidget(self.pdbx_funding_organization_two, 1, 1)

        grid.addWidget(QtGui.QLabel("Grant Number"), 2, 0)
        self.pdbx_grant_number_two = QtGui.QLineEdit()
        self.pdbx_grant_number_two.setText("")
        self.pdbx_grant_number_two.setFixedWidth(700)
        grid.addWidget(self.pdbx_grant_number_two, 2, 1)

        grid.addWidget(QtGui.QLabel("Country"), 3, 0)
        self.pdbx_grant_country_two = QtGui.QComboBox()
        for item in XChemMain.pdbx_country():
            self.pdbx_grant_country_two.addItem(item)
        grid.addWidget(self.pdbx_grant_country_two, 3, 1)

        frame.setLayout(grid)
        vb.addWidget(frame)

        frame = QtGui.QFrame()
        frame.setFrameShape(QtGui.QFrame.StyledPanel)

        grid = QtGui.QGridLayout()

        grid.addWidget(QtGui.QLabel("Funding Organization"), 1, 0)
        self.pdbx_funding_organization_three = QtGui.QLineEdit()
        self.pdbx_funding_organization_three.setText("")
        self.pdbx_funding_organization_three.setFixedWidth(700)
        grid.addWidget(self.pdbx_funding_organization_three, 1, 1)

        grid.addWidget(QtGui.QLabel("Grant Number"), 2, 0)
        self.pdbx_grant_number_three = QtGui.QLineEdit()
        self.pdbx_grant_number_three.setText("")
        self.pdbx_grant_number_three.setFixedWidth(700)
        grid.addWidget(self.pdbx_grant_number_three, 2, 1)

        grid.addWidget(QtGui.QLabel("Country"), 3, 0)
        self.pdbx_grant_country_three = QtGui.QComboBox()
        for item in XChemMain.pdbx_country():
            self.pdbx_grant_country_three.addItem(item)
        grid.addWidget(self.pdbx_grant_country_three, 3, 1)

        frame.setLayout(grid)
        vb.addWidget(frame)

        vb.addStretch(1)

        deposit_tab_dict["Funding"][1].addLayout(vb)

        vbox.addWidget(deposit_tab_widget)

        hbox = QtGui.QHBoxLayout()
        button = QtGui.QPushButton("Load\nFile")
        button.clicked.connect(self.load_deposit_config_file)
        hbox.addWidget(button)
        button = QtGui.QPushButton("Save\nFile")
        button.clicked.connect(self.save_deposit_config_file)
        hbox.addWidget(button)
        button = QtGui.QPushButton("Load from\nDatabase")
        button.clicked.connect(self.load_deposit_from_database)
        button.setEnabled(False)
        hbox.addWidget(button)
        button = QtGui.QPushButton("Save to\nDatabase")
        button.clicked.connect(self.save_deposit_to_database)
        hbox.addWidget(button)

        vbox.addLayout(hbox)
        depositDataLayout.addLayout(vbox, 0, 0)

        depositData.exec_()

    def save_deposit_config_file(self):
        self.update_deposit_dict()
        file_name = str(
            QtGui.QFileDialog.getSaveFileName(
                self.window, "Save file", self.current_directory
            )
        )
        # make sure that the file always has .deposit extension
        if str(file_name).rfind(".") != -1:
            file_name = file_name[: file_name.rfind(".")] + ".deposit"
        else:
            file_name = file_name + ".deposit"
        pickle.dump(self.deposit_dict, open(file_name, "wb"))

    def update_database_with_pdb_codes(self):
        self.work_thread = XChemDeposit.import_PDB_IDs(
            str(self.pdb_code_entry.toPlainText()),
            os.path.join(self.database_directory, self.data_source_file),
            self.xce_logfile,
        )
        self.explorer_active = 1
        self.connect(
            self.work_thread,
            QtCore.SIGNAL("update_progress_bar"),
            self.update_progress_bar,
        )
        self.connect(
            self.work_thread,
            QtCore.SIGNAL("update_status_bar(QString)"),
            self.update_status_bar,
        )
        self.connect(
            self.work_thread, QtCore.SIGNAL("finished()"), self.thread_finished
        )
        self.work_thread.start()

    def update_database_with_labelInfo(self):
        for n, l in enumerate(self.labelList):
            label = str(l[0].text())
            description = str(l[1].text())

            self.db.execute_statement(
                "update labelTable set Label='%s',Description='%s' where ID=%s"
                % (label, description, str(n + 1))
            )

    def load_deposit_config_file(self):
        file_name_temp = QtGui.QFileDialog.getOpenFileNameAndFilter(
            self.window, "Open file", self.current_directory, "*.deposit"
        )
        file_name = tuple(file_name_temp)[0]
        self.deposit_dict = pickle.load(open(file_name, "rb"))

        for key in self.get_deposit_dict_template():
            if key not in self.deposit_dict:
                self.update_log.warning("field not in .deposit file: " + str(key))
                self.deposit_dict[key] = ""
        self.update_deposit_input()

    def load_deposit_from_database(self):
        print("hallo")

    def save_deposit_to_database(self):
        self.update_deposit_dict()
        msgBox = QtGui.QMessageBox()
        msgBox.setText(
            "*** WARNING ***\n"
            "Are you sure you want to update the database?\n"
            "This will overwrite previous entries!"
        )
        msgBox.addButton(QtGui.QPushButton("Yes"), QtGui.QMessageBox.YesRole)
        msgBox.addButton(QtGui.QPushButton("No"), QtGui.QMessageBox.RejectRole)
        reply = msgBox.exec_()
        if reply == 0:
            self.work_thread = XChemDeposit.update_depositTable(
                self.deposit_dict,
                os.path.join(self.database_directory, self.data_source_file),
                self.xce_logfile,
            )
            self.explorer_active = 1
            self.connect(
                self.work_thread,
                QtCore.SIGNAL("update_progress_bar"),
                self.update_progress_bar,
            )
            self.connect(
                self.work_thread,
                QtCore.SIGNAL("update_status_bar(QString)"),
                self.update_status_bar,
            )
            self.connect(
                self.work_thread, QtCore.SIGNAL("finished()"), self.thread_finished
            )
            self.work_thread.start()

    def update_deposit_input(self):
        try:
            self.contact_author_PI_salutation.setText(
                self.deposit_dict["contact_author_PI_salutation"]
            )
            self.contact_author_PI_first_name.setText(
                self.deposit_dict["contact_author_PI_first_name"]
            )
            self.contact_author_PI_last_name.setText(
                self.deposit_dict["contact_author_PI_last_name"]
            )
            self.contact_author_PI_middle_name.setText(
                self.deposit_dict["contact_author_PI_middle_name"]
            )
            index = self.contact_author_PI_role.findText(
                self.deposit_dict["contact_author_PI_role"], QtCore.Qt.MatchFixedString
            )
            self.contact_author_PI_role.setCurrentIndex(index)
            index = self.contact_author_PI_organization_type.findText(
                self.deposit_dict["contact_author_PI_organization_type"],
                QtCore.Qt.MatchFixedString,
            )
            self.contact_author_PI_organization_type.setCurrentIndex(index)
            self.contact_author_PI_organization_name.setText(
                self.deposit_dict["contact_author_PI_organization_name"]
            )
            self.contact_author_PI_email.setText(
                self.deposit_dict["contact_author_PI_email"]
            )
            self.contact_author_PI_address.setText(
                self.deposit_dict["contact_author_PI_address"]
            )
            self.contact_author_PI_city.setText(
                self.deposit_dict["contact_author_PI_city"]
            )
            self.contact_author_PI_State_or_Province.setText(
                self.deposit_dict["contact_author_PI_State_or_Province"]
            )
            self.contact_author_PI_Zip_Code.setText(
                self.deposit_dict["contact_author_PI_Zip_Code"]
            )
            self.contact_author_PI_Country.setText(
                self.deposit_dict["contact_author_PI_Country"]
            )
            self.contact_author_PI_phone_number.setText(
                self.deposit_dict["contact_author_PI_phone_number"]
            )
            self.contact_author_PI_ORCID.setText(
                self.deposit_dict["contact_author_PI_ORCID"]
            )

            self.contact_author_salutation.setText(
                self.deposit_dict["contact_author_salutation"]
            )
            self.contact_author_first_name.setText(
                self.deposit_dict["contact_author_first_name"]
            )
            self.contact_author_last_name.setText(
                self.deposit_dict["contact_author_last_name"]
            )
            self.contact_author_middle_name.setText(
                self.deposit_dict["contact_author_middle_name"]
            )
            index = self.contact_author_role.findText(
                self.deposit_dict["contact_author_role"], QtCore.Qt.MatchFixedString
            )
            self.contact_author_role.setCurrentIndex(index)
            index = self.contact_author_organization_type.findText(
                self.deposit_dict["contact_author_organization_type"],
                QtCore.Qt.MatchFixedString,
            )
            self.contact_author_organization_type.setCurrentIndex(index)
            self.contact_author_organization_name.setText(
                self.deposit_dict["contact_author_organization_name"]
            )
            self.contact_author_email.setText(self.deposit_dict["contact_author_email"])
            self.contact_author_address.setText(
                self.deposit_dict["contact_author_address"]
            )
            self.contact_author_city.setText(self.deposit_dict["contact_author_city"])
            self.contact_author_State_or_Province.setText(
                self.deposit_dict["contact_author_State_or_Province"]
            )
            self.contact_author_Zip_Code.setText(
                self.deposit_dict["contact_author_Zip_Code"]
            )
            self.contact_author_Country.setText(
                self.deposit_dict["contact_author_Country"]
            )
            self.contact_author_phone_number.setText(
                self.deposit_dict["contact_author_phone_number"]
            )
            self.contact_author_ORCID.setText(self.deposit_dict["contact_author_ORCID"])
            index = self.Release_status_for_coordinates.findText(
                self.deposit_dict["Release_status_for_coordinates"],
                QtCore.Qt.MatchFixedString,
            )
            self.Release_status_for_coordinates.setCurrentIndex(index)
            index = self.Release_status_for_sequence.findText(
                self.deposit_dict["Release_status_for_sequence"],
                QtCore.Qt.MatchFixedString,
            )
            self.Release_status_for_sequence.setCurrentIndex(index)

            self.group_deposition_title.setText(
                self.deposit_dict["group_deposition_title"]
            )
            self.group_description.setText(self.deposit_dict["group_description"])

            self.structure_title.setText(self.deposit_dict["structure_title"])
            self.structure_title_apo.setText(self.deposit_dict["structure_title_apo"])

            for n, name in enumerate(
                self.deposit_dict["structure_author_name"].split(";")
            ):
                self.structure_author_name_List[n].setText(name)

            self.primary_citation_id.setText(self.deposit_dict["primary_citation_id"])
            self.primary_citation_journal_abbrev.setText(
                self.deposit_dict["primary_citation_journal_abbrev"]
            )
            self.primary_citation_title.setText(
                self.deposit_dict["primary_citation_title"]
            )
            self.primary_citation_year.setText(
                self.deposit_dict["primary_citation_year"]
            )
            self.primary_citation_journal_volume.setText(
                self.deposit_dict["primary_citation_journal_volume"]
            )
            self.primary_citation_page_first.setText(
                self.deposit_dict["primary_citation_page_first"]
            )
            self.primary_citation_page_last.setText(
                self.deposit_dict["primary_citation_page_last"]
            )

            for n, name in enumerate(
                self.deposit_dict["primary_citation_author_name"].split(";")
            ):
                self.primary_citation_author_name_List[n].setText(name)

            # entity 1

            self.molecule_name.setText(self.deposit_dict["molecule_name"])
            self.fragment_name_one_specific_mutation.setText(
                self.deposit_dict["fragment_name_one_specific_mutation"]
            )
            index = self.Source_organism_scientific_name.findText(
                self.deposit_dict["Source_organism_scientific_name"],
                QtCore.Qt.MatchFixedString,
            )
            self.Source_organism_scientific_name.setCurrentIndex(index)

            self.Source_organism_gene.setText(self.deposit_dict["Source_organism_gene"])
            self.Source_organism_strain.setText(
                self.deposit_dict["Source_organism_strain"]
            )
            index = self.Expression_system_scientific_name.findText(
                self.deposit_dict["Expression_system_scientific_name"],
                QtCore.Qt.MatchFixedString,
            )
            self.Expression_system_scientific_name.setCurrentIndex(index)

            self.Expression_system_strain.setText(
                self.deposit_dict["Expression_system_strain"]
            )
            self.Expression_system_vector_type.setText(
                self.deposit_dict["Expression_system_vector_type"]
            )
            self.Expression_system_plasmid_name.setText(
                self.deposit_dict["Expression_system_plasmid_name"]
            )
            self.Manipulated_source_details.setText(
                self.deposit_dict["Manipulated_source_details"]
            )

            self.molecule_chain_one.setText(self.deposit_dict["molecule_chain_one"])
            # entity 2
            self.molecule_name_two.setText(self.deposit_dict["molecule_name_two"])
            self.fragment_name_two_specific_mutation.setText(
                self.deposit_dict["fragment_name_two_specific_mutation"]
            )
            index = self.Source_organism_scientific_name_two.findText(
                self.deposit_dict["Source_organism_scientific_name_two"],
                QtCore.Qt.MatchFixedString,
            )
            self.Source_organism_scientific_name_two.setCurrentIndex(index)
            self.Source_organism_gene_two.setText(
                self.deposit_dict["Source_organism_gene_two"]
            )
            self.Source_organism_strain_two.setText(
                self.deposit_dict["Source_organism_strain_two"]
            )
            index = self.Expression_system_scientific_name_two.findText(
                self.deposit_dict["Expression_system_scientific_name_two"],
                QtCore.Qt.MatchFixedString,
            )
            self.Expression_system_scientific_name_two.setCurrentIndex(index)

            self.Expression_system_strain_two.setText(
                self.deposit_dict["Expression_system_strain_two"]
            )
            self.Expression_system_vector_type_two.setText(
                self.deposit_dict["Expression_system_vector_type_two"]
            )
            self.Expression_system_plasmid_name_two.setText(
                self.deposit_dict["Expression_system_plasmid_name_two"]
            )
            self.Manipulated_source_details_two.setText(
                self.deposit_dict["Manipulated_source_details_two"]
            )
            self.molecule_chain_two.setText(self.deposit_dict["molecule_chain_two"])
            self.molecule_one_letter_sequence_uniprot_id_two.setText(
                self.deposit_dict["molecule_two_letter_sequence_uniprot_id"]
            )
            self.molecule_one_letter_sequence_two.setText(
                self.deposit_dict["molecule_two_letter_sequence"]
            )

            self.structure_keywords.setText(self.deposit_dict["structure_keywords"])
            self.biological_assembly_chain_number.setText(
                self.deposit_dict["biological_assembly_chain_number"]
            )
            self.molecule_one_letter_sequence_uniprot_id.setText(
                self.deposit_dict["molecule_one_letter_sequence_uniprot_id"]
            )
            self.molecule_one_letter_sequence.setText(
                self.deposit_dict["molecule_one_letter_sequence"]
            )
            self.SG_project_name.setText(self.deposit_dict["SG_project_name"])
            self.full_name_of_SG_center.setText(
                self.deposit_dict["full_name_of_SG_center"]
            )

            index = self.crystallization_method.findText(
                self.deposit_dict["crystallization_method"], QtCore.Qt.MatchFixedString
            )
            self.crystallization_method.setCurrentIndex(index)

            self.crystallization_pH.setText(self.deposit_dict["crystallization_pH"])
            self.crystallization_temperature.setText(
                self.deposit_dict["crystallization_temperature"]
            )
            self.crystallization_details.setText(
                self.deposit_dict["crystallization_details"]
            )
            index = self.radiation_source.findText(
                self.deposit_dict["radiation_source"], QtCore.Qt.MatchFixedString
            )
            self.radiation_source.setCurrentIndex(index)

            index = self.radiation_source_type.findText(
                self.deposit_dict["radiation_source_type"], QtCore.Qt.MatchFixedString
            )
            self.radiation_source_type.setCurrentIndex(index)

            self.radiation_wavelengths.setText(
                self.deposit_dict["radiation_wavelengths"]
            )
            index = self.radiation_detector.findText(
                self.deposit_dict["radiation_detector"], QtCore.Qt.MatchFixedString
            )
            self.radiation_detector.setCurrentIndex(index)

            index = self.radiation_detector_type.findText(
                self.deposit_dict["radiation_detector_type"], QtCore.Qt.MatchFixedString
            )
            self.radiation_detector_type.setCurrentIndex(index)

            self.data_collection_date.setText(self.deposit_dict["data_collection_date"])
            self.data_collection_temperature.setText(
                self.deposit_dict["data_collection_temperature"]
            )
            self.data_collection_protocol.setText(
                self.deposit_dict["data_collection_protocol"]
            )

            self.pdbx_starting_model.setText(self.deposit_dict["pdbx_starting_model"])
            index = self.data_integration_software.findText(
                self.deposit_dict["data_integration_software"],
                QtCore.Qt.MatchFixedString,
            )
            self.data_integration_software.setCurrentIndex(index)
            index = self.phasing_software.findText(
                self.deposit_dict["phasing_software"], QtCore.Qt.MatchFixedString
            )
            self.phasing_software.setCurrentIndex(index)

            self.pdbx_funding_organization_one.setText(
                self.deposit_dict["pdbx_funding_organization_one"]
            )
            self.pdbx_grant_number_one.setText(
                self.deposit_dict["pdbx_grant_number_one"]
            )
            index = self.pdbx_grant_country_one.findText(
                self.deposit_dict["pdbx_grant_country_one"], QtCore.Qt.MatchFixedString
            )
            self.pdbx_grant_country_one.setCurrentIndex(index)

            self.pdbx_funding_organization_two.setText(
                self.deposit_dict["pdbx_funding_organization_two"]
            )
            self.pdbx_grant_number_two.setText(
                self.deposit_dict["pdbx_grant_number_two"]
            )
            index = self.pdbx_grant_country_two.findText(
                self.deposit_dict["pdbx_grant_country_two"], QtCore.Qt.MatchFixedString
            )
            self.pdbx_grant_country_two.setCurrentIndex(index)

            self.pdbx_funding_organization_three.setText(
                self.deposit_dict["pdbx_funding_organization_three"]
            )
            self.pdbx_grant_number_three.setText(
                self.deposit_dict["pdbx_grant_number_three"]
            )
            index = self.pdbx_grant_country_three.findText(
                self.deposit_dict["pdbx_grant_country_three"],
                QtCore.Qt.MatchFixedString,
            )
            self.pdbx_grant_country_three.setCurrentIndex(index)

        except ValueError as e:
            self.update_log.error("file is not a valid .deposit file: " + str(e))

    def update_deposit_dict(self):
        pdbx_funding_ordinal_one = ""
        pdbx_funding_organization_one = ""
        pdbx_grant_number_one = ""
        pdbx_grant_country_one = ""
        if str(self.pdbx_funding_organization_one.text()).replace(" ", "") != "":
            pdbx_funding_ordinal_one = "1"
            pdbx_funding_organization_one = str(
                self.pdbx_funding_organization_one.text()
            )
            pdbx_grant_number_one = str(self.pdbx_grant_number_one.text())
            pdbx_grant_country_one = str(self.pdbx_grant_country_one.currentText())

        pdbx_funding_ordinal_two = ""
        pdbx_funding_organization_two = ""
        pdbx_grant_number_two = ""
        pdbx_grant_country_two = ""
        if str(self.pdbx_funding_organization_two.text()).replace(" ", "") != "":
            pdbx_funding_ordinal_two = "2"
            pdbx_funding_organization_two = str(
                self.pdbx_funding_organization_two.text()
            )
            pdbx_grant_number_two = str(self.pdbx_grant_number_two.text())
            pdbx_grant_country_two = str(self.pdbx_grant_country_two.currentText())

        pdbx_funding_ordinal_three = ""
        pdbx_funding_organization_three = ""
        pdbx_grant_number_three = ""
        pdbx_grant_country_three = ""
        if str(self.pdbx_funding_organization_three.text()).replace(" ", "") != "":
            pdbx_funding_ordinal_three = "3"
            pdbx_funding_organization_three = str(
                self.pdbx_funding_organization_three.text()
            )
            pdbx_grant_number_three = str(self.pdbx_grant_number_three.text())
            pdbx_grant_country_three = str(self.pdbx_grant_country_three.currentText())

        self.deposit_dict = {
            "contact_author_PI_salutation": str(
                self.contact_author_PI_salutation.text()
            ),
            "contact_author_PI_first_name": str(
                self.contact_author_PI_first_name.text()
            ),
            "contact_author_PI_last_name": str(self.contact_author_PI_last_name.text()),
            "contact_author_PI_middle_name": str(
                self.contact_author_PI_middle_name.text()
            ),
            "contact_author_PI_role": str(self.contact_author_PI_role.currentText()),
            "contact_author_PI_organization_type": str(
                self.contact_author_PI_organization_type.currentText()
            ),
            "contact_author_PI_organization_name": str(
                self.contact_author_PI_organization_name.text()
            ),
            "contact_author_PI_email": str(self.contact_author_PI_email.text()),
            "contact_author_PI_address": str(self.contact_author_PI_address.text()),
            "contact_author_PI_city": str(self.contact_author_PI_city.text()),
            "contact_author_PI_State_or_Province": str(
                self.contact_author_PI_State_or_Province.text()
            ),
            "contact_author_PI_Zip_Code": str(self.contact_author_PI_Zip_Code.text()),
            "contact_author_PI_Country": str(self.contact_author_PI_Country.text()),
            "contact_author_PI_phone_number": str(
                self.contact_author_PI_phone_number.text()
            ),
            "contact_author_PI_ORCID": str(self.contact_author_PI_ORCID.text()),
            "contact_author_salutation": str(self.contact_author_salutation.text()),
            "contact_author_first_name": str(self.contact_author_first_name.text()),
            "contact_author_last_name": str(self.contact_author_last_name.text()),
            "contact_author_middle_name": str(self.contact_author_middle_name.text()),
            "contact_author_role": str(self.contact_author_role.currentText()),
            "contact_author_organization_type": str(
                self.contact_author_organization_type.currentText()
            ),
            "contact_author_organization_name": str(
                self.contact_author_organization_name.text()
            ),
            "contact_author_email": str(self.contact_author_email.text()),
            "contact_author_address": str(self.contact_author_address.text()),
            "contact_author_city": str(self.contact_author_city.text()),
            "contact_author_State_or_Province": str(
                self.contact_author_State_or_Province.text()
            ),
            "contact_author_Zip_Code": str(self.contact_author_Zip_Code.text()),
            "contact_author_Country": str(self.contact_author_Country.text()),
            "contact_author_phone_number": str(self.contact_author_phone_number.text()),
            "contact_author_ORCID": str(self.contact_author_ORCID.text()),
            "Release_status_for_coordinates": str(
                self.Release_status_for_coordinates.currentText()
            ),
            "Release_status_for_sequence": str(
                self.Release_status_for_sequence.currentText()
            ),
            "group_deposition_title": str(self.group_deposition_title.text()),
            "group_description": str(self.group_description.text()),
            "structure_title": str(self.structure_title.text()),
            "structure_title_apo": str(self.structure_title_apo.text()),
            "primary_citation_id": str(self.primary_citation_id.text()),
            "primary_citation_journal_abbrev": str(
                self.primary_citation_journal_abbrev.text()
            ),
            "primary_citation_title": str(self.primary_citation_title.text()),
            "primary_citation_year": str(self.primary_citation_year.text()),
            "primary_citation_journal_volume": str(
                self.primary_citation_journal_volume.text()
            ),
            "primary_citation_page_first": str(self.primary_citation_page_first.text()),
            "primary_citation_page_last": str(self.primary_citation_page_last.text()),
            # entity 1
            "molecule_name": str(self.molecule_name.text()),
            "Source_organism_scientific_name": str(
                self.Source_organism_scientific_name.currentText()
            ),
            "Source_organism_gene": str(self.Source_organism_gene.text()),
            "Source_organism_strain": str(self.Source_organism_strain.text()),
            "Expression_system_scientific_name": str(
                self.Expression_system_scientific_name.currentText()
            ),
            "Expression_system_strain": str(self.Expression_system_strain.text()),
            "Expression_system_plasmid_name": str(
                self.Expression_system_plasmid_name.text()
            ),
            "Expression_system_vector_type": str(
                self.Expression_system_vector_type.text()
            ),
            "Manipulated_source_details": str(self.Manipulated_source_details.text()),
            "fragment_name_one_specific_mutation": str(
                self.fragment_name_one_specific_mutation.text()
            ),
            "molecule_chain_one": str(self.molecule_chain_one.text()),
            # entity 2
            "molecule_name_two": str(self.molecule_name_two.text()),
            "Source_organism_scientific_name_two": str(
                self.Source_organism_scientific_name_two.currentText()
            ),
            "Source_organism_gene_two": str(self.Source_organism_gene_two.text()),
            "Source_organism_strain_two": str(self.Source_organism_strain_two.text()),
            "Expression_system_scientific_name_two": str(
                self.Expression_system_scientific_name_two.currentText()
            ),
            "Expression_system_strain_two": str(
                self.Expression_system_strain_two.text()
            ),
            "Expression_system_plasmid_name_two": str(
                self.Expression_system_plasmid_name_two.text()
            ),
            "Expression_system_vector_type_two": str(
                self.Expression_system_vector_type_two.text()
            ),
            "Manipulated_source_details_two": str(
                self.Manipulated_source_details_two.text()
            ),
            "fragment_name_two_specific_mutation": str(
                self.fragment_name_two_specific_mutation.text()
            ),
            "molecule_chain_two": str(self.molecule_chain_two.text()),
            "structure_keywords": str(self.structure_keywords.text()),
            "biological_assembly_chain_number": str(
                self.biological_assembly_chain_number.text()
            ),
            "molecule_one_letter_sequence_uniprot_id": str(
                self.molecule_one_letter_sequence_uniprot_id.text()
            ),
            "molecule_two_letter_sequence_uniprot_id": str(
                self.molecule_one_letter_sequence_uniprot_id_two.text()
            ),
            "SG_project_name": str(self.SG_project_name.text()),
            "full_name_of_SG_center": str(self.full_name_of_SG_center.text()),
            "molecule_one_letter_sequence": str(
                self.molecule_one_letter_sequence.toPlainText()
            )
            .replace(" ", "")
            .replace("\n", "")
            .replace("\r", ""),
            "molecule_two_letter_sequence": str(
                self.molecule_one_letter_sequence_two.toPlainText()
            )
            .replace(" ", "")
            .replace("\n", "")
            .replace("\r", ""),
            "crystallization_method": str(self.crystallization_method.currentText()),
            "crystallization_pH": str(self.crystallization_pH.text()),
            "crystallization_temperature": str(self.crystallization_temperature.text()),
            "crystallization_details": str(self.crystallization_details.text()),
            "radiation_source": str(self.radiation_source.currentText()),
            "radiation_source_type": str(self.radiation_source_type.currentText()),
            "radiation_wavelengths": str(self.radiation_wavelengths.text()),
            "radiation_detector": str(self.radiation_detector.currentText()),
            "radiation_detector_type": str(self.radiation_detector_type.currentText()),
            "data_collection_date": str(self.data_collection_date.text()),
            "data_collection_temperature": str(self.data_collection_temperature.text()),
            "data_collection_protocol": str(self.data_collection_protocol.text()),
            "pdbx_starting_model": str(self.pdbx_starting_model.text()),
            "data_integration_software": str(
                self.data_integration_software.currentText()
            ),
            "phasing_software": str(self.phasing_software.currentText()),
            "pdbx_funding_ordinal_one": pdbx_funding_ordinal_one,
            "pdbx_funding_organization_one": pdbx_funding_organization_one,
            "pdbx_grant_number_one": pdbx_grant_number_one,
            "pdbx_grant_country_one": pdbx_grant_country_one,
            "pdbx_funding_ordinal_two": pdbx_funding_ordinal_two,
            "pdbx_funding_organization_two": pdbx_funding_organization_two,
            "pdbx_grant_number_two": pdbx_grant_number_two,
            "pdbx_grant_country_two": pdbx_grant_country_two,
            "pdbx_funding_ordinal_three": pdbx_funding_ordinal_three,
            "pdbx_funding_organization_three": pdbx_funding_organization_three,
            "pdbx_grant_number_three": pdbx_grant_number_three,
            "pdbx_grant_country_three": pdbx_grant_country_three,
        }

        structure_author_name = ""
        for widget in self.structure_author_name_List:
            structure_author_name += str(widget.text()) + ";"
        self.deposit_dict["structure_author_name"] = structure_author_name[:-1]

        primary_citation_author_name = ""
        for widget in self.primary_citation_author_name_List:
            primary_citation_author_name += str(widget.text()) + ";"
        self.deposit_dict[
            "primary_citation_author_name"
        ] = primary_citation_author_name[:-1]

    def get_deposit_dict_template(self):
        deposit_dict_template = {
            "contact_author_PI_salutation": None,
            "contact_author_PI_first_name": None,
            "contact_author_PI_last_name": None,
            "contact_author_PI_middle_name": None,
            "contact_author_PI_role": None,
            "contact_author_PI_organization_type": None,
            "contact_author_PI_organization_name": None,
            "contact_author_PI_email": None,
            "contact_author_PI_address": None,
            "contact_author_PI_city": None,
            "contact_author_PI_State_or_Province": None,
            "contact_author_PI_Zip_Code": None,
            "contact_author_PI_Country": None,
            "contact_author_PI_phone_number": None,
            "contact_author_PI_ORCID": None,
            "contact_author_salutation": None,
            "contact_author_first_name": None,
            "contact_author_last_name": None,
            "contact_author_middle_name": None,
            "contact_author_role": None,
            "contact_author_organization_type": None,
            "contact_author_organization_name": None,
            "contact_author_email": None,
            "contact_author_address": None,
            "contact_author_city": None,
            "contact_author_State_or_Province": None,
            "contact_author_Zip_Code": None,
            "contact_author_Country": None,
            "contact_author_phone_number": None,
            "contact_author_ORCID": None,
            "Release_status_for_coordinates": None,
            "Release_status_for_sequence": None,
            "group_deposition_title": None,
            "group_description": None,
            "structure_title": None,
            "structure_title_apo": None,
            "primary_citation_id": None,
            "primary_citation_journal_abbrev": None,
            "primary_citation_title": None,
            "primary_citation_year": None,
            "primary_citation_journal_volume": None,
            "primary_citation_page_first": None,
            "primary_citation_page_last": None,
            # entity 1
            "molecule_name": None,
            "Source_organism_scientific_name": None,
            "Source_organism_gene": None,
            "Source_organism_strain": None,
            "Expression_system_scientific_name": None,
            "Expression_system_strain": None,
            "Expression_system_plasmid_name": None,
            "Expression_system_vector_type": None,
            "Manipulated_source_details": None,
            "fragment_name_one_specific_mutation": None,
            "molecule_chain_one": None,
            # entity 2
            "molecule_name_two": None,
            "Source_organism_scientific_name_two": None,
            "Source_organism_gene_two": None,
            "Source_organism_strain_two": None,
            "Expression_system_scientific_name_two": None,
            "Expression_system_strain_two": None,
            "Expression_system_plasmid_name_two": None,
            "Expression_system_vector_type_two": None,
            "Manipulated_source_details_two": None,
            "fragment_name_two_specific_mutation": None,
            "molecule_chain_two": None,
            "structure_keywords": None,
            "biological_assembly_chain_number": None,
            "molecule_one_letter_sequence_uniprot_id": None,
            "molecule_two_letter_sequence_uniprot_id": None,
            "SG_project_name": None,
            "full_name_of_SG_center": None,
            "molecule_one_letter_sequence": None,
            "molecule_two_letter_sequence": None,
            "crystallization_method": None,
            "crystallization_pH": None,
            "crystallization_temperature": None,
            "crystallization_details": None,
            "radiation_source": None,
            "radiation_source_type": None,
            "radiation_wavelengths": None,
            "radiation_detector": None,
            "radiation_detector_type": None,
            "data_collection_date": None,
            "data_collection_temperature": None,
            "data_collection_protocol": None,
            "pdbx_starting_model": None,
            "data_integration_software": None,
            "phasing_software": None,
            "structure_author_name": None,
            "primary_citation_author_name": None,
            "pdbx_funding_organization_one": "",
            "pdbx_grant_number_one": "",
            "pdbx_grant_country_one": "",
            "pdbx_funding_organization_two": "",
            "pdbx_grant_number_two": "",
            "pdbx_grant_country_two": "",
            "pdbx_funding_organization_three": "",
            "pdbx_grant_number_three": "",
            "pdbx_grant_country_three": "",
        }

        return deposit_dict_template

    def set_primary_citation_as_structure_authors(self, state):
        if state == QtCore.Qt.Checked:
            for n, entry in enumerate(self.structure_author_name_List):
                self.primary_citation_author_name_List[n].setText(str(entry.text()))
        else:
            for n, entry in enumerate(self.primary_citation_author_name_List):
                entry.setText("")

    def set_xce_logfile(self):
        file_name = str(
            QtGui.QFileDialog.getSaveFileName(
                self.window, "Save file", self.current_directory
            )
        )
        self.xce_logfile = str(file_name)
        self.xce_logfile_label.setText(str(self.xce_logfile))
        if (
            self.xce_logfile == ""
            or self.xce_logfile[self.xce_logfile.rfind("/") + 1 :] == ""
        ):
            print("==> XCE: invalid file format")
        else:
            XChemLog.startLog(self.xce_logfile).create_logfile(self.xce_version)
            self.update_log = XChemLog.updateLog(self.xce_logfile)

    def set_second_cif_file(self):
        filepath_temp = QtGui.QFileDialog.getOpenFileNameAndFilter(
            self.window, "Select CIF File", self.initial_model_directory, "*.cif"
        )
        filepath = str(tuple(filepath_temp)[0])
        self.second_cif_file = str(filepath)
        self.second_cif_file_label.setText(str(self.second_cif_file))
        self.update_log.insert(
            "user selected %s as CIF file for merging into ligand CIF files"
            % self.second_cif_file
        )

    def select_datasource_columns_to_display(self):
        columns_to_show = QtGui.QMessageBox()
        columns_to_showLayout = columns_to_show.layout()
        columns_in_data_source = self.db.return_column_list()
        try:
            columns_in_data_source = self.db.return_column_list()
        except AttributeError:
            print("==> XCE: please select a datasource file")
            self.status_bar.showMessage("please select a datasource file")
            return

        column_dict = {}
        vbox = QtGui.QVBoxLayout()
        number_of_entries = len(columns_in_data_source)
        columns_shown_in_dialog_column = 15
        grid = QtGui.QGridLayout()
        x = 0
        y = 0
        columns_to_ignore = self.db.columns_not_to_display()
        for entries_added in range(number_of_entries):
            if not columns_in_data_source[entries_added][1] in columns_to_ignore:
                data_source_column = QtGui.QCheckBox(
                    columns_in_data_source[entries_added][1]
                )
                column_dict[entries_added] = data_source_column
                if (
                    columns_in_data_source[entries_added][1]
                    in self.overview_datasource_table_columns
                ):
                    data_source_column.setChecked(True)
                grid.addWidget(data_source_column, y, x)
                y += 1
            if y == columns_shown_in_dialog_column:
                y = 0
                x += 1
        vbox.addLayout(grid)
        columns_to_showLayout.addLayout(vbox, 0, 0)

        columns_to_show.addButton(QtGui.QPushButton("OK"), QtGui.QMessageBox.YesRole)
        columns_to_show.addButton(
            QtGui.QPushButton("Cancel"), QtGui.QMessageBox.RejectRole
        )
        reply = columns_to_show.exec_()
        if reply == 0:
            columns_to_show_list = ["Sample ID"]
            for key in column_dict:
                if column_dict[key].isChecked():
                    columns_to_show_list.append(columns_in_data_source[key][1])
            self.overview_datasource_table_columns = columns_to_show_list
            self.populate_and_update_datasource_table()

    def update_header_and_data_from_datasource(self):
        self.update_log.insert(
            "getting information for all samples from data source..."
        )
        self.db = XChemDB.data_source(
            os.path.join(self.database_directory, self.data_source_file)
        )
        self.update_log.insert("creating missing columns in data source")
        self.db.create_missing_columns()
        self.update_log.insert("load header and data from data source")
        self.header, self.data = self.db.load_samples_from_data_source()
        self.update_log.insert("get all samples in data source")

        self.xtal_db_dict = {}
        sampleID_column = 0
        for n, entry in enumerate(self.header):
            if entry == "CrystalName":
                sampleID_column = n
                break
        for line in self.data:
            if str(line[sampleID_column]) != "":
                db_dict = {}
                for n, entry in enumerate(line):
                    if n != sampleID_column:
                        db_dict[str(self.header[n])] = str(entry)
                self.xtal_db_dict[str(line[sampleID_column])] = db_dict

        print(("==> XCE: found " + str(len(self.xtal_db_dict)) + " samples"))

    def datasource_menu_save_samples(self):
        print("hallo")

    def datasource_menu_export_csv_file(self):
        file_name = str(
            QtGui.QFileDialog.getSaveFileName(
                self.window, "Save file", self.database_directory
            )
        )
        if file_name.rfind(".") != -1:
            file_name = file_name[: file_name.rfind(".")] + ".csv"
        else:
            file_name = file_name + ".csv"
        self.db.export_to_csv_file(file_name)

    def datasource_menu_import_csv_file(self):
        if self.data_source_set:
            file_name = QtGui.QFileDialog.getOpenFileName(
                self.window, "Open file", self.database_directory
            )
            self.db.import_csv_file(file_name)
        else:
            self.update_status_bar("Please load a data source file first")

    def datasource_menu_update_datasource(self):
        self.work_thread = XChemThread.synchronise_db_and_filesystem(
            self.initial_model_directory,
            os.path.join(self.database_directory, self.data_source_file),
            self.panddas_directory,
            self.xce_logfile,
            "project_directory",
        )
        self.connect(
            self.work_thread,
            QtCore.SIGNAL("update_progress_bar"),
            self.update_progress_bar,
        )
        self.connect(
            self.work_thread,
            QtCore.SIGNAL("update_status_bar(QString)"),
            self.update_status_bar,
        )
        self.connect(
            self.work_thread, QtCore.SIGNAL("finished()"), self.thread_finished
        )
        self.connect(
            self.work_thread,
            QtCore.SIGNAL("datasource_menu_reload_samples"),
            self.datasource_menu_reload_samples,
        )
        self.work_thread.start()

    def export_data_for_WONKA(self):
        self.update_log.insert("exporting CSV file for input into WONKA")
        self.db.export_csv_for_WONKA()

    def populate_reference_combobox(self, combobox):
        combobox.clear()
        for reference_file in self.reference_file_list:
            combobox.addItem(reference_file[0])

    def populate_refinement_outcome_combobox(self, combobox):
        combobox.clear()
        for stage in self.refinement_stage:
            combobox.addItem(stage)

    def populate_target_selection_combobox(self, combobox):
        combobox.clear()
        for target in self.target_list:
            combobox.addItem(target)

    def combo_selected(self, text):
        self.map_url = str(
            self.panddas_directory
            + "/analyses/html_summaries/pandda_map_"
            + text
            + ".html"
        )
        self.pandda_maps_html.load(QtCore.QUrl(self.map_url))
        self.pandda_maps_html.show()

    def add_map_html(self):
        self.map_list = glob.glob(
            str(self.panddas_directory + "/analyses/html_summaries/pandda_map_*.html")
        )
        self.list_options = []
        for i in range(0, len(self.map_list)):
            string = self.map_list[i]
            string = string.replace("/analyses/html_summaries/pandda_map_", "")
            string = string.replace(".html", "")
            string = string.replace(self.panddas_directory, "")
            self.list_options.append(string)
        self.pandda_map_list.clear()
        for i in range(0, len(self.list_options)):
            self.pandda_map_list.addItem(self.list_options[i])
        self.connect(
            self.pandda_map_list,
            QtCore.SIGNAL("activated(QString)"),
            self.combo_selected,
        )

    def open_config_file(self):
        file_name_temp = QtGui.QFileDialog.getOpenFileNameAndFilter(
            self.window, "Open file", self.current_directory, "*.conf"
        )
        file_name = tuple(file_name_temp)[0]

        try:
            pickled_settings = pickle.load(open(file_name, "rb"))
        except Exception:
            print("==> XCE: failed to open config file...")

        key_list = {
            "initial_model_directory": "initial_model_directory",
            "panddas_directory": "panddas_directory",
            "html_export_directory": "html_export_directory",
            "group_deposit_directory": "group_deposit_directory",
            "database_directory": "database_directory",
            "datasets_summary_file": "datasets_summary",
            "ccp4_scratch_directory": "ccp4_scratch",
            "allowed_unitcell_difference_percent": "unitcell_difference",
            "acceptable_low_resolution_limit_for_data": "too_low_resolution_data",
        }

        for current_key in key_list:
            try:
                command = str(
                    "self."
                    + current_key
                    + " = pickled_settings['"
                    + key_list[current_key]
                    + "']"
                )
                exec(command)
                command = str(
                    'self.settings["'
                    + key_list[current_key]
                    + '"]= self.'
                    + current_key
                )
                exec(command)
                print(("==> XCE: found " + key_list[current_key]))
            except Exception:
                print(
                    (
                        "==> XCE: WARNING: Failed to find settings for: "
                        + key_list[current_key]
                        + " Error type: "
                        + str(sys.exc_info()[0])
                    )
                )
                exec(str(current_key + " = ''"))
                continue

        try:
            pickled_settings = pickle.load(open(file_name, "rb"))
            if pickled_settings["beamline_directory"] != self.beamline_directory:
                self.beamline_directory = pickled_settings["beamline_directory"]
                self.target_list, self.visit_list = XChemMain.get_target_and_visit_list(
                    self.beamline_directory, self.read_agamemnon.isChecked()
                )
                self.settings["beamline_directory"] = self.beamline_directory
                self.populate_target_selection_combobox(self.target_selection_combobox)

            layout_functions.pandda_html(self)
            self.show_pandda_html_summary()

            self.html_export_directory_label.setText(self.html_export_directory)

            self.group_deposition_directory_label.setText(self.group_deposit_directory)

            self.datasets_summary_file_label.setText(self.datasets_summary_file)

            self.data_source_file = pickled_settings["data_source"]
            if self.data_source_file != "":
                self.settings["data_source"] = os.path.join(
                    self.database_directory, self.data_source_file
                )
                # this is probably not necessary
                if os.path.isfile(self.settings["data_source"]):
                    write_enabled = self.check_write_permissions_of_data_source()
                    if not write_enabled:
                        self.data_source_file_label.setText("")
                        self.data_source_set = False
                    else:
                        self.data_source_file_label.setText(
                            os.path.join(self.database_directory, self.data_source_file)
                        )
                        self.data_source_set = True
                        self.db = XChemDB.data_source(
                            os.path.join(self.database_directory, self.data_source_file)
                        )
                        self.datasource_menu_reload_samples()

            reference_directory_temp = pickled_settings["reference_directory"]
            if reference_directory_temp != self.reference_directory:
                self.reference_directory = reference_directory_temp
                self.settings["reference_directory"] = self.reference_directory
                self.update_reference_files(" ")
                for xtal in self.initial_model_dimple_dict:
                    reference_file_selection_combobox = self.initial_model_dimple_dict[
                        xtal
                    ][1]
                    self.populate_reference_combobox(reference_file_selection_combobox)

            self.initial_model_directory_label.setText(self.initial_model_directory)
            self.panddas_directory_label.setText(self.panddas_directory)
            self.pandda_output_data_dir_entry.setText(self.panddas_directory)
            self.reference_directory_label.setText(self.reference_directory)
            self.beamline_directory_label.setText(self.beamline_directory)
            self.ccp4_scratch_directory_label.setText(self.ccp4_scratch_directory)
            self.reference_file_list = self.get_reference_file_list(" ")
            self.pandda_input_data_dir_entry.setText(
                os.path.join(self.initial_model_directory, "*")
            )

            self.update_all_tables()

        except KeyError:
            self.update_status_bar("Sorry, this is not a XChemExplorer config file!")
            self.update_log.insert("Sorry, this is not a XChemExplorer config file!")

        except Exception as exception:
            print(("Unexpected error:", sys.exc_info()[0]))
            raise exception

    def save_config_file(self):
        file_name = str(
            QtGui.QFileDialog.getSaveFileName(
                self.window, "Save file", self.current_directory
            )
        )
        # make sure that the file always has .conf extension
        if str(file_name).rfind(".") != -1:
            file_name = file_name[: file_name.rfind(".")] + ".conf"
        else:
            file_name = file_name + ".conf"
        pickle.dump(self.settings, open(file_name, "wb"))

    def update_reference_files(self, reference_root):
        self.reference_file_list = self.get_reference_file_list(reference_root)
        self.populate_reference_combobox(self.reference_file_selection_combobox)
        self.populate_reference_combobox(self.pandda_reference_file_selection_combobox)

    def run_dimple_on_selected_autoprocessing_file(self, instruction):
        job_list = []
        for xtal in sorted(self.initial_model_dimple_dict):
            if self.initial_model_dimple_dict[xtal][0].isChecked():
                db_dict = self.xtal_db_dict[xtal]

                # the if statement below is so convoluted, so that it is compatible
                # with older data source files

                if (
                    os.path.isfile(
                        os.path.join(
                            db_dict["ProjectDirectory"],
                            xtal,
                            db_dict["DataProcessingPathToMTZfile"],
                            db_dict["DataProcessingMTZfileName"],
                        )
                    )
                    or os.path.isfile(
                        os.path.join(
                            db_dict["ProjectDirectory"],
                            xtal,
                            db_dict["DataProcessingPathToMTZfile"],
                        )
                    )
                    or os.path.isfile(
                        os.path.join(
                            db_dict["DataProcessingPathToMTZfile"],
                            db_dict["DataProcessingMTZfileName"],
                        )
                    )
                    or os.path.isfile(
                        os.path.join(db_dict["DataProcessingPathToMTZfile"])
                    )
                ):
                    if os.path.isfile(
                        os.path.join(
                            db_dict["DataProcessingPathToMTZfile"],
                            db_dict["DataProcessingMTZfileName"],
                        )
                    ):
                        mtzin = os.path.join(
                            db_dict["DataProcessingPathToMTZfile"],
                            db_dict["DataProcessingMTZfileName"],
                        )
                    elif os.path.isfile(
                        os.path.join(db_dict["DataProcessingPathToMTZfile"])
                    ):
                        mtzin = os.path.join(db_dict["DataProcessingPathToMTZfile"])
                    elif os.path.isfile(
                        os.path.join(
                            db_dict["ProjectDirectory"],
                            xtal,
                            db_dict["DataProcessingPathToMTZfile"],
                            db_dict["DataProcessingMTZfileName"],
                        )
                    ):
                        mtzin = os.path.join(
                            db_dict["ProjectDirectory"],
                            xtal,
                            db_dict["DataProcessingPathToMTZfile"],
                            db_dict["DataProcessingMTZfileName"],
                        )
                    elif os.path.isfile(
                        os.path.join(
                            db_dict["ProjectDirectory"],
                            xtal,
                            db_dict["DataProcessingPathToMTZfile"],
                        )
                    ):
                        mtzin = os.path.join(
                            db_dict["ProjectDirectory"],
                            xtal,
                            db_dict["DataProcessingPathToMTZfile"],
                        )

                    reference_file = str(
                        self.initial_model_dimple_dict[xtal][1].currentText()
                    )

                    reference_file_pdb = os.path.join(
                        self.reference_directory, reference_file + ".pdb"
                    )

                    if not os.path.isfile(reference_file_pdb):
                        continue

                    if os.path.isfile(
                        os.path.join(self.reference_directory, reference_file + ".mtz")
                    ):
                        reference_file_mtz = " -R " + os.path.join(
                            self.reference_directory, reference_file + ".mtz"
                        )
                    else:
                        reference_file_mtz = ""

                    if os.path.isfile(
                        os.path.join(self.reference_directory, reference_file + ".cif")
                    ):
                        reference_file_cif = " --libin " + os.path.join(
                            self.reference_directory, reference_file + ".cif"
                        )
                    else:
                        reference_file_cif = ""

                    job_list.append(
                        [
                            xtal,
                            "dimple_rerun_on_selected_file",
                            mtzin,
                            reference_file_pdb,
                            reference_file_mtz,
                            reference_file_cif,
                        ]
                    )
                else:
                    print(
                        (
                            "WARNING: "
                            + xtal
                            + " has not been submitted to dimple"
                            + " because no files were found: "
                        )
                    )
                    if not os.path.isfile(
                        os.path.join(
                            db_dict["ProjectDirectory"],
                            xtal,
                            db_dict["DataProcessingPathToMTZfile"],
                            db_dict["DataProcessingMTZfileName"],
                        )
                    ):
                        print(
                            (
                                "    "
                                + str(
                                    os.path.join(
                                        db_dict["ProjectDirectory"],
                                        xtal,
                                        db_dict["DataProcessingPathToMTZfile"],
                                        db_dict["DataProcessingMTZfileName"],
                                    )
                                )
                                + " is missing"
                            )
                        )
                    if not os.path.isfile(
                        os.path.join(
                            db_dict["ProjectDirectory"],
                            xtal,
                            db_dict["DataProcessingPathToMTZfile"],
                        )
                    ):
                        print(
                            (
                                "    "
                                + str(
                                    os.path.join(
                                        db_dict["ProjectDirectory"],
                                        xtal,
                                        db_dict["DataProcessingPathToMTZfile"],
                                    )
                                )
                                + " is missing"
                            )
                        )
                    if not os.path.isfile(
                        os.path.join(db_dict["DataProcessingPathToMTZfile"])
                    ):
                        print(
                            (
                                "    "
                                + str(
                                    os.path.join(db_dict["DataProcessingPathToMTZfile"])
                                    + " is missing"
                                )
                            )
                        )

        if job_list:
            self.update_log.insert(
                "trying to run DIMPLE on SELECTED auto-processing files"
            )
            self.check_before_running_dimple(job_list, instruction)

    def remove_selected_dimple_files(self, instruction):
        if "dimple" in instruction.lower():
            pipeline = "dimple"
        elif "pipedream" in instruction.lower():
            pipeline = "pipedream"
        elif "phenix" in instruction.lower():
            pipeline = "phenix.ligand_pipeline"

        job_list = []
        for xtal in sorted(self.initial_model_dimple_dict):
            if self.initial_model_dimple_dict[xtal][0].isChecked():
                job_list.append(xtal)

        if job_list:
            msgBox = QtGui.QMessageBox()
            msgBox.setText(
                "Do you really want to delete {0!s} {1!s} files?".format(
                    len(job_list), self.preferences["initial_refinement_pipeline"]
                )
            )
            msgBox.addButton(QtGui.QPushButton("Go"), QtGui.QMessageBox.YesRole)
            msgBox.addButton(QtGui.QPushButton("Cancel"), QtGui.QMessageBox.RejectRole)
            reply = msgBox.exec_()

            if reply == 0:
                self.status_bar.showMessage(
                    "preparing to remove {0!s} files".format(pipeline)
                )
                self.update_log.insert(
                    "preparing to remove {0!s} files".format(pipeline)
                )
                self.work_thread = XChemThread.remove_selected_dimple_files(
                    job_list,
                    self.initial_model_directory,
                    self.xce_logfile,
                    self.database_directory,
                    self.data_source_file,
                    pipeline,
                )
                self.explorer_active = 1
                self.connect(
                    self.work_thread, QtCore.SIGNAL("finished()"), self.thread_finished
                )
                self.connect(
                    self.work_thread,
                    QtCore.SIGNAL("update_progress_bar"),
                    self.update_progress_bar,
                )
                self.connect(
                    self.work_thread,
                    QtCore.SIGNAL("update_status_bar(QString)"),
                    self.update_status_bar,
                )
                self.connect(
                    self.work_thread, QtCore.SIGNAL("finished()"), self.thread_finished
                )
                self.connect(
                    self.work_thread,
                    QtCore.SIGNAL("datasource_menu_reload_samples"),
                    self.datasource_menu_reload_samples,
                )
                self.work_thread.start()

    def set_results_from_selected_pipeline(self, instruction):
        if "dimple" in instruction.lower():
            pipeline = "dimple"
        elif "pipedream" in instruction.lower():
            pipeline = "pipedream"
        elif "phenix" in instruction.lower():
            pipeline = "phenix.ligand_pipeline"

        self.update_log.warning("selecting initial refinement results from " + pipeline)

        job_list = []
        for xtal in sorted(self.initial_model_dimple_dict):
            if self.initial_model_dimple_dict[xtal][0].isChecked():
                job_list.append(xtal)

        self.work_thread = XChemThread.set_results_from_selected_pipeline(
            job_list,
            self.initial_model_directory,
            self.xce_logfile,
            self.database_directory,
            self.data_source_file,
            pipeline,
        )
        self.explorer_active = 1
        self.connect(
            self.work_thread, QtCore.SIGNAL("finished()"), self.thread_finished
        )
        self.connect(
            self.work_thread,
            QtCore.SIGNAL("update_progress_bar"),
            self.update_progress_bar,
        )
        self.connect(
            self.work_thread,
            QtCore.SIGNAL("update_status_bar(QString)"),
            self.update_status_bar,
        )
        self.connect(
            self.work_thread, QtCore.SIGNAL("finished()"), self.thread_finished
        )
        self.connect(
            self.work_thread,
            QtCore.SIGNAL("datasource_menu_reload_samples"),
            self.datasource_menu_reload_samples,
        )
        self.work_thread.start()

    def check_before_running_dimple(self, job_list, instruction):
        msgBox = QtGui.QMessageBox()
        msgBox.setText(
            "Do you really want to run {0!s} {1!s} jobs?\nNote: we will not run more"
            " than {2!s} at once on the cluster!".format(
                len(job_list),
                self.preferences["initial_refinement_pipeline"],
                self.preferences["max_queue_jobs"],
            )
        )
        msgBox.addButton(QtGui.QPushButton("Go"), QtGui.QMessageBox.YesRole)
        msgBox.addButton(QtGui.QPushButton("Cancel"), QtGui.QMessageBox.RejectRole)
        reply = msgBox.exec_()

        if reply == 0:
            if "dimple" in instruction.lower():
                pipeline = "dimple"
            elif "pipedream" in instruction.lower():
                pipeline = "pipedream"
            elif "phenix" in instruction.lower():
                pipeline = "phenix.ligand_pipeline"

            self.status_bar.showMessage(
                "preparing {0!s} DIMPLE jobs".format(len(job_list))
            )
            self.update_log.insert(
                "preparing to run {0!s} DIMPLE jobs".format(len(job_list))
            )
            self.update_log.insert("preparing input files for DIMPLE...")
            self.work_thread = XChemThread.run_dimple_on_all_autoprocessing_files_new(
                job_list,
                self.initial_model_directory,
                self.external_software,
                self.ccp4_scratch_directory,
                self.database_directory,
                self.data_source_file,
                self.max_queue_jobs,
                self.xce_logfile,
                self.preferences["dimple_twin_mode"],
                pipeline,
                get_token(),
            )
            self.explorer_active = 1
            self.connect(
                self.work_thread, QtCore.SIGNAL("finished()"), self.thread_finished
            )
            self.connect(
                self.work_thread,
                QtCore.SIGNAL("update_progress_bar"),
                self.update_progress_bar,
            )
            self.connect(
                self.work_thread,
                QtCore.SIGNAL("update_status_bar(QString)"),
                self.update_status_bar,
            )
            self.connect(
                self.work_thread, QtCore.SIGNAL("finished()"), self.thread_finished
            )
            self.connect(
                self.work_thread,
                QtCore.SIGNAL("datasource_menu_reload_samples"),
                self.datasource_menu_reload_samples,
            )
            self.work_thread.start()

    def update_all_tables(self):
        self.update_log.insert("checking for new reference files")
        self.update_status_bar("checking for new reference files")
        self.reference_file_list = self.get_reference_file_list(" ")
        self.update_log.insert("updating Overview table")
        self.update_status_bar("updating Overview table")
        self.populate_and_update_datasource_table()
        self.update_log.insert("updating Maps table")
        self.update_status_bar("updating Maps table")
        self.create_maps_table()
        self.update_log.insert("updating PANDDA table")
        self.update_status_bar("updating PANDDA table")
        self.populate_pandda_analyse_input_table()
        self.update_log.insert("updating REFINEMENT table")
        self.update_status_bar("updating REFINEMENT table")
        self.populate_and_update_refinement_table()
        self.update_status_bar("idle")
        self.update_summary_plot()

    def change_allowed_unitcell_difference_percent(self, text):
        try:
            self.allowed_unitcell_difference_percent = int(text)
            self.settings[
                "unitcell_difference"
            ] = self.allowed_unitcell_difference_percent
            self.update_log.insert(
                "changing max allowed unit cell difference between reference and xtal"
                " to {0!s} percent".format(self.allowed_unitcell_difference_percent)
            )
        except ValueError:
            if str(text).find(".") != -1:
                self.allowed_unitcell_difference_percent = int(
                    str(text)[: str(text).find(".")]
                )
                self.settings[
                    "unitcell_difference"
                ] = self.allowed_unitcell_difference_percent
                self.update_log.insert(
                    "changing max allowed unit cell difference between reference and"
                    " xtal to {0!s} percent".format(
                        self.allowed_unitcell_difference_percent
                    )
                )
            else:
                pass

    def change_max_queue_jobs(self, text):
        try:
            self.max_queue_jobs = int(text)
            self.settings["max_queue_jobs"] = self.max_queue_jobs
            self.update_log.insert(
                "changing max number of jobs running simultaneously on DLS cluster to"
                " {0!s}".format(self.max_queue_jobs)
            )
        except ValueError:
            if str(text).find(".") != -1:
                self.max_queue_jobs = int(str(text)[: str(text).find(".")])
                self.settings["max_queue_jobs"] = self.max_queue_jobs
                self.update_log.insert(
                    "changing max number of jobs running simultaneously on DLS cluster"
                    " to {0!s}".format(self.max_queue_jobs)
                )
            else:
                pass

    def change_acceptable_low_resolution_limit(self, text):
        try:
            self.acceptable_low_resolution_limit_for_data = float(text)
            self.settings[
                "too_low_resolution_data"
            ] = self.acceptable_low_resolution_limit_for_data
        except ValueError:
            pass

    def change_filename_root(self, text):
        self.filename_root = str(text)
        self.settings["filename_root"] = self.filename_root

    def button_clicked(self):
        if not self.data_source_set:
            print("sender text bit")
            if self.sender().text() == "Create New Data\nSource (SQLite)":
                file_name = str(
                    QtGui.QFileDialog.getSaveFileName(
                        self.window, "Save file", self.database_directory
                    )
                )
                # make sure that the file always has .sqlite extension
                if file_name.rfind(".") != -1:
                    file_name = file_name[: file_name.rfind(".")] + ".sqlite"
                else:
                    file_name = file_name + ".sqlite"
                self.db = XChemDB.data_source(file_name)
                print("==> XCE: creating new data source")
                self.db.create_empty_data_source_file()
                self.db.create_missing_columns()
                if self.data_source_file == "":
                    self.database_directory = file_name[: file_name.rfind("/")]
                    self.data_source_file = file_name[file_name.rfind("/") + 1 :]
                    self.data_source_file_label.setText(
                        os.path.join(self.database_directory, self.data_source_file)
                    )
                    self.settings["database_directory"] = self.database_directory
                    self.settings["data_source"] = self.data_source_file
                    self.data_source_set = True
            else:
                self.no_data_source_selected()
                print("No datasource selected")
                pass

        # first find out which of the 'Run' or 'Status' buttons is sending
        for item in self.workflow_widget_dict:
            for widget in self.workflow_widget_dict[item]:
                if widget == self.sender():
                    # get index of item in self.workflow; Note this index should be the
                    # same as the index of the self.main_tab_widget which belongs to
                    # this task
                    task_index = self.workflow.index(item)
                    instruction = str(self.workflow_widget_dict[item][0].currentText())
                    print(instruction)
                    action = str(self.sender().text())
                    if self.main_tab_widget.currentIndex() == task_index:
                        if self.explorer_active == 0 and self.data_source_set is True:
                            if action == "Run":
                                self.prepare_and_run_task(instruction)
                            elif action == "Status":
                                self.get_status_of_workflow_milestone(instruction)
                                if os.path.exists(
                                    str(self.panddas_directory + "/pandda.done")
                                ):
                                    self.pandda_status = "Finished!"
                                    self.pandda_status_label.setStyleSheet(
                                        "color: green"
                                    )
                                if os.path.exists(
                                    str(self.panddas_directory + "/pandda.running")
                                ):
                                    self.pandda_status = "Running..."
                                    self.pandda_status_label.setStyleSheet(
                                        "color: orange"
                                    )
                                if os.path.exists(
                                    str(self.panddas_directory + "/pandda.errored")
                                ):
                                    self.pandda_status = (
                                        "Error encountered..."
                                        " please check the log files for pandda!"
                                    )
                                    self.pandda_status_label.setStyleSheet("color: red")
                                self.pandda_status_label.setText(
                                    str("STATUS: " + self.pandda_status)
                                )
                    else:
                        self.need_to_switch_main_tab(task_index)

    def get_status_of_workflow_milestone(self, instruction):
        # first update all tables
        self.datasource_menu_reload_samples()

        cluster_dict = XChemMain.get_jobs_running_on_cluster(
            self.external_software, self.xce_logfile, get_token()
        )

        self.update_log.insert("getting status updates...")

        self.status_bar.showMessage(
            "please check terminal window for further information"
        )

        self.update_log.insert(
            "{0!s} samples are currently in database".format(
                str(len(self.xtal_db_dict))
            )
        )

        if "DIMPLE" in instruction:
            XChemMain.print_cluster_status_message(
                "dimple", cluster_dict, self.xce_logfile
            )

        elif "Create CIF/PDB/PNG file" in instruction:
            XChemMain.print_acedrg_status(self.xce_logfile, self.xtal_db_dict)
            XChemMain.print_cluster_status_message(
                "acedrg", cluster_dict, self.xce_logfile
            )

        elif instruction.startswith("Run xia2 on selected datasets"):
            XChemMain.print_cluster_status_message(
                "xia2", cluster_dict, self.xce_logfile
            )

        elif "pandda" in instruction.lower():
            XChemMain.print_cluster_status_message(
                "pandda", cluster_dict, self.xce_logfile
            )

        elif "coot" in instruction.lower():
            XChemMain.print_cluster_status_message(
                "refmac", cluster_dict, self.xce_logfile
            )

    def prepare_and_run_task(self, instruction):
        if instruction == "Get New Results from Autoprocessing":
            self.rescore = False
            self.check_for_new_autoprocessing_results()

        elif instruction == "Rescore Datasets":
            self.rescore = True
            self.select_best_autoprocessing_result()

        elif instruction == "Run DIMPLE on selected MTZ files":
            self.run_dimple_on_selected_autoprocessing_file(instruction)

        elif instruction == "Run PIPEDREAM on selected MTZ files":
            self.run_dimple_on_selected_autoprocessing_file(instruction)

        elif instruction == "Run PHENIX.LIGAND_PIPELINE on selected MTZ files":
            self.run_dimple_on_selected_autoprocessing_file(instruction)

        elif instruction == "Remove selected DIMPLE files":
            self.remove_selected_dimple_files(instruction)

        elif instruction == "Remove selected PIPEDREAM files":
            self.remove_selected_dimple_files(instruction)

        elif instruction == "Remove selected PHENIX.LIGAND_PIPELINE files":
            self.remove_selected_dimple_files(instruction)

        elif instruction == "Set DIMPLE output":
            self.set_results_from_selected_pipeline(instruction)

        elif instruction == "Set PIPEDREAM output":
            self.set_results_from_selected_pipeline(instruction)

        elif instruction == "Set PHENIX.LIGAND_PIPELINE output":
            self.set_results_from_selected_pipeline(instruction)

        elif instruction == "Create CIF/PDB/PNG file of SELECTED compounds":
            self.create_cif_pdb_png_files("SELECTED")

        elif instruction == "Merge ligand CIF file with selected compounds":
            self.merge_cif_files("merge")

        elif instruction == "Restore original CIF file of selected compounds":
            self.merge_cif_files("restore")

        elif instruction == "Fit ligands into maps after initial refinement":
            self.fit_ligands_into_dimple_maps()

        elif instruction == "pandda.analyse":
            self.run_pandda_analyse("production_run")

        elif instruction == "pandda.analyse (PanDDA2)":
            self.run_pandda_analyse("production_run_pandda_two")

        elif instruction == "pre-run for ground state model":
            self.run_pandda_analyse("pre_run")

        elif instruction == "pandda.inspect":
            self.run_pandda_inspect()

        elif instruction == "run pandda.inspect at home":
            self.run_pandda_inspect_at_home()

        elif instruction == "Export NEW PANDDA models":
            update_datasource_only = False
            which_models = "new"
            self.run_pandda_export(update_datasource_only, which_models)

        elif instruction == "Export ALL PANDDA models":
            update_datasource_only = False
            which_models = "all"
            self.run_pandda_export(update_datasource_only, which_models)

        elif instruction == "Export SELECTED PANDDA models":
            update_datasource_only = False
            which_models = "selected"
            self.run_pandda_export(update_datasource_only, which_models)

        elif instruction == "refine ALL bound-state models with BUSTER":
            self.run_refine_bound_state_with_buster("all")

        elif instruction == "refine NEW bound-state models with BUSTER":
            self.run_refine_bound_state_with_buster("new")

        elif (
            instruction == "refine ALL bound-state models with BUSTER (no sanity check)"
        ):
            self.run_refine_bound_state_with_buster("allnocheck")

        elif (
            instruction == "refine NEW bound-state models with BUSTER (no sanity check)"
        ):
            self.run_refine_bound_state_with_buster("newnocheck")

        elif instruction == "cluster datasets":
            self.cluster_datasets_for_pandda()

        elif instruction == "Update datasource with results from pandda.inspect":
            update_datasource_only = True
            which_models = "all"
            self.run_pandda_export(update_datasource_only, which_models)

        elif instruction == "Show HTML summary":
            self.show_pandda_html_summary()

        elif instruction == "Event Map -> SF":
            self.convert_event_maps_to_SF()

        elif instruction == "apo -> mmcif":
            self.convert_apo_to_mmcif()

        elif instruction == "check modelled ligands":
            self.compare_modelled_ligands_and_panddaTable()

        elif (
            instruction.startswith("Open COOT")
            or instruction == "Build ground state model"
        ):
            if not self.coot_running:
                self.update_log.insert("starting coot...")
                if instruction == "Open COOT":
                    interface = "new"
                elif instruction == "Open COOT - REFMAC refinement -":
                    interface = "new"
                elif instruction == "Open COOT - test -":
                    interface = "test"
                elif instruction == "Open COOT for old PanDDA":
                    interface = "panddaV1"
                elif instruction == "Build ground state model":
                    interface = "reference"
                elif instruction == "Open COOT - BUSTER refinement -":
                    interface = "buster"
                elif instruction == "Open COOT - dimple_twin -":
                    interface = "dimple_twin"
                else:
                    interface = "old"
                self.work_thread = XChemThread.start_COOT(self.settings, interface)
                self.connect(
                    self.work_thread, QtCore.SIGNAL("finished()"), self.thread_finished
                )
                self.work_thread.start()

        elif instruction == "Update Deposition Table":
            self.update_deposition_table()

    def run_pandda_analyse(self, run):
        pandda_params = {
            "data_dir": str(self.pandda_input_data_dir_entry.text()),
            "out_dir": str(self.pandda_output_data_dir_entry.text()),
            "submit_mode": str(
                self.pandda_submission_mode_selection_combobox.currentText()
            ),
            "nproc": str(self.pandda_nproc_entry.text()),
            "min_build_datasets": str(self.pandda_min_build_dataset_entry.text()),
            "pdb_style": str(self.pandda_pdb_style_entry.text()),
            "mtz_style": str(self.pandda_mtz_style_entry.text()),
            "sort_event": str(self.pandda_sort_event_combobox.currentText()),
            "average_map": str(self.pandda_calc_map_combobox.currentText()),
            "max_new_datasets": str(self.pandda_max_new_datasets_entry.text()),
            "grid_spacing": str(self.pandda_grid_spacing_entry.text()),
            "keyword_arguments": str(self.pandda_keyword_arguments_entry.text()),
            "pandda_dir_structure": str(self.pandda_input_data_dir_entry.text()),
            "perform_diffraction_data_scaling": str(self.wilson_checkbox.isChecked()),
            "filter_pdb": str(
                self.pandda_reference_file_selection_combobox.currentText()
            ),
            "reference_dir": self.reference_directory,
            "appendix": "",
            "N_datasets": len(
                glob.glob(os.path.join(self.initial_model_directory, "*", "dimple.pdb"))
            ),
            "write_mean_map": "interesting",
            "pandda_table": self.pandda_analyse_data_table,
        }

        if run == "pre_run":
            msgBox = QtGui.QMessageBox()
            msgBoxLayout = msgBox.layout()
            vbox = QtGui.QVBoxLayout()
            vbox.addWidget(
                QtGui.QLabel(XChemToolTips.pandda_pre_run(self.reference_directory))
            )
            hbox = QtGui.QHBoxLayout()
            hbox.addWidget(QtGui.QLabel("appendix:"))
            appendix = QtGui.QLineEdit()
            appendix.setText("pre")
            appendix.setFixedWidth(200)
            hbox.addWidget(appendix)
            vbox.addLayout(hbox)

            msgBoxLayout.addLayout(vbox, 0, 0)
            msgBox.addButton(QtGui.QPushButton("Go"), QtGui.QMessageBox.YesRole)
            msgBox.addButton(QtGui.QPushButton("Cancel"), QtGui.QMessageBox.RejectRole)
            reply = msgBox.exec_()
            if reply == 0:
                pandda_params["appendix"] = str(appendix.text())
                pandda_params["max_new_datasets"] = "100"
                pandda_params["N_datasets"] = 100
                pandda_params["write_mean_map"] = "all"
            else:
                return None

        self.update_log.insert("preparing pandda.analyse input script")
        if run == "production_run_pandda_two":
            self.work_thread = XChemPANDDA.run_pandda_two_analyse(
                pandda_params,
                self.xce_logfile,
                os.path.join(self.database_directory, self.data_source_file),
                self.external_software,
                get_token(),
            )
        else:
            self.work_thread = XChemPANDDA.run_pandda_analyse(
                pandda_params,
                self.xce_logfile,
                os.path.join(self.database_directory, self.data_source_file),
                self.external_software,
                get_token(),
            )
        self.connect(
            self.work_thread, QtCore.SIGNAL("finished()"), self.thread_finished
        )
        self.work_thread.start()

    def cluster_datasets_for_pandda(self):
        pandda_params = {
            "out_dir": str(self.pandda_output_data_dir_entry.text()),
            "pdb_style": str(self.pandda_pdb_style_entry.text()),
            "mtz_style": str(self.pandda_mtz_style_entry.text()),
        }
        self.update_log.insert("starting giant.cluster_mtzs_and_pdbs")
        self.work_thread = XChemPANDDA.giant_cluster_datasets(
            self.initial_model_directory,
            pandda_params,
            self.xce_logfile,
            os.path.join(self.database_directory, self.data_source_file),
            self.run_pandda_analyse,
        )
        self.explorer_active = 1
        self.connect(
            self.work_thread,
            QtCore.SIGNAL("update_progress_bar"),
            self.update_progress_bar,
        )
        self.connect(
            self.work_thread,
            QtCore.SIGNAL("update_status_bar(QString)"),
            self.update_status_bar,
        )
        self.connect(
            self.work_thread,
            QtCore.SIGNAL("datasource_menu_reload_samples"),
            self.datasource_menu_reload_samples,
        )
        self.connect(
            self.work_thread, QtCore.SIGNAL("finished()"), self.thread_finished
        )
        self.work_thread.start()

    def run_pandda_inspect(self):
        self.settings["panddas_directory"] = str(
            self.pandda_output_data_dir_entry.text()
        )
        print("==> XCE: starting pandda.inspect")
        self.work_thread = XChemThread.start_pandda_inspect(
            self.settings, self.xce_logfile
        )
        self.connect(
            self.work_thread, QtCore.SIGNAL("finished()"), self.thread_finished
        )
        self.work_thread.start()

    def run_pandda_inspect_at_home(self):
        self.work_thread = XChemPANDDA.run_pandda_inspect_at_home(
            self.panddas_directory, self.xce_logfile
        )
        self.connect(
            self.work_thread, QtCore.SIGNAL("finished()"), self.thread_finished
        )
        self.work_thread.start()
        self.connect(
            self.work_thread,
            QtCore.SIGNAL("update_progress_bar"),
            self.update_progress_bar,
        )
        self.connect(
            self.work_thread,
            QtCore.SIGNAL("update_status_bar(QString)"),
            self.update_status_bar,
        )
        self.connect(
            self.work_thread, QtCore.SIGNAL("finished()"), self.thread_finished
        )

    def convert_event_maps_to_SF(self):
        self.update_log.insert(
            "converting all event maps in {0!s} to mtz files".format(
                self.initial_model_directory
            )
        )
        self.work_thread = XChemPANDDA.find_event_map_for_ligand(
            self.initial_model_directory, self.xce_logfile, self.external_software
        )

        self.explorer_active = 1
        self.connect(
            self.work_thread,
            QtCore.SIGNAL("update_progress_bar"),
            self.update_progress_bar,
        )
        self.connect(
            self.work_thread,
            QtCore.SIGNAL("update_status_bar(QString)"),
            self.update_status_bar,
        )
        self.connect(
            self.work_thread, QtCore.SIGNAL("finished()"), self.thread_finished
        )
        self.work_thread.start()

    def convert_apo_to_mmcif(self):
        self.work_thread = XChemPANDDA.convert_apo_structures_to_mmcif(
            self.panddas_directory, self.xce_logfile
        )

        self.explorer_active = 1
        self.connect(
            self.work_thread,
            QtCore.SIGNAL("update_progress_bar"),
            self.update_progress_bar,
        )
        self.connect(
            self.work_thread,
            QtCore.SIGNAL("update_status_bar(QString)"),
            self.update_status_bar,
        )
        self.connect(
            self.work_thread, QtCore.SIGNAL("finished()"), self.thread_finished
        )
        self.work_thread.start()

    def compare_modelled_ligands_and_panddaTable(self):
        self.update_log.insert(
            "checking agreement of ligands in refine.pdb and entries in panddaTable"
        )
        self.work_thread = XChemPANDDA.check_number_of_modelled_ligands(
            self.initial_model_directory,
            self.xce_logfile,
            os.path.join(self.database_directory, self.data_source_file),
        )
        self.explorer_active = 1
        self.connect(
            self.work_thread,
            QtCore.SIGNAL("update_progress_bar"),
            self.update_progress_bar,
        )
        self.connect(
            self.work_thread,
            QtCore.SIGNAL("update_status_bar(QString)"),
            self.update_status_bar,
        )
        self.connect(
            self.work_thread, QtCore.SIGNAL("finished()"), self.thread_finished
        )
        self.connect(
            self.work_thread, QtCore.SIGNAL("show_error_dict"), self.show_error_dict
        )
        self.work_thread.start()

    def run_pandda_export(self, update_datasource_only, which_models):
        pandda_params = {
            "data_dir": str(self.pandda_input_data_dir_entry.text()),
            "out_dir": str(self.pandda_output_data_dir_entry.text()),
            "submit_mode": str(
                self.pandda_submission_mode_selection_combobox.currentText()
            ),
            "nproc": str(self.pandda_nproc_entry.text()),
            "min_build_datasets": str(self.pandda_min_build_dataset_entry.text()),
            "pdb_style": str(self.pandda_pdb_style_entry.text()),
            "mtz_style": str(self.pandda_mtz_style_entry.text()),
            "sort_event": str(self.pandda_sort_event_combobox.currentText()),
            "average_map": str(self.pandda_calc_map_combobox.currentText()),
            "max_new_datasets": str(self.pandda_max_new_datasets_entry.text()),
            "grid_spacing": str(self.pandda_grid_spacing_entry.text()),
            "pandda_dir_structure": str(self.pandda_input_data_dir_entry.text()),
            "perform_diffraction_data_scaling": str(self.wilson_checkbox.isChecked()),
            "filter_pdb": str(
                self.pandda_reference_file_selection_combobox.currentText()
            ),
            "reference_dir": self.reference_directory,
            "appendix": "",
            "N_datasets": len(
                glob.glob(os.path.join(self.initial_model_directory, "*", "dimple.pdb"))
            ),
            "write_mean_map": "interesting",
            "pandda_table": self.pandda_analyse_data_table,
        }

        self.settings["panddas_directory"] = str(
            self.pandda_output_data_dir_entry.text()
        )
        if update_datasource_only:
            self.update_log.insert(
                "updating data source with results from pandda.inspect"
            )
        else:
            self.update_log.insert(
                "exporting PANDDA models,"
                " updating data source and launching inital refinement for new models"
            )

        start_thread = False
        if which_models == "all":
            self.update_log.insert(
                "exporting ALL models!"
                " *** WARNING *** This may overwrite previous refinements!!!"
            )
            msgBox = QtGui.QMessageBox()
            msgBox.setText(
                "*** WARNING ***\n"
                "This will overwrite all your manual selections!\n"
                "Do you want to continue?"
            )
            msgBox.addButton(QtGui.QPushButton("Yes"), QtGui.QMessageBox.YesRole)
            msgBox.addButton(QtGui.QPushButton("No"), QtGui.QMessageBox.RejectRole)
            reply = msgBox.exec_()
            if reply == 0:
                if update_datasource_only:
                    self.update_log.insert("will update panddaTable in database only")
                else:
                    self.update_log.insert("will export ALL models!")
                start_thread = True
            else:
                start_thread = False
        else:
            self.update_log.insert("exporting new models only")
            start_thread = True

        if start_thread:
            self.work_thread = XChemPANDDA.run_pandda_export(
                self.panddas_directory,
                os.path.join(self.database_directory, self.data_source_file),
                self.initial_model_directory,
                self.xce_logfile,
                which_models,
                pandda_params,
                get_token(),
            )
            self.connect(
                self.work_thread, QtCore.SIGNAL("finished()"), self.thread_finished
            )
            self.work_thread.start()

    def run_refine_bound_state_with_buster(self, which_models):
        start_thread = True
        if start_thread:
            self.work_thread = XChemPANDDA.export_and_refine_ligand_bound_models(
                self.panddas_directory,
                os.path.join(self.database_directory, self.data_source_file),
                self.initial_model_directory,
                self.xce_logfile,
                which_models,
                get_token(),
            )
            self.connect(
                self.work_thread, QtCore.SIGNAL("finished()"), self.thread_finished
            )
            self.work_thread.start()

    def show_pandda_html_summary(self):
        self.pandda_initial_html.load(QtCore.QUrl(self.pandda_initial_html_file))
        self.pandda_initial_html.show()
        self.pandda_analyse_html.load(QtCore.QUrl(self.pandda_analyse_html_file))
        self.pandda_analyse_html.show()
        self.add_map_html()
        self.pandda_inspect_html.load(QtCore.QUrl(self.pandda_inspect_html_file))
        self.pandda_inspect_html.show()

    def create_cif_pdb_png_files(self, todo):
        tmp = self.db.execute_statement(
            "select CrystalName,CompoundCode,CompoundSmiles from mainTable where"
            " CrystalName is not '' and CompoundSmiles is not '' and CompoundSmiles"
            " is not NULL;"
        )
        compound_list = []
        for item in tmp:
            if str(item[1]) == "" or str(item[1]) == "NULL":
                compoundID = "compound"
            else:
                compoundID = str(item[1])

            if todo == "ALL":
                compound_list.append([str(item[0]), compoundID, str(item[2])])
            elif todo == "NEW":
                if not os.path.isfile(
                    os.path.join(
                        self.initial_model_directory, str(item[0]), compoundID + ".cif"
                    )
                ):
                    compound_list.append([str(item[0]), compoundID, str(item[2])])
            elif todo == "SELECTED":
                if str(item[0]) in self.initial_model_dimple_dict:
                    if self.initial_model_dimple_dict[str(item[0])][0].isChecked():
                        compound_list.append([str(item[0]), compoundID, str(item[2])])

        if compound_list:
            self.update_log.insert(
                "trying to create cif and pdb files for "
                + str(len(compound_list))
                + " compounds using ACEDRG..."
            )
            self.explorer_active = 1
            self.work_thread = XChemThread.create_png_and_cif_of_compound(
                self.external_software,
                self.initial_model_directory,
                compound_list,
                self.database_directory,
                self.data_source_file,
                todo,
                self.ccp4_scratch_directory,
                self.xce_logfile,
                self.max_queue_jobs,
                self.restraints_program,
                get_token(),
            )
            self.connect(
                self.work_thread,
                QtCore.SIGNAL("update_progress_bar"),
                self.update_progress_bar,
            )
            self.connect(
                self.work_thread,
                QtCore.SIGNAL("update_status_bar(QString)"),
                self.update_status_bar,
            )
            self.connect(
                self.work_thread, QtCore.SIGNAL("finished()"), self.thread_finished
            )
            self.connect(
                self.work_thread,
                QtCore.SIGNAL("datasource_menu_reload_samples"),
                self.datasource_menu_reload_samples,
            )
            self.work_thread.start()

    def fit_ligands_into_dimple_maps(self):
        tmp = self.db.execute_statement(
            "select CrystalName,CompoundCode,CompoundSmiles from mainTable where"
            " CrystalName is not '' and CompoundSmiles is not '' and CompoundSmiles"
            " is not NULL;"
        )
        compound_list = []
        for item in tmp:
            if str(item[1]) == "" or str(item[1]) == "NULL":
                compoundID = "compound"
            else:
                compoundID = str(item[1])

            if str(item[0]) in self.initial_model_dimple_dict:
                if self.initial_model_dimple_dict[str(item[0])][0].isChecked():
                    compound_list.append([str(item[0]), compoundID, str(item[2])])

        if compound_list:
            self.update_log.insert(
                "trying to auto-fitting into inital maps for "
                + str(len(compound_list))
                + " compounds..."
            )
            self.explorer_active = 1
            self.work_thread = XChemThread.fit_ligands(
                self.external_software,
                self.initial_model_directory,
                compound_list,
                self.database_directory,
                self.data_source_file,
                self.ccp4_scratch_directory,
                self.xce_logfile,
                self.max_queue_jobs,
                get_token(),
            )
            self.connect(
                self.work_thread,
                QtCore.SIGNAL("update_progress_bar"),
                self.update_progress_bar,
            )
            self.connect(
                self.work_thread,
                QtCore.SIGNAL("update_status_bar(QString)"),
                self.update_status_bar,
            )
            self.connect(
                self.work_thread, QtCore.SIGNAL("finished()"), self.thread_finished
            )
            self.connect(
                self.work_thread,
                QtCore.SIGNAL("datasource_menu_reload_samples"),
                self.datasource_menu_reload_samples,
            )
            self.work_thread.start()

    def merge_cif_files(self, todo):
        start_thread = False
        if todo == "merge":
            self.update_log.insert(
                "trying to merge %s with ligand restraint files in project directory"
                % self.second_cif_file
            )
        elif todo == "restore":
            self.update_log.insert("restoring original CIF files")
            start_thread = True

        if todo == "merge":
            if os.path.isfile(str(self.second_cif_file)):
                self.update_log.insert(
                    "checking compound code of second CIF file (%s)"
                    % self.second_cif_file
                )
                self.update_log.insert("Note: LIG and DRG are not allowed!")
                import iotbx.cif

                cif_model = iotbx.cif.reader(file_path=self.second_cif_file).model()
                cif_block = cif_model["comp_list"]
                ligID = cif_block["_chem_comp.id"]
                self.update_log.insert(
                    "found the following compound codes in the supplied CIF file: %s"
                    % str(list(ligID))
                )
                if "LIG" in list(ligID) or "DRG" in list(ligID):
                    self.update_log.error(
                        "please change compound code to something other than LIG or DRG"
                    )
                    start_thread = False
                else:
                    start_thread = True
            else:
                self.update_log.error(XChemToolTips.second_cif_file_not_exists())
                start_thread = False

            if start_thread:
                msgBox = QtGui.QMessageBox()
                msgBox.setText(XChemToolTips.second_cif_file_info(self.second_cif_file))
                msgBox.addButton(QtGui.QPushButton("OK"), QtGui.QMessageBox.YesRole)
                msgBox.addButton(
                    QtGui.QPushButton("Cancel"), QtGui.QMessageBox.RejectRole
                )
                reply = msgBox.exec_()
                if reply == 0:
                    start_thread = True
                else:
                    start_thread = False
            else:
                self.status_bar.showMessage(
                    "Error. Please check terminal window for further information"
                )

        tmp = self.db.execute_statement(
            "select CrystalName,CompoundCode from mainTable where CrystalName is not ''"
            " and CompoundSmiles is not '' and CompoundSmiles is not NULL;"
        )
        compound_list = []
        for item in tmp:
            xtal = str(item[0])
            compoundID = str(item[1])
            if compoundID == "" or compoundID == "NULL":
                self.update_log.warning(
                    "%s: no compound ID in database; skipping..." % xtal
                )
            else:
                if str(item[0]) in self.initial_model_dimple_dict:
                    if self.initial_model_dimple_dict[str(item[0])][0].isChecked():
                        self.update_log.warning(
                            "%s: %s is flagged for merging" % (xtal, compoundID)
                        )
                        compound_list.append([xtal, compoundID])

        if compound_list == []:
            self.update_log.error(
                "Either no compound ID information in database or no sample selected!"
            )
            start_thread = False

        if start_thread:
            self.explorer_active = 1
            self.work_thread = XChemThread.merge_cif_files(
                self.initial_model_directory,
                self.xce_logfile,
                self.second_cif_file,
                compound_list,
                todo,
            )
            self.connect(
                self.work_thread,
                QtCore.SIGNAL("update_progress_bar"),
                self.update_progress_bar,
            )
            self.connect(
                self.work_thread,
                QtCore.SIGNAL("update_status_bar(QString)"),
                self.update_status_bar,
            )
            self.connect(
                self.work_thread, QtCore.SIGNAL("finished()"), self.thread_finished
            )
            self.connect(
                self.work_thread,
                QtCore.SIGNAL("datasource_menu_reload_samples"),
                self.datasource_menu_reload_samples,
            )
            self.work_thread.start()

    def update_deposition_table(self):
        # check if PanDDA models are ready for deposition

        depositChecks = XChemDeposit.update_deposition_table(
            os.path.join(self.database_directory, self.data_source_file)
        )

        toDeposit, mismatch = depositChecks.PanDDA_models_to_deposit()

        if mismatch != {}:
            self.update_log.insert(
                "The following samples contain ligand that are not ready for"
                " deposition:"
            )
            for entry in mismatch:
                self.update_log.insert(
                    entry[0]
                    + " -> site: "
                    + entry[1]
                    + " @ "
                    + entry[2]
                    + " => "
                    + entry[4]
                )
            self.update_log.insert("You need to change this before you can continue!")
            return None

        for xtal in toDeposit:
            self.db.update_insert_depositTable(xtal, {})

    def need_to_switch_main_tab(self, task_index):
        msgBox = QtGui.QMessageBox()
        msgBox.setText("Need to switch main tab before you can launch this job")
        msgBox.addButton(QtGui.QPushButton("Yes"), QtGui.QMessageBox.YesRole)
        msgBox.addButton(QtGui.QPushButton("No"), QtGui.QMessageBox.RejectRole)
        reply = msgBox.exec_()
        if reply == 0:
            self.main_tab_widget.setCurrentIndex(task_index)

    def check_write_permissions_of_data_source(self):
        write_enabled = True
        if not os.access(
            os.path.join(self.database_directory, self.data_source_file), os.W_OK
        ):
            QtGui.QMessageBox.warning(
                self.window,
                "Data Source Problem",
                "\nData Source is Read-Only\n",
                QtGui.QMessageBox.Cancel,
                QtGui.QMessageBox.NoButton,
                QtGui.QMessageBox.NoButton,
            )
            write_enabled = False
        return write_enabled

    def no_data_source_selected(self):
        QtGui.QMessageBox.warning(
            self.window,
            "Data Source Problem",
            ("Please set or create a data source file\n")
            + ("Options:\n")
            + ("1. Use an existing file:\n")
            + ("- Settings -> Select Data Source File\n")
            + ("2. Create a new file\n")
            + ("- Data Source -> Create New Data\nSource (SQLite)"),
            QtGui.QMessageBox.Cancel,
            QtGui.QMessageBox.NoButton,
            QtGui.QMessageBox.NoButton,
        )

    def update_progress_bar(self, progress):
        self.progress_bar.setValue(progress)

    def update_status_bar(self, message):
        self.status_bar.showMessage(message)

    def thread_finished(self):
        self.explorer_active = 0
        self.update_progress_bar(0)
        self.update_status_bar("idle")

    def show_error_dict(self, errorDict):
        text = ""
        for key in errorDict:
            text += "{0!s}:\n".format(key)
            for entry in errorDict[key]:
                text += "  - " + entry + "\n"
        msgBox = QtGui.QMessageBox()
        msgBox.setText(text)
        msgBox.exec_()

    def create_widgets_for_autoprocessing_results_only(self, data_dict):
        self.status_bar.showMessage(
            "Building details table for data processing results"
        )
        self.data_collection_dict = data_dict

        column_name = [
            "Program",
            "Resolution\nOverall",
            "Resolution\n[Mn<I/sig(I)> = 2.0]",
            "DataProcessing\nSpaceGroup",
            "Mn<I/sig(I)>\nHigh",
            "Rmerge\nLow",
            "Completeness\nOverall",
            "DataProcessing\nUnitCell",
            "DataProcessing\nRfree",
            "DataProcessing\nScore",
        ]

        # need to do this because db_dict keys are SQLite column names
        diffraction_data_column_name = XChemDB.data_source(
            os.path.join(self.database_directory, self.data_source_file)
        ).translate_xce_column_list_to_sqlite(column_name)

        for xtal in sorted(self.data_collection_dict):
            # column 2: data collection date
            # this one should always be there; it may need updating in case another run
            # appears first find latest run
            tmp = []
            for entry in self.data_collection_dict[xtal]:
                if entry[0] == "image":
                    tmp.append(
                        [entry[3], datetime.strptime(entry[3], "%Y-%m-%d %H:%M:%S")]
                    )

            # first check if it does already exist
            if xtal not in self.data_collection_column_three_dict:
                # generate all the widgets which can later be appended and add them to
                # the dictionary table with data processing results for each pipeline
                data_collection_table = QtGui.QTableWidget()
                selection_changed_by_user = False
                self.data_collection_column_three_dict[xtal] = [
                    data_collection_table,
                    selection_changed_by_user,
                ]
            else:
                data_collection_table = self.data_collection_column_three_dict[xtal][0]
                selection_changed_by_user = self.data_collection_column_three_dict[
                    xtal
                ][1]

            data_collection_table.setVerticalScrollBarPolicy(
                QtCore.Qt.ScrollBarAlwaysOff
            )
            data_collection_table.setColumnCount(len(column_name))
            font = QtGui.QFont()
            font.setPointSize(8)
            data_collection_table.setFont(font)
            data_collection_table.setHorizontalHeaderLabels(column_name)
            data_collection_table.horizontalHeader().setFont(font)
            data_collection_table.setSelectionBehavior(
                QtGui.QAbstractItemView.SelectRows
            )

            ############################################################################
            # crystal images
            # first check there are new images that are not displayed yet; i.e. they
            # are not in the self.data_collection_image_dict
            if xtal not in self.data_collection_image_dict:
                # OK this is the first time
                self.data_collection_image_dict[xtal] = []

            # sort crystal images by timestamp
            # a) get only image entries from self.data_collection_dict
            tmp = []
            for entry in self.data_collection_dict[xtal]:
                if entry[0] == "image":
                    tmp.append(entry)

            # b) sort by the previously assigned run number
            for entry in sorted(tmp, key=lambda x: x[6]):
                run_number = entry[6]
                images_already_in_table = False
                for image in self.data_collection_image_dict[xtal]:
                    if run_number == image[0]:
                        images_already_in_table = True
                        break
                if not images_already_in_table:
                    # not if there is a run, but images are for whatever reason not
                    # present in self.data_collection_dict then use image not available
                    # from $XChemExplorer_DIR/image/IMAGE_NOT_AVAILABLE.png
                    self.data_collection_image_dict[xtal].append(
                        [entry[6], entry[1], entry[2], entry[3], entry[5]]
                    )

            ############################################################################
            # initialize dataset_outcome_dict for xtal
            if xtal not in self.dataset_outcome_dict:
                self.dataset_outcome_dict[xtal] = []
                # dataset outcome buttons

            ############################################################################
            # table for data processing results
            # check if results from particular pipeline are already in table
            row_position = data_collection_table.rowCount()
            if xtal not in self.data_collection_table_dict:
                self.data_collection_table_dict[xtal] = []
            logfile_list = []
            for entry in self.data_collection_dict[xtal]:
                if entry[0] == "logfile":
                    logfile_list.append(entry)
            # sort by aimless_index and so make sure
            for entry in sorted(logfile_list, key=lambda x: x[7]):
                entry_already_in_table = False  # that aimless_index == row
                cell_text = QtGui.QTableWidgetItem()
                cell_text.setTextAlignment(
                    QtCore.Qt.AlignCenter | QtCore.Qt.AlignCenter
                )
                db_dict = entry[6]
                for logfile in self.data_collection_table_dict[xtal]:
                    if (
                        entry[1] == logfile[1]
                        and entry[2] == logfile[2]
                        and entry[3] == logfile[3]
                        and entry[4] == logfile[4]
                    ):
                        entry_already_in_table = True
                        # might have to update Rfree column
                        for column, header in enumerate(diffraction_data_column_name):
                            if header == "DataProcessing\nRfree":
                                # entry[7]==aimless_index, i.e. row number
                                cell_text.setText(str(db_dict[header[1]]))
                                data_collection_table.setItem(
                                    entry[7], column, cell_text
                                )
                                break
                        break
                if not entry_already_in_table:
                    data_collection_table.insertRow(row_position)
                    for column, header in enumerate(diffraction_data_column_name):
                        cell_text = QtGui.QTableWidgetItem()
                        try:
                            cell_text.setText(str(db_dict[header[1]]))
                        except KeyError:
                            # this may happen if not score exists
                            cell_text.setText("0")
                        data_collection_table.setItem(row_position, column, cell_text)
                    data_collection_table.setRowHeight(row_position, 20)
                    row_position += 1

                self.data_collection_table_dict[xtal].append(
                    ["logfile", entry[1], entry[2], entry[3], entry[4]]
                )  # 'logfile' is just added to have
                # same index numbers between lists
            data_collection_table.cellClicked.connect(
                self.user_update_selected_autoproc_datasets_summary_table
            )

            # select best resolution file + set data collection outcome
            # the assumption is that index in data_collection_dict and row number are
            # identical
            # the assumption for data collection outcome is that as long as a logfile
            # is found, it's a success
            for entry in self.data_collection_dict[xtal]:
                if entry[0] == "logfile":
                    index = entry[7]
                    best_file = entry[8]
                    if best_file:
                        # we change the selection only if the user did not touch it,
                        # assuming that he/she knows best if not
                        # selection_changed_by_user:
                        data_collection_table.selectRow(index)

        self.populate_datasets_summary_table()

    def find_suitable_reference_file(self, db_dict):
        reference_file = []
        dummy = ["...", "", "", "", 0, "0"]
        reference_file.append([dummy, 999])
        for reference in self.reference_file_list:
            # first we need one in the same pointgroup
            if reference[5] == db_dict["DataProcessingPointGroup"]:
                try:
                    difference = (
                        math.fabs(
                            1
                            - (
                                float(db_dict["DataProcessingUnitCellVolume"])
                                / float(reference[4])
                            )
                        )
                        * 100
                    )
                    reference_file.append([reference, difference])
                except ValueError:
                    continue
        return reference_file

    def create_maps_table(self):
        column_name = self.db.translate_xce_column_list_to_sqlite(
            self.maps_table_columns
        )

        for xtal in sorted(self.xtal_db_dict):
            new_xtal = False
            db_dict = self.xtal_db_dict[xtal]
            if str(db_dict["DataCollectionOutcome"]).lower().startswith("success"):
                reference_file = self.find_suitable_reference_file(db_dict)
                smallest_uc_difference = min(reference_file, key=lambda x: x[1])
                row = self.maps_table.rowCount()
                if xtal not in self.initial_model_dimple_dict:
                    self.maps_table.insertRow(row)
                    current_row = row
                    new_xtal = True
                else:
                    for table_row in range(row):
                        if self.maps_table.item(table_row, 0).text() == xtal:
                            current_row = table_row
                            break
                for column, header in enumerate(column_name):
                    if header[0] == "Sample ID":
                        cell_text = QtGui.QTableWidgetItem()
                        cell_text.setText(str(xtal))
                        cell_text.setTextAlignment(
                            QtCore.Qt.AlignCenter | QtCore.Qt.AlignCenter
                        )
                        self.maps_table.setItem(current_row, column, cell_text)
                    elif header[0] == "Select":
                        if new_xtal:
                            run_dimple = QtGui.QCheckBox()
                            run_dimple.toggle()
                            self.maps_table.setCellWidget(
                                current_row, column, run_dimple
                            )
                            run_dimple.setChecked(False)
                    elif header[0] == "Reference\nSpaceGroup":
                        cell_text = QtGui.QTableWidgetItem()
                        cell_text.setText(str(smallest_uc_difference[0][1]))
                        cell_text.setTextAlignment(
                            QtCore.Qt.AlignCenter | QtCore.Qt.AlignCenter
                        )
                        self.maps_table.setItem(current_row, column, cell_text)
                    elif header[0] == "Difference\nUC Volume (%)":
                        cell_text = QtGui.QTableWidgetItem()
                        smallest_uc_difference = min(reference_file, key=lambda x: x[1])
                        cell_text.setText(
                            str(round(float(smallest_uc_difference[1]), 1))
                        )
                        cell_text.setTextAlignment(
                            QtCore.Qt.AlignCenter | QtCore.Qt.AlignCenter
                        )
                        self.maps_table.setItem(current_row, column, cell_text)
                    elif header[0] == "Reference File":
                        if new_xtal:
                            reference_file_selection_combobox = QtGui.QComboBox()
                            self.populate_reference_combobox(
                                reference_file_selection_combobox
                            )
                            if (
                                float(smallest_uc_difference[1])
                                < self.allowed_unitcell_difference_percent
                            ):
                                index = reference_file_selection_combobox.findText(
                                    str(smallest_uc_difference[0][0]),
                                    QtCore.Qt.MatchFixedString,
                                )
                                reference_file_selection_combobox.setCurrentIndex(index)
                            else:
                                reference_file_selection_combobox.setCurrentIndex(0)
                            self.maps_table.setCellWidget(
                                current_row, column, reference_file_selection_combobox
                            )
                        else:
                            reference_file_selection_combobox = (
                                self.initial_model_dimple_dict[xtal][1]
                            )
                            self.populate_reference_combobox(
                                reference_file_selection_combobox
                            )
                            if (
                                float(smallest_uc_difference[1])
                                < self.allowed_unitcell_difference_percent
                            ):
                                index = reference_file_selection_combobox.findText(
                                    str(smallest_uc_difference[0][0]),
                                    QtCore.Qt.MatchFixedString,
                                )
                                reference_file_selection_combobox.setCurrentIndex(index)
                            else:
                                reference_file_selection_combobox.setCurrentIndex(0)
                    else:
                        cell_text = QtGui.QTableWidgetItem()
                        cell_text.setText(str(db_dict[header[1]]))
                        cell_text.setTextAlignment(
                            QtCore.Qt.AlignCenter | QtCore.Qt.AlignCenter
                        )
                        if header[0] == "Dimple\nStatus":
                            if str(db_dict[header[1]]) == "running":
                                cell_text.setBackground(QtGui.QColor(100, 230, 150))
                            elif str(db_dict[header[1]]) == "pending":
                                cell_text.setBackground(QtGui.QColor(20, 100, 230))
                            elif str(db_dict[header[1]]) == "started":
                                cell_text.setBackground(QtGui.QColor(230, 240, 110))
                            elif str(db_dict[header[1]]) == "finished":
                                cell_text.setBackground(QtGui.QColor(255, 255, 255))
                        if header[0] == "Compound\nStatus":
                            if str(db_dict[header[1]]) == "running":
                                cell_text.setBackground(QtGui.QColor(100, 230, 150))
                            elif str(db_dict[header[1]]) == "pending":
                                cell_text.setBackground(QtGui.QColor(20, 100, 230))
                            elif str(db_dict[header[1]]) == "started":
                                cell_text.setBackground(QtGui.QColor(230, 240, 110))
                            elif str(db_dict[header[1]]) == "restraints generated":
                                cell_text.setBackground(QtGui.QColor(255, 255, 255))
                            elif str(db_dict[header[1]]) == "restraints failed":
                                cell_text.setBackground(QtGui.QColor(255, 0, 0))
                            elif str(db_dict[header[1]]) == "missing smiles":
                                cell_text.setBackground(QtGui.QColor(240, 150, 20))
                        self.maps_table.setItem(current_row, column, cell_text)
            if new_xtal:
                self.initial_model_dimple_dict[xtal] = [
                    run_dimple,
                    reference_file_selection_combobox,
                ]

    def preferences_data_to_copy_combobox_changed(self, i):
        text = str(self.preferences_data_to_copy_combobox.currentText())
        for item in self.preferences_data_to_copy:
            if item[0] == text:
                self.preferences["processed_data_to_copy"] = item[1]
                break

    def preferences_selection_mechanism_combobox_changed(self, i):
        text = str(self.preferences_selection_mechanism_combobox.currentText())
        self.preferences["dataset_selection_mechanism"] = text
        self.update_log.insert("setting datasets selection mechanism to " + text)

    def preferences_restraints_generation_combobox_changed(self):
        text = str(self.preferences_restraints_generation_combobox.currentText())
        self.restraints_program = text
        self.update_log.insert(
            "will use {0!s} for generation of ligand coordinates and restraints".format(
                text
            )
        )

    def refinement_outcome_combobox_changed(self):
        for xtal in self.refinement_table_dict:
            if self.sender() == self.refinement_table_dict[xtal]:
                db_dict = {}
                db_dict["RefinementOutcome"] = str(self.sender().currentText())
                db_dict["RefinementOutcomePerson"] = getpass.getuser()
                db_dict["RefinementOutcomeDate"] = datetime.strftime(
                    datetime.now(), "%Y-%m-%d_%H-%M-%S.%f"
                )[:-4]
                self.db.create_or_remove_missing_records_in_depositTable(
                    self.xce_logfile, xtal, "ligand_bound", db_dict
                )

    def get_reference_file_list(self, reference_root):
        # check available reference files
        reference_file_list = []
        dummy = ["...", "", "", "", 0, "0"]
        reference_file_list.append(dummy)
        if os.path.isfile(
            os.path.join(self.reference_directory, reference_root + ".pdb")
        ):
            pdb_reference = parse().PDBheader(
                os.path.join(self.reference_directory, reference_root + ".pdb")
            )
            spg_reference = pdb_reference["SpaceGroup"]
            unitcell_reference = pdb_reference["UnitCell"]
            lattice_reference = pdb_reference["Lattice"]
            unitcell_volume_reference = pdb_reference["UnitCellVolume"]
            pointgroup_reference = pdb_reference["PointGroup"]
            reference_file_list.append(
                [
                    reference_root,
                    spg_reference,
                    unitcell_reference,
                    lattice_reference,
                    unitcell_volume_reference,
                    pointgroup_reference,
                ]
            )
        else:
            for files in glob.glob(self.reference_directory + "/*"):
                if files.endswith(".pdb"):
                    reference_root = files[files.rfind("/") + 1 : files.rfind(".")]

                    if os.path.isfile(
                        os.path.join(self.reference_directory, reference_root + ".pdb")
                    ):
                        pdb_reference = parse().PDBheader(
                            os.path.join(
                                self.reference_directory, reference_root + ".pdb"
                            )
                        )
                        spg_reference = pdb_reference["SpaceGroup"]
                        unitcell_reference = pdb_reference["UnitCell"]
                        lattice_reference = pdb_reference["Lattice"]
                        unitcell_volume_reference = pdb_reference["UnitCellVolume"]
                        pointgroup_reference = pdb_reference["PointGroup"]
                        reference_file_list.append(
                            [
                                reference_root,
                                spg_reference,
                                unitcell_reference,
                                lattice_reference,
                                unitcell_volume_reference,
                                pointgroup_reference,
                            ]
                        )
        for n, file in enumerate(reference_file_list):
            self.update_log.insert("reference file {0!s}: {1!s}".format(n, file))
        return reference_file_list

    def dataset_outcome_combobox_change_outcome(self, text):
        outcome = str(text)
        xtal = ""
        for key in self.dataset_outcome_combobox_dict:
            if self.dataset_outcome_combobox_dict[key] == self.sender():
                xtal = key
                self.update_log.insert(
                    "user changed data collection outcome of {0!s} to {1!s}".format(
                        xtal, outcome
                    )
                )
                break
        self.dataset_outcome_dict[xtal] = outcome
        if xtal != "":
            # need to also update if not yet done
            self.update_log.insert(
                "updating dataset outcome in datasource for {0!s}".format(xtal)
            )
            update_dict = {"DataCollectionOutcome": outcome}
            self.db.update_insert_data_source(xtal, update_dict)

    def set_run_dimple_flag(self, state):
        if state == QtCore.Qt.Checked:
            for key in self.initial_model_dimple_dict:
                self.initial_model_dimple_dict[key][0].setChecked(True)
        else:
            for key in self.initial_model_dimple_dict:
                self.initial_model_dimple_dict[key][0].setChecked(False)

    ####################################################################################
    #
    #
    #
    # => new data collection summary table
    # > start

    def get_sample_list_from_table(self, table):
        sampleList = []
        allRows = table.rowCount()
        for row in range(0, allRows):
            sample_id = str(table.item(row, 0).text())
            sampleList.append(sample_id)
        return sorted(sampleList)

    def get_row_of_sample_in_table(self, table, xtal):
        allRows = table.rowCount()
        sampleRow = allRows
        for n, row in enumerate(range(0, allRows)):
            sample_id = str(table.item(row, 0).text())
            if sample_id == xtal:
                sampleRow = n
                break
        return sampleRow

    def update_row_in_table(self, sample, row, db_dict, table, columns_to_show):
        xtal = str(sample)
        column_name = self.db.translate_xce_column_list_to_sqlite(columns_to_show)

        for column, header in enumerate(column_name):
            if header[0] == "Sample ID":
                cell_text = QtGui.QTableWidgetItem()
                cell_text.setText(str(xtal))
                cell_text.setTextAlignment(
                    QtCore.Qt.AlignCenter | QtCore.Qt.AlignCenter
                )
                table.setItem(row, column, cell_text)

            elif header[0] == "DataCollection\nOutcome":
                if xtal not in self.dataset_outcome_combobox_dict:
                    dataset_outcome_combobox = QtGui.QComboBox()
                    for outcomeItem in self.dataset_outcome:
                        dataset_outcome_combobox.addItem(outcomeItem)
                    dataset_outcome_combobox.activated[str].connect(
                        self.dataset_outcome_combobox_change_outcome
                    )
                    self.dataset_outcome_combobox_dict[xtal] = dataset_outcome_combobox
                    table.setCellWidget(row, column, dataset_outcome_combobox)
                index = self.dataset_outcome_combobox_dict[xtal].findText(
                    str(db_dict["DataCollectionOutcome"]), QtCore.Qt.MatchFixedString
                )
                self.dataset_outcome_combobox_dict[xtal].setCurrentIndex(index)

            elif header[0].startswith("img"):
                if os.path.isfile(db_dict[header[1]]):
                    pixmap = QtGui.QPixmap(db_dict[header[1]])
                else:
                    pixmap = QtGui.QPixmap(
                        os.path.join(
                            os.getenv("XChemExplorer_DIR"),
                            "xce",
                            "image",
                            "IMAGE_NOT_AVAILABLE.png",
                        )
                    )
                image = QtGui.QLabel()
                image.resize(128, 80)
                image.setPixmap(pixmap.scaled(image.size(), QtCore.Qt.KeepAspectRatio))
                table.setCellWidget(row, column, image)

            elif header[0] == "Select":
                checkbox = QtGui.QCheckBox()
                checkbox.toggle()
                if table == self.deposition_table_apo:
                    if xtal not in self.deposition_table_apo_dict:
                        self.deposition_table_apo_dict[xtal] = checkbox
                if table == self.deposition_table_bound:
                    if xtal not in self.deposition_table_bound_dict:
                        self.deposition_table_bound_dict[xtal] = checkbox
                table.setCellWidget(row, column, checkbox)
                checkbox.setChecked(False)

            else:
                cell_text = QtGui.QTableWidgetItem()
                # in case data collection failed for whatever reason
                try:
                    cell_text.setText(str(db_dict[header[1]]))
                except KeyError:  # older pkl files may not have all the columns
                    cell_text.setText("n/a")
                cell_text.setTextAlignment(
                    QtCore.Qt.AlignCenter | QtCore.Qt.AlignCenter
                )
                table.setItem(row, column, cell_text)
                print(
                    (
                        "row: {0!s}  column: {1!s}  value: {2!s}  header: {3!s}".format(
                            row, column, cell_text, header[0]
                        )
                    )
                )
                print(("column_name {0!s}".format(column_name)))

    def populate_datasets_summary_table_NEW(self):
        self.status_bar.showMessage(
            "Building summary table for data processing results; be patient this may"
            " take a while"
        )

        # get information about all samples collected during the current visit
        visit, beamline = XChemMain.getVisitAndBeamline(self.beamline_directory)
        if self.read_agamemnon.isChecked():
            visit = []
            for v in glob.glob(
                os.path.join(
                    self.beamline_directory[: self.beamline_directory.rfind("-") + 1]
                    + "*"
                )
            ):
                visit.append(v[v.rfind("/") + 1 :])

        self.update_log.insert(
            "reading information about collected crystals from database..."
        )
        collectedXtalsDict = self.db.xtals_collected_during_visit_as_dict(visit)

        # instead of using dictionaries, query table of which crystals are in table
        samples_in_table = self.get_sample_list_from_table(self.datasets_summary_table)
        for xtal in sorted(collectedXtalsDict):
            if xtal not in samples_in_table:
                row = self.datasets_summary_table.rowCount()
                self.datasets_summary_table.insertRow(row)
            else:
                row = self.get_row_of_sample_in_table(self.datasets_summary_table, xtal)
            db_dict = collectedXtalsDict[xtal]
            self.update_row_in_table(
                xtal,
                row,
                db_dict,
                self.datasets_summary_table,
                self.datasets_summary_table_columns,
            )

        self.datasets_summary_table.resizeRowsToContents()
        self.datasets_summary_table.resizeColumnsToContents()

        self.status_bar.showMessage("updating Overview table")

        self.status_bar.showMessage("idle")

    def get_selected_row(self, table):
        indexes = table.selectionModel().selectedRows()
        for index in sorted(indexes):
            selected_row = index.row()
        return selected_row

    def show_results_from_all_pipelines(self):
        selected_row = self.get_selected_row(self.datasets_summary_table)
        xtal = self.datasets_summary_table.item(selected_row, 0).text()
        # get details of currently selected autoprocessing result
        selectedResultDict = self.db.get_db_dict_for_sample(xtal)

        dbList = self.db.all_autoprocessing_results_for_xtal_as_dict(xtal)

        self.make_data_collection_table()
        # needs to be created here, otherwise the cellClicked function
        dialog = QtGui.QDialog()
        # will reference it before it exists
        for db_dict in dbList:
            if (
                str(db_dict["DataProcessingSpaceGroup"]).lower() == "null"
                or str(db_dict["DataProcessingSpaceGroup"]).lower() == "none"
            ):
                continue
            row = self.data_collection_table.rowCount()
            self.data_collection_table.insertRow(row)
            self.update_row_in_table(
                xtal,
                row,
                db_dict,
                self.data_collection_table,
                self.data_collection_table_columns,
            )
            if (
                selectedResultDict["DataCollectionVisit"]
                == db_dict["DataCollectionVisit"]
                and selectedResultDict["DataCollectionRun"]
                == db_dict["DataCollectionRun"]
                and selectedResultDict["DataProcessingProgram"]
                == db_dict["DataProcessingProgram"]
                and selectedResultDict["DataProcessingScore"]
                == db_dict["DataProcessingScore"]
            ):
                self.current_row = row
                self.data_collection_table.selectRow(row)
        self.data_collection_table.cellClicked.connect(
            lambda: self.select_different_autoprocessing_result(dialog)
        )
        self.data_collection_table_popup(dialog)

    def make_data_collection_table(self):
        # this creates a new table widget every time
        # more elegant would be to delete or reset an existing widget...
        self.data_collection_table = QtGui.QTableWidget()
        self.data_collection_table.setVerticalScrollBarPolicy(
            QtCore.Qt.ScrollBarAlwaysOff
        )
        self.data_collection_table.setColumnCount(
            len(self.data_collection_table_columns)
        )
        font = QtGui.QFont()
        font.setPointSize(8)
        self.data_collection_table.setFont(font)
        self.data_collection_table.setHorizontalHeaderLabels(
            self.data_collection_table_columns
        )
        self.data_collection_table.horizontalHeader().setFont(font)
        self.data_collection_table.setSelectionBehavior(
            QtGui.QAbstractItemView.SelectRows
        )
        self.data_collection_table.setMinimumWidth(1000)
        self.data_collection_table.setMinimumHeight(500)

    def data_collection_table_popup(self, dialog):
        dialog_layout = QtGui.QGridLayout(dialog)
        dialog_layout.addWidget(self.data_collection_table)
        dialog.exec_()

    def select_different_autoprocessing_result(self, dialog):
        selected_row = self.get_selected_row(self.data_collection_table)
        if selected_row != self.current_row:
            xtal = self.data_collection_table.item(selected_row, 0).text()
            visit = self.data_collection_table.item(selected_row, 1).text()
            run = self.data_collection_table.item(selected_row, 2).text()
            autoproc = self.data_collection_table.item(selected_row, 3).text()
            score = self.data_collection_table.item(selected_row, 12).text()
            for q in range(13):
                try:
                    print(
                        (
                            "--> {0!s}: {1!s}".format(
                                q,
                                self.data_collection_table.item(selected_row, q).text(),
                            )
                        )
                    )
                except AttributeError:
                    print(("--> {0!s}: None".format(q)))
            dbDict = self.db.get_db_dict_for_visit_run_autoproc_score(
                xtal, visit, run, autoproc, score
            )
            dbDict["DataProcessingAutoAssigned"] = "False"
            self.update_log.insert(
                "%s: changing selected autoprocessing result to %s %s %s"
                % (xtal, visit, run, autoproc)
            )
            # xtal is QString -> str(xtal)
            XChemMain.linkAutoProcessingResult(
                str(xtal), dbDict, self.initial_model_directory, self.xce_logfile
            )
            self.update_log.insert("%s: updating row in Datasets table" % xtal)
            self.db.update_data_source(str(xtal), dbDict)
            self.update_log.insert(
                "%s: getting updated information from DB mainTable" % xtal
            )
            dbDict = self.db.get_db_dict_for_sample(xtal)
            row = self.get_row_of_sample_in_table(self.datasets_summary_table, xtal)
            self.update_row_in_table(
                xtal,
                row,
                dbDict,
                self.datasets_summary_table,
                self.datasets_summary_table_columns,
            )
        else:
            print("nothing to change")
        dialog.done(1)

    # < end
    ####################################################################################

    def user_update_selected_autoproc_datasets_summary_table(self):
        for key in self.data_collection_column_three_dict:
            if self.data_collection_column_three_dict[key][0] == self.sender():
                self.update_log.insert("here: " + self.sender())
                self.update_log.insert(
                    "herere" + str(self.data_collection_column_three_dict)
                )
                dbTmp = self.xtal_db_dict[key]
                stage = dbTmp["RefinementOutcome"].split()[0]
                print(("===>", key, stage))
                if int(stage) > 2:
                    msgBox = QtGui.QMessageBox()
                    msgBox.setText(
                        "*** WARNING ***\n"
                        "%s is currently %s\n"
                        "It will disappear from the Refinement table,\n"
                        "when you refresh it next time.\n"
                        "Do you want to continue?" % (key, dbTmp["RefinementOutcome"])
                    )
                    msgBox.addButton(QtGui.QPushButton("No"), QtGui.QMessageBox.YesRole)
                    msgBox.addButton(
                        QtGui.QPushButton("Yes"), QtGui.QMessageBox.RejectRole
                    )
                    reply = msgBox.exec_()
                    if reply == 0:
                        self.update_log.insert(
                            "will not change data processing selection"
                        )
                        # restore previous selection
                        for n, entry in enumerate(self.data_collection_dict[key]):
                            print(("==>", n))
                            if entry[0] == "logfile":
                                if entry[8]:
                                    print(("===> found:", n))
                                    self.data_collection_column_three_dict[key][
                                        0
                                    ].selectRow(n)
                        break

                indexes = self.sender().selectionModel().selectedRows()
                selected_processing_result = 1000000
                for index in sorted(indexes):
                    selected_processing_result = index.row()
                # the user changed the selection, i.e. no automated selection will
                # update it
                self.update_log.insert("user changed selection")
                self.data_collection_column_three_dict[key][1] = True
                # need to also update if not yet done
                user_already_changed_selection = False
                for n, entry in enumerate(self.data_collection_dict[key]):
                    if entry[0] == "user_changed_selection":
                        user_already_changed_selection = True
                    if entry[0] == "logfile":
                        db_dict = entry[6]
                        db_dict["DataProcessingAutoAssigned"] = "False"
                        if entry[7] == selected_processing_result:
                            db_dict_current = entry[6]
                            program = db_dict["DataProcessingProgram"]
                            visit = db_dict["DataCollectionVisit"]
                            run = db_dict["DataCollectionRun"]
                            self.update_log.insert(
                                "user changed data processing files for {0!s}"
                                " to visit={1!s}, run={2!s}, program={3!s}".format(
                                    key, visit, run, program
                                )
                            )
                            entry[8] = True
                        else:
                            entry[8] = False

                        entry[6] = db_dict
                        self.data_collection_dict[key][n] = entry
                if not user_already_changed_selection:
                    self.data_collection_dict[key].append(["user_changed_selection"])
                XChemMain.change_links_to_selected_data_collection_outcome(
                    key,
                    self.data_collection_dict,
                    self.data_collection_column_three_dict,
                    self.dataset_outcome_dict,
                    self.initial_model_directory,
                    os.path.join(self.database_directory, self.data_source_file),
                    self.xce_logfile,
                )

                # update 'Datasets' table
                column_name = XChemDB.data_source(
                    os.path.join(self.database_directory, self.data_source_file)
                ).translate_xce_column_list_to_sqlite(
                    self.datasets_summary_table_columns
                )
                rows_in_table = self.datasets_summary_table.rowCount()
                for row in range(rows_in_table):
                    if self.datasets_summary_table.item(row, 0).text() == key:
                        for column, header in enumerate(column_name):
                            if header[0] == "Sample ID":
                                continue
                            elif header[0] == "DataCollection\nOutcome":
                                continue
                            elif header[0].startswith("img"):
                                continue
                            elif header[0].startswith("Show"):
                                continue
                            else:
                                cell_text = QtGui.QTableWidgetItem()
                                try:
                                    cell_text.setText(str(db_dict_current[header[1]]))
                                    cell_text.setTextAlignment(
                                        QtCore.Qt.AlignCenter | QtCore.Qt.AlignCenter
                                    )
                                    self.datasets_summary_table.setItem(
                                        row, column, cell_text
                                    )
                                except KeyError:
                                    pass

    def populate_and_update_datasource_table(self):
        self.overview_datasource_table.setColumnCount(
            len(self.overview_datasource_table_columns)
        )

        # first get a list of all the samples that are already in the table and which
        # will be updated
        samples_in_table = []
        current_row = self.overview_datasource_table.rowCount()
        for row in range(current_row):
            sampleID = str(
                self.overview_datasource_table.item(row, 0).text()
            )  # this must be the case
            samples_in_table.append(sampleID)

        columns_to_show = self.get_columns_to_show(
            self.overview_datasource_table_columns
        )
        sample_id_column = self.get_columns_to_show(["Sample ID"])

        for row in self.data:
            if (
                str(row[sample_id_column[0]]).lower() == "none"
                or str(row[sample_id_column[0]]).replace(" ", "") == ""
            ):
                # do not show rows where sampleID is null
                continue
            else:
                if not str(row[sample_id_column[0]]) in samples_in_table:
                    # insert row, this is a new sample
                    x = self.overview_datasource_table.rowCount()
                    self.overview_datasource_table.insertRow(x)
                else:
                    # find row of this sample in data_source_table
                    for present_rows in range(
                        self.overview_datasource_table.rowCount()
                    ):
                        if str(row[sample_id_column[0]]) == str(
                            self.overview_datasource_table.item(present_rows, 0).text()
                        ):
                            x = present_rows
                            break
            for y, item in enumerate(columns_to_show):
                cell_text = QtGui.QTableWidgetItem()
                if row[item] is None:
                    cell_text.setText("")
                else:
                    cell_text.setText(str(row[item]))
                # assumption is that column 0 is always sampleID
                if self.overview_datasource_table_columns[y] == "Sample ID":
                    # and this field cannot be changed
                    cell_text.setFlags(QtCore.Qt.ItemIsEnabled)
                cell_text.setTextAlignment(
                    QtCore.Qt.AlignCenter | QtCore.Qt.AlignCenter
                )
                self.overview_datasource_table.setItem(x, y, cell_text)
        self.overview_datasource_table.setHorizontalHeaderLabels(
            self.overview_datasource_table_columns
        )

    def kill_other_pandda_options(self):
        for i in range(0, self.pandda_analyse_data_table.rowCount()):
            checkbox1 = self.pandda_analyse_data_table.cellWidget(i, 7)
            checkbox2 = self.pandda_analyse_data_table.cellWidget(i, 8)
            checkbox3 = self.pandda_analyse_data_table.cellWidget(i, 9)
            if checkbox1.isChecked():
                checkbox2.setChecked(False)
                checkbox3.setChecked(False)
            if checkbox1.isChecked() and checkbox2.isChecked() or checkbox3.isChecked():
                checkbox1.setChecked(False)
            if checkbox2.isChecked() or checkbox3.isChecked():
                checkbox1.setChecked(False)

    def populate_pandda_analyse_input_table(self):
        column_name = self.db.translate_xce_column_list_to_sqlite(
            self.pandda_table_columns
        )
        print(column_name)
        for xtal in sorted(self.xtal_db_dict):
            new_xtal = False
            db_dict = self.xtal_db_dict[xtal]
            if os.path.isfile(db_dict["DimplePathToPDB"]):
                row = self.pandda_analyse_data_table.rowCount()
                if xtal not in self.pandda_analyse_input_table_dict:
                    self.pandda_analyse_data_table.insertRow(row)
                    current_row = row
                    new_xtal = True
                else:
                    for table_row in range(row):
                        if (
                            self.pandda_analyse_data_table.item(table_row, 0).text()
                            == xtal
                        ):
                            current_row = table_row
                            break
                for column, header in enumerate(column_name):
                    if header[0] == "Exclude":
                        deselect_button = QtGui.QCheckBox()
                        deselect_button.stateChanged.connect(
                            self.kill_other_pandda_options
                        )
                        self.pandda_analyse_data_table.setCellWidget(
                            current_row, column, deselect_button
                        )

                    elif header[0] == "Ignore":
                        deselect_button = QtGui.QCheckBox()
                        deselect_button.stateChanged.connect(
                            self.kill_other_pandda_options
                        )
                        self.pandda_analyse_data_table.setCellWidget(
                            current_row, column, deselect_button
                        )

                    elif header[0] == "Export":
                        deselect_button = QtGui.QCheckBox()
                        deselect_button.stateChanged.connect(
                            self.kill_other_pandda_options
                        )
                        self.pandda_analyse_data_table.setCellWidget(
                            current_row, column, deselect_button
                        )

                    elif header[0] == "Sample ID":
                        cell_text = QtGui.QTableWidgetItem()
                        cell_text.setText(str(xtal))
                        cell_text.setTextAlignment(
                            QtCore.Qt.AlignCenter | QtCore.Qt.AlignCenter
                        )
                        self.pandda_analyse_data_table.setItem(
                            current_row, column, cell_text
                        )
                    else:
                        cell_text = QtGui.QTableWidgetItem()
                        cell_text.setText(str(db_dict[header[1]]))
                        if header[0] == "PanDDA\nStatus":
                            if str(db_dict[header[1]]) == "running":
                                cell_text.setBackground(QtGui.QColor(100, 230, 150))
                            elif str(db_dict[header[1]]) == "pending":
                                cell_text.setBackground(QtGui.QColor(20, 100, 230))
                            elif str(db_dict[header[1]]) == "started":
                                cell_text.setBackground(QtGui.QColor(230, 240, 110))
                            elif str(db_dict[header[1]]) == "finished":
                                cell_text.setBackground(QtGui.QColor(255, 255, 255))
                            elif "problem" in str(db_dict[header[1]]):
                                cell_text.setBackground(QtGui.QColor(255, 0, 0))
                        cell_text.setTextAlignment(
                            QtCore.Qt.AlignCenter | QtCore.Qt.AlignCenter
                        )
                        self.pandda_analyse_data_table.setItem(
                            current_row, column, cell_text
                        )
            if new_xtal:
                self.pandda_analyse_input_table_dict[xtal] = []

    def select_sample_for_pandda(self, option):
        indexes = self.pandda_analyse_data_table.selectionModel().selectedRows()
        if option == "deselect":
            for index in sorted(indexes):
                self.pandda_analyse_data_table.cellWidget(index.row(), 6).setChecked(
                    False
                )
                self.pandda_analyse_data_table.cellWidget(index.row(), 7).setChecked(
                    False
                )
                self.pandda_analyse_data_table.cellWidget(index.row(), 8).setChecked(
                    False
                )
        else:
            for index in sorted(indexes):
                self.pandda_analyse_data_table.cellWidget(index.row(), 6).setChecked(
                    False
                )
                self.pandda_analyse_data_table.cellWidget(index.row(), 7).setChecked(
                    False
                )
                self.pandda_analyse_data_table.cellWidget(index.row(), 8).setChecked(
                    False
                )
                if option == "ignore":
                    checkbox = self.pandda_analyse_data_table.cellWidget(index.row(), 6)
                if option == "char":
                    checkbox = self.pandda_analyse_data_table.cellWidget(index.row(), 7)
                if option == "zmap":
                    checkbox = self.pandda_analyse_data_table.cellWidget(index.row(), 8)

                checkbox.setChecked(True)
            self.kill_other_pandda_options()

    def populate_and_update_refinement_table(self):
        column_name = self.db.translate_xce_column_list_to_sqlite(
            self.refinement_table_columns
        )
        for xtal in sorted(self.xtal_db_dict):
            new_xtal = False
            db_dict = self.xtal_db_dict[xtal]
            try:
                stage = int(str(db_dict["RefinementOutcome"]).split()[0])
                refinementStage = db_dict["RefinementOutcome"]
            except ValueError:
                stage = 0
            except IndexError:
                stage = 0

            if stage >= 3 and stage < 7:
                row = self.refinement_table.rowCount()
                if xtal not in self.refinement_table_dict:
                    self.refinement_table.insertRow(row)
                    current_row = row
                    new_xtal = True
                else:
                    for table_row in range(row):
                        if self.refinement_table.item(table_row, 0).text() == xtal:
                            current_row = table_row
                            break
                for column, header in enumerate(column_name):
                    if header[0] == "Sample ID":
                        cell_text = QtGui.QTableWidgetItem()
                        cell_text.setText(str(xtal))
                        cell_text.setTextAlignment(
                            QtCore.Qt.AlignCenter | QtCore.Qt.AlignCenter
                        )
                        self.refinement_table.setItem(current_row, column, cell_text)

                    elif header[0] == "Refinement\nOutcome":
                        if new_xtal:
                            refinement_outcome_combobox = QtGui.QComboBox()
                            self.populate_refinement_outcome_combobox(
                                refinement_outcome_combobox
                            )
                            self.refinement_table.setCellWidget(
                                current_row, column, refinement_outcome_combobox
                            )
                        else:
                            refinement_outcome_combobox = self.refinement_table_dict[
                                xtal
                            ]
                        index = refinement_outcome_combobox.findText(
                            refinementStage, QtCore.Qt.MatchFixedString
                        )
                        refinement_outcome_combobox.setCurrentIndex(index)
                        refinement_outcome_combobox.currentIndexChanged.connect(
                            self.refinement_outcome_combobox_changed
                        )

                    elif header[0] == "buster-reports":
                        buster_report = db_dict["RefinementBusterReportHTML"]
                        ref_name = buster_report.split("/")[
                            len(buster_report.split("/")) - 2
                        ]
                        buster_report_link = QtGui.QLabel(
                            '<a href="{0!s}">{1!s}</a>'.format(buster_report, ref_name)
                        )
                        buster_report_link.setOpenExternalLinks(True)
                        self.refinement_table.setCellWidget(
                            current_row, column, buster_report_link
                        )
                    else:
                        cell_text = QtGui.QTableWidgetItem()
                        cell_text.setText(str(db_dict[header[1]]))
                        if header[0] == "Refinement\nStatus":
                            if str(db_dict[header[1]]) == "running":
                                cell_text.setBackground(QtGui.QColor(100, 230, 150))
                            elif str(db_dict[header[1]]) == "pending":
                                cell_text.setBackground(QtGui.QColor(20, 100, 230))
                            elif str(db_dict[header[1]]) == "started":
                                cell_text.setBackground(QtGui.QColor(230, 240, 110))
                            elif str(db_dict[header[1]]) == "finished":
                                cell_text.setBackground(QtGui.QColor(255, 255, 255))
                            elif "problem" in str(db_dict[header[1]]):
                                cell_text.setBackground(QtGui.QColor(255, 0, 0))
                        cell_text.setTextAlignment(
                            QtCore.Qt.AlignCenter | QtCore.Qt.AlignCenter
                        )
                        self.refinement_table.setItem(current_row, column, cell_text)
            if new_xtal:
                self.refinement_table_dict[xtal] = refinement_outcome_combobox

        self.refinement_table.resizeColumnsToContents()
        self.refinement_table.resizeRowsToContents()

    def get_columns_to_show(self, column_list):
        # maybe I coded some garbage before, but I need to find out which column name
        # in the data source corresponds to the actually displayed column name in the
        # table reason being that the unique column ID for DB may not be nice to look at
        columns_to_show = []
        for column in column_list:
            # first find out what the column name in the header is:
            column_name = ""
            for name in self.all_columns_in_data_source:
                if column == name[1]:
                    column_name = name[0]
            for n, all_column in enumerate(self.header):
                if column_name == all_column:
                    columns_to_show.append(n)
                    break
        return columns_to_show

    def quit_xce(self):
        # save pkl file
        if self.data_collection_dict != {}:
            if os.path.isfile(self.datasets_summary_file):
                self.update_log.insert("saving results to PKL file")
                pickle.dump(
                    self.data_collection_dict, open(self.datasets_summary_file, "wb")
                )
        self.update_log.insert("quitting XCE... bye,bye!")
        QtGui.qApp.quit()
