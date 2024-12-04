import os
import subprocess
import sys
from PyQt4 import QtCore
from xce.lib import XChemDB, XChemLog, XChemMain, XChemUtils


class setup:
    def __init__(self):
        pass

    def openFile(self, file):
        if sys.platform == "linux2":
            subprocess.call(["xdg-open", file])
        else:
            os.startfile(file)

    def set_xce_logfile(self, xce_object):
        XChemLog.startLog(xce_object.xce_logfile).create_logfile(xce_object.xce_version)
        xce_object.update_log = XChemLog.updateLog(xce_object.xce_logfile)

    def settings(self, xce_object):
        # set XCE version
        xce_object.xce_version = "v2.0.1"

        # general settings
        xce_object.allowed_unitcell_difference_percent = 12
        xce_object.acceptable_low_resolution_limit_for_data = 3.5
        xce_object.filename_root = "${samplename}"
        xce_object.data_source_set = False
        xce_object.max_queue_jobs = 100

        # directory settings

        # set current directory and direct log to it
        xce_object.current_directory = os.getcwd()
        xce_object.xce_logfile = os.path.join(xce_object.current_directory, "xce.log")

        # if in the correct place, set the various directories
        if xce_object.current_directory.startswith("/dls/labxchem"):
            #        if 'labxchem' in xce_object.current_directory:
            if (
                len(xce_object.current_directory.split("/")) >= 9
                and xce_object.current_directory.split("/")[6] == "processing"
                and xce_object.current_directory.split("/")[8] == "processing"
            ):
                xce_object.labxchem_directory = "/" + os.path.join(
                    *xce_object.current_directory.split("/")[1:8]
                )  # need splat operator: *
                xce_object.labxchem_directory_current = "/" + os.path.join(
                    *xce_object.current_directory.split("/")[1:9]
                )
                # labxchem_directory_current is where they actually have write
                # permission
            else:
                xce_object.labxchem_directory = "/" + os.path.join(
                    *xce_object.current_directory.split("/")[1:6]
                )  # need splat operator: *
                xce_object.labxchem_directory_current = "/" + os.path.join(
                    *xce_object.current_directory.split("/")[1:7]
                )  # need splat operator: *
            xce_object.beamline_directory = os.path.join(
                xce_object.labxchem_directory, "processing", "beamline"
            )
            if os.path.isdir(
                os.path.join(
                    xce_object.labxchem_directory,
                    "processing",
                    "analysis",
                    "model_building",
                )
            ):
                xce_object.initial_model_directory = os.path.join(
                    xce_object.labxchem_directory,
                    "processing",
                    "analysis",
                    "model_building",
                )
            else:
                xce_object.initial_model_directory = os.path.join(
                    xce_object.labxchem_directory,
                    "processing",
                    "analysis",
                    "initial_model",
                )
            xce_object.reference_directory = os.path.join(
                xce_object.labxchem_directory, "processing", "reference"
            )
            xce_object.database_directory = os.path.join(
                xce_object.labxchem_directory, "processing", "database"
            )
            xce_object.panddas_directory = os.path.join(
                xce_object.labxchem_directory, "processing", "analysis", "panddas"
            )

            xce_object.datasets_summary_file = None
            xce_object.data_source_file = ""
            xce_object.html_export_directory = os.path.join(
                xce_object.labxchem_directory, "processing", "html"
            )
            xce_object.group_deposit_directory = os.path.join(
                xce_object.labxchem_directory, "processing", "group_deposition"
            )
            if os.path.isfile(
                os.path.join(
                    xce_object.labxchem_directory,
                    "processing",
                    "database",
                    "soakDBDataFile.sqlite",
                )
            ):
                xce_object.data_source_file = "soakDBDataFile.sqlite"
                xce_object.database_directory = os.path.join(
                    xce_object.labxchem_directory, "processing", "database"
                )
                xce_object.data_source_set = True
                xce_object.db = XChemDB.data_source(
                    os.path.join(
                        xce_object.database_directory, xce_object.data_source_file
                    )
                )
                xce_object.db.create_missing_columns()

            xce_object.ccp4_scratch_directory = os.path.join(
                xce_object.labxchem_directory, "processing", "tmp"
            )

            directory_list = [
                xce_object.beamline_directory,
                os.path.join(xce_object.labxchem_directory, "processing", "analysis"),
                xce_object.initial_model_directory,
                xce_object.panddas_directory,
                xce_object.reference_directory,
                xce_object.database_directory,
                xce_object.ccp4_scratch_directory,
                xce_object.html_export_directory,
                xce_object.group_deposit_directory,
            ]

            for directory in directory_list:
                if not os.path.isdir(directory):
                    os.mkdir(directory)

        # otherwise, use the current working directory
        else:
            xce_object.labxchem_directory_current = xce_object.current_directory
            xce_object.beamline_directory = xce_object.current_directory
            xce_object.initial_model_directory = xce_object.current_directory
            xce_object.reference_directory = xce_object.current_directory
            xce_object.database_directory = xce_object.current_directory
            xce_object.data_source_file = ""
            xce_object.ccp4_scratch_directory = os.getenv("CCP4_SCR")
            xce_object.panddas_directory = xce_object.current_directory
            xce_object.datasets_summary_file = ""
            xce_object.group_deposit_directory = xce_object.current_directory

        # deposition

        xce_object.deposit_dict = {}

        # internal lists and dictionaries

        xce_object.visit_list = []
        xce_object.target = ""
        xce_object.dataset_outcome_combobox_dict = {}
        xce_object.data_collection_dict = {}
        xce_object.xtal_db_dict = {}
        xce_object.pandda_analyse_input_table_dict = {}
        # contains toggle button if dimple should be run
        xce_object.initial_model_dimple_dict = {}
        xce_object.reference_file_list = []
        xce_object.all_columns_in_data_source = XChemDB.data_source(
            os.path.join(xce_object.database_directory, xce_object.data_source_file)
        ).return_column_list()

        xce_object.dataset_outcome_dict = {}  # contains the dataset outcome buttons
        xce_object.data_collection_table_dict = {}  # contains the dataset table
        xce_object.data_collection_image_dict = {}
        xce_object.data_collection_column_three_dict = {}
        xce_object.refinement_table_dict = {}
        xce_object.timer_to_check_for_new_data_collection = QtCore.QTimer()

        xce_object.agamemnon = False
        (
            xce_object.target_list,
            xce_object.visit_list,
        ) = XChemMain.get_target_and_visit_list(xce_object.beamline_directory, False)

        # internal switches and flags

        xce_object.explorer_active = 0
        xce_object.coot_running = 0
        xce_object.gdaLogInstructions = [0, False]

        xce_object.dataset_outcome = [
            "success",
            "Failed - centring failed",
            "Failed - no diffraction",
            "Failed - processing",
            "Failed - loop empty",
            "Failed - loop broken",
            "Failed - low resolution",
            "Failed - no X-rays",
            "Failed - unknown",
        ]

        xce_object.refinement_stage = [
            "0 - All Datasets",
            "1 - Analysis Pending",
            "2 - PANDDA model",
            "3 - In Refinement",
            "4 - CompChem ready",
            "5 - Deposition ready",
            "6 - Deposited",
            "7 - Analysed & Rejected",
        ]

        self.set_xce_logfile(xce_object)

        # external software packages
        xce_object.update_log = XChemLog.updateLog(xce_object.xce_logfile)
        xce_object.update_log.insert("new session started")
        xce_object.diffraction_data_directory = xce_object.current_directory
        xce_object.html_export_directory = os.getcwd()
        xce_object.external_software = XChemUtils.external_software(
            xce_object.xce_logfile
        ).check()

        xce_object.second_cif_file = None

        restraints_program_candidates = ["acedrg", "phenix.elbow", "grade"]

        xce_object.restraints_program = ""
        for restraints_program_candidate in restraints_program_candidates:
            if restraints_program_candidate in xce_object.external_software.keys():
                xce_object.restraints_program = restraints_program_candidate
                xce_object.update_log.insert(
                    "will use {0!s} for"
                    " generation of ligand coordinates and restraints".format(
                        restraints_program_candidate
                    )
                )
                break
        if xce_object.restraints_program is None:
            xce_object.update_log.warning(
                "No program for generation of ligand coordinates and restraints"
                " available!"
            )

    def preferences(self, xce_object):
        # preferences

        xce_object.preferences_data_to_copy = [
            ["aimless logiles and merged mtz only", "mtz_log_only"],
        ]

        xce_object.preferences_selection_mechanism = [
            "IsigI*Comp*UniqueRefl",
            "highest_resolution",
            "lowest_Rfree",
            "dials - only",
            "xia2 3dii - only",
            "autoProc - only",
            "autoProc_staraniso - only",
        ]

        xce_object.allowed_unitcell_difference_percent = 12
        xce_object.acceptable_low_resolution_limit_for_data = 3.5
        xce_object.filename_root = "${samplename}"
        xce_object.max_queue_jobs = 100
        xce_object.dimple_twin_mode = False

        xce_object.preferences = {
            "processed_data_to_copy": "mtz_log_only",
            "dataset_selection_mechanism": "IsigI*Comp*UniqueRefl",
            "allowed_unitcell_difference_percent": 12,
            "acceptable_low_resolution_limit_for_data": 3.5,
            "acceptable_low_resolution_Rmerge": 0.1,
            "filename_root": "${samplename}",
            "max_queue_jobs": 100,
            "dimple_twin_mode": False,
            "initial_refinement_pipeline": "dimple",
        }

        # settings

        xce_object.settings = {
            "current_directory": xce_object.current_directory,
            "beamline_directory": xce_object.beamline_directory,
            "datasets_summary": xce_object.datasets_summary_file,
            "initial_model_directory": xce_object.initial_model_directory,
            "panddas_directory": xce_object.panddas_directory,
            "reference_directory": xce_object.reference_directory,
            "database_directory": xce_object.database_directory,
            "data_source": os.path.join(
                xce_object.database_directory, xce_object.data_source_file
            ),
            "ccp4_scratch": xce_object.ccp4_scratch_directory,
            "unitcell_difference": xce_object.allowed_unitcell_difference_percent,
            "too_low_resolution_data": (
                xce_object.acceptable_low_resolution_limit_for_data
            ),
            "filename_root": xce_object.filename_root,
            "preferences": xce_object.preferences,
            "xce_logfile": xce_object.xce_logfile,
            "max_queue_jobs": xce_object.max_queue_jobs,
            "diffraction_data_directory": xce_object.diffraction_data_directory,
            "html_export_directory": xce_object.html_export_directory,
            "group_deposit_directory": xce_object.group_deposit_directory,
            "dimple_twin_mode": False,
            "agamemnon": False,
        }

    def tables(self, xce_object):
        # Table column settings

        # functions that use tables.overview_datasource_table_columns:
        #
        # select_datasource_columns_to_display()
        # - dropdown in datasource top menu (select columns to show)
        # populate_data_source_table()
        # - appears to be completely unused, so commented out
        # populate_and_update_datasource_table()
        # - used within select_datasource_columns_to_display and
        # update_all_tables()

        xce_object.overview_datasource_table_columns = [
            "Sample ID",
            "Compound ID",
            "Smiles",
            "Visit",
            "Resolution\n[Mn<I/sig(I)> = 1.5]",
            "Refinement\nRfree",
            "Data Collection\nDate",
            "Puck",
            "PuckPosition",
            "Ligand\nConfidence",
        ]

        # functions that use tables.datasets_summary_table_columns:
        #
        # populate_datasets_summary_table()
        # - appears in create_widgets_for_autoprocessing_results_only()
        # user_update_selected_autoproc_datasets_summary_table()
        # - appears in create_widgets_for_autoprocessing_results_only()

        xce_object.datasets_summary_table_columns = [
            "Sample ID",
            "Resolution\nHigh",
            "DataProcessing\nSpaceGroup",
            "DataProcessing\nRfree",
            "SoakDB\nBarcode",
            "GDA\nBarcode",
            "Rmerge\nLow",
            "auto-assigned",
            "DataCollection\nOutcome",
            "img1",
            "img2",
            "img3",
            "img4",
        ]

        # functions that use tables.data_collection_table_columns:
        #
        # show_results_from_all_pipelines()
        # - appears in populate_datasets_summary_table()

        xce_object.data_collection_table_columns = [
            "Sample ID",
            "Visit",
            "Run",
            "Program",
            "Resolution\nOverall",
            "Resolution\nHigh",
            "DataProcessing\nSpaceGroup",
            "Mn<I/sig(I)>\nHigh",
            "Rmerge\nLow",
            "Completeness\nOverall",
            "DataProcessing\nUnitCell",
            "DataProcessing\nRfree",
            "DataProcessing\nScore",
        ]

        # functions that use tables.maps_table_columns:
        #
        # - appears in create_maps_table()

        xce_object.maps_table_columns = [
            "Sample ID",
            "Select",
            "Compound ID",
            "Smiles",
            "Resolution\nHigh",
            "Dimple\nRcryst",
            "Dimple\nRfree",
            "DataProcessing\nSpaceGroup",
            "Reference\nSpaceGroup",
            "Difference\nUC Volume (%)",
            "Reference File",
            "DataProcessing\nUnitCell",
            "Dimple\nStatus",
            "Compound\nStatus",
            "LastUpdated",
        ]

        # functions that use tables.pandda_table_columns:
        #
        # populate_pandda_analyse_input_table()
        # - appears in update_all_tables()

        xce_object.pandda_table_columns = [
            "Sample ID",
            "Export\nSelected",
            "Refinement\nSpace Group",
            "Resolution\n[Mn<I/sig(I)> = 1.5]",
            "Dimple\nRcryst",
            "Dimple\nRfree",
            "Crystal Form\nName",
            "Ignore\ncompletely",
            "Exclude from\n characterisation\n(binds)",
            "Exclude from\n z-map analysis\n(does not bind)",
        ]

        # functions that use tables.refinement_table_columns:
        #
        # populate_and_update_refinement_table()
        # - appears in update_all_tables

        xce_object.refinement_table_columns = [
            "Sample ID",
            "Compound ID",
            "Refinement\nSpace Group",
            "Refinement\nResolution",
            "Refinement\nRcryst",
            "Refinement\nRfree",
            "Refinement\nOutcome",
            "buster-reports",
            "Ligand CC",
            "Refinement\nStatus",
        ]

    def top_menu_dict(self, xce_object):
        xce_object.menu_dict = {
            "A: file": [
                "&File",
                [
                    ["Open Config File", "Ctrl+O", xce_object.open_config_file],
                    ["Save Config File", "Ctrl+S", xce_object.save_config_file],
                    ["Quit", "Ctrl+Q", xce_object.quit_xce],
                ],
            ],
            "B: datasource": [
                "&Datasource",
                [
                    [
                        "Reload Samples From Datasource",
                        "",
                        xce_object.datasource_menu_reload_samples,
                    ],
                    [
                        "Save Samples to Datasource",
                        "",
                        xce_object.datasource_menu_save_samples,
                    ],
                    [
                        "Import CSV file into Datasource",
                        "",
                        xce_object.datasource_menu_import_csv_file,
                    ],
                    [
                        "Export CSV file from Datasource",
                        "",
                        xce_object.datasource_menu_export_csv_file,
                    ],
                    [
                        "Update datasource from file system",
                        "",
                        xce_object.datasource_menu_update_datasource,
                    ],
                    [
                        "Select columns to show",
                        "",
                        xce_object.select_datasource_columns_to_display,
                    ],
                    [
                        "Create New Datasource (SQLite)",
                        "",
                        xce_object.create_new_data_source,
                    ],
                    ["Export CSV for wonka", "", xce_object.export_data_for_WONKA],
                ],
            ],
            "C: preferences": [
                "&Preferences",
                [["Edit preferences", "", xce_object.show_preferences]],
            ],
            "D: deposition": [
                "&Deposition",
                [
                    ["Edit information", "", xce_object.deposition_data],
                    ["Export to HTML", "", xce_object.export_to_html],
                    [
                        "Export to HTML - CompChem",
                        "",
                        xce_object.export_to_html_CompChem,
                    ],
                    [
                        "Export to HTML - deposition ready",
                        "",
                        xce_object.export_to_html_deposition_ready,
                    ],
                    ["Update DB with PDB codes", "", xce_object.enter_pdb_codes],
                    ["Check SMILES", "", xce_object.check_smiles_in_db_and_pdb],
                ],
            ],
            "F: help": [
                "&Help",
                [
                    [
                        "Open XCE manual",
                        "",
                        lambda: setup().openFile(
                            "/dls/science/groups/i04-1/software/"
                            "XCE_manual_2018-11-09.pdf"
                        ),
                    ],
                    [
                        "Open XCE tutorial",
                        "",
                        lambda: setup().openFile(
                            "/dls/science/groups/i04-1/software/docs/XChemExplorer.pdf"
                        ),
                    ],
                    [
                        "Troubleshooting",
                        "",
                        lambda: setup().openFile(
                            "/dls/science/groups/i04-1/software/xce_troubleshooting.pdf"
                        ),
                    ],
                ],
            ],
            "G: labels": [
                "&Labels",
                [["Edit label information", "", xce_object.add_label_information]],
            ],
        }

    def bottom_box_buttons(self, xce_object):
        self.dropdown_items(xce_object)

        xce_object.datasource_button_dict = {
            "datasource_button": [
                r"Update Tables\nFrom Datasource",
                [
                    [
                        "XChemToolTips.update_from_datasource_button_tip()",
                        # tooltip
                        "QPushButton { padding: 1px; margin: 1px; "
                        "background: rgb(197,197,197) }",
                        # stylesheet
                        "xce_object.headlineLabelfont",  # font
                        "xce_object.datasource_menu_reload_samples",
                    ]
                    # action
                ],
            ]
        }

        xce_object.dataset_task_run_button_dict = {
            "dataset_run_button": [
                r"Run",
                [
                    [
                        "XChemToolTips.dataset_task_run_button_tip()",  # tooltip
                        # stylesheet
                        "QPushButton { padding: 1px; margin: 1px }",
                        "",  # font
                        "xce_object.button_clicked",
                    ]  # action
                ],
            ]
        }

        xce_object.dataset_task_status_button_dict = {
            "dataset_status_button": [
                r"Status",
                [
                    [
                        "XChemToolTips.dataset_task_status_button_tip()",
                        # tooltip
                        "QPushButton { padding: 1px; margin: 1px }",
                        # stylesheet
                        "",  # font
                        "xce_object.button_clicked",
                    ]  # action
                ],
            ]
        }

        xce_object.map_cif_file_task_run_button_dict = {
            "dataset_run_button": [
                r"Run",
                [
                    [
                        "XChemToolTips.map_cif_file_task_run_button_tip()",
                        # tooltip
                        "QPushButton { padding: 1px; margin: 1px }",
                        # stylesheet
                        "",  # font
                        "xce_object.button_clicked",
                    ]  # action
                ],
            ]
        }

        xce_object.map_cif_file_task_status_button_dict = {
            "dataset_status_button": [
                r"Status",
                [
                    [
                        "XChemToolTips.map_cif_file_task_status_button_tip()",
                        # tooltip
                        "QPushButton { padding: 1px; margin: 1px }",
                        # stylesheet
                        "",  # font
                        "xce_object.button_clicked",
                    ]  # action
                ],
            ]
        }

        xce_object.panddas_file_task_run_button_dict = {
            "dataset_run_button": [
                r"Run",
                [
                    [
                        "XChemToolTips.panddas_file_task_run_button_tip()",
                        # tooltip
                        "QPushButton { padding: 1px; margin: 1px }",
                        # stylesheet
                        "",  # font
                        "xce_object.button_clicked",
                    ]  # action
                ],
            ]
        }

        xce_object.panddas_file_task_status_button_dict = {
            "dataset_status_button": [
                r"Status",
                [
                    [
                        "XChemToolTips.panddas_file_task_status_button_tip()",
                        # tooltip
                        "QPushButton { padding: 1px; margin: 1px }",
                        # stylesheet
                        "",  # font
                        "xce_object.button_clicked",
                    ]  # action
                ],
            ]
        }

        xce_object.refine_file_task_run_button_dict = {
            "dataset_run_button": [
                r"Run",
                [
                    [
                        "XChemToolTips.refine_file_task_run_button_tip()",
                        # tooltip
                        "QPushButton { padding: 1px; margin: 1px }",
                        # stylesheet
                        "",  # font
                        "xce_object.button_clicked",
                    ]  # action
                ],
            ]
        }

        xce_object.refine_file_task_status_button_dict = {
            "dataset_status_button": [
                r"Status",
                [
                    [
                        "XChemToolTips.refine_file_task_status_button_tip()",
                        # tooltip
                        "QPushButton { padding: 1px; margin: 1px }",
                        # stylesheet
                        "",  # font
                        "xce_object.button_clicked",
                    ]  # action
                ],
            ]
        }

    def dropdown_items(self, xce_object):
        xce_object.dataset_tasks = [
            "Get New Results from Autoprocessing",
            "Rescore Datasets",
            "Run xia2 on selected datasets",
            "Run xia2 on selected datasets - overwrite",
        ]

        xce_object.map_cif_file_tasks = [
            "Run DIMPLE on selected MTZ files",
            "Remove selected DIMPLE files",
            "Set DIMPLE output",
            "------------------------------------------",
            "Create CIF/PDB/PNG file of SELECTED compounds",
            "Merge ligand CIF file with selected compounds",
            "Restore original CIF file of selected compounds",
            "Fit ligands into maps after initial refinement",
            "------------------------------------------",
            "Run PIPEDREAM on selected MTZ files",
            "Remove selected PIPEDREAM files",
            "Set PIPEDREAM output",
            "------------------------------------------",
            "Run PHENIX.LIGAND_PIPELINE on selected MTZ files",
            "Remove selected PHENIX.LIGAND_PIPELINE files",
            "Set PHENIX.LIGAND_PIPELINE output",
        ]

        xce_object.panddas_file_tasks = [
            "pre-run for ground state model",
            "Build ground state model",
            "------------------------------------------",
            "pandda.analyse",
            "pandda.analyse (PanDDA2)",
            "pandda.inspect",
            "run pandda.inspect at home",
            "------------------------------------------",
            "Export NEW PANDDA models",
            "Export ALL PANDDA models",
            "Export SELECTED PANDDA models",
            "------------------------------------------",
            "refine ALL bound-state models with BUSTER",
            "refine NEW bound-state models with BUSTER",
            "------------------------------------------",
            "refine ALL bound-state models with BUSTER (no sanity check)",
            "refine NEW bound-state models with BUSTER (no sanity check)",
            "------------------------------------------",
            "Show HTML summary",
            "cluster datasets",
            "Event Map -> SF",
            "apo -> mmcif",
            "check modelled ligands",
        ]

        xce_object.refine_file_tasks = [
            "Open COOT - REFMAC refinement -",
            "Open COOT - BUSTER refinement -",
            "Open COOT - dimple_twin -",
            "",
        ]
