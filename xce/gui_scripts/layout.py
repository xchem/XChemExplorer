import os

from PyQt4 import QtGui

from xce.gui_scripts.datasets_tab import DatasetsTab
from xce.gui_scripts.deposition_tab import DepositionTab
from xce.gui_scripts.maps_tab import MapsTab
from xce.gui_scripts.overview_tab import OverviewTab
from xce.gui_scripts.pandda_tab import PanddaTab
from xce.gui_scripts.refinement_tab import RefinementTab
from xce.gui_scripts.settings_preferences import setup
from xce.gui_scripts.settings_tab import SettingsTab
from xce.gui_scripts import layout_functions


class LayoutObjects:
    def initialise_menu_bar(self, xce_object):
        ################################################################################
        #                                                                              #
        #                             MENU BAR - TOP OF GUI                            #
        #                                                                              #
        ################################################################################

        # initiate menu widget
        menu_bar = QtGui.QMenuBar()

        # import menu bar dictionary
        setup().top_menu_dict(xce_object)

        # create menu from menu dictionary
        menu_bar = layout_functions.setup_menubar(
            xce_object, menu_bar, xce_object.menu_dict
        )

        # END OF MENU BAR - CODE BELOW: stuff removed from apo structure stuff that
        # appears might have a funky consequence - work out later.
        return menu_bar

    # function containing setup for bottom boxes
    def initialise_bottom_boxes(self, xce_object):
        icons_directory = os.path.join((os.getenv("XChemExplorer_DIR")), "xce/icons")

        # import all buttons
        setup().bottom_box_buttons(xce_object)

        # setup datasource button
        update_from_datasource_button = layout_functions.setup_push_button(
            xce_object, xce_object.datasource_button_dict
        )

        ################################################################################
        #                                                                              #
        #                                 DATASETS BOX                                 #
        #                                                                              #
        ################################################################################

        # setup the run button with push button function
        xce_object.dataset_task_run_button = layout_functions.setup_push_button(
            xce_object, xce_object.dataset_task_run_button_dict
        )

        # setup the task button with push button function
        xce_object.dataset_task_status_button = layout_functions.setup_push_button(
            xce_object, xce_object.dataset_task_status_button_dict
        )

        # array of both button xce_objects to apply to bottom box layout
        dataset_buttons = [
            xce_object.dataset_task_run_button,
            xce_object.dataset_task_status_button,
        ]

        # label for the bottom box layout
        dataset_label = str(
            "<html><img src='"
            + str(icons_directory)
            + "/004-database-configuration.png'></html> Datasets"
        )

        # return the frame and combobox from the bottom box setup function
        (
            frame_dataset_task,
            xce_object.dataset_tasks_combobox,
        ) = layout_functions.bottom_box_setup(
            xce_object,
            dataset_label,
            xce_object.dataset_tasks,
            "XChemToolTips." "dataset_task_tip()",
            dataset_buttons,
            "background: " "rgb(240, 255, 140);",
        )

        # define the combobox and buttons in dictionary key to determine behaviour
        xce_object.workflow_widget_dict["Datasets"] = [
            xce_object.dataset_tasks_combobox,
            xce_object.dataset_task_run_button,
            xce_object.dataset_task_status_button,
        ]

        ################################################################################
        #                                                                              #
        #                             MAPS & RESTRAINTS BOX                            #
        #                                                                              #
        ################################################################################
        # settings for the run button

        # setup the run button with push button function
        xce_object.map_cif_file_task_run_button = layout_functions.setup_push_button(
            xce_object, xce_object.map_cif_file_task_run_button_dict
        )

        # setup the task button with push button function
        xce_object.map_cif_file_task_status_button = layout_functions.setup_push_button(
            xce_object, xce_object.map_cif_file_task_status_button_dict
        )

        # array of both button xce_objects to apply to bottom box layout
        map_cif_file_buttons = [
            xce_object.map_cif_file_task_run_button,
            xce_object.map_cif_file_task_status_button,
        ]

        # label for the bottom box layout
        map_cif_file_label = str(
            "<html><img src='"
            + str(icons_directory)
            + "/003-internet.png'></html> Maps & Restraints"
        )

        # return the frame and combobox from the bottom box setup function
        (
            frame_map_cif_file_task,
            xce_object.map_cif_file_tasks_combobox,
        ) = layout_functions.bottom_box_setup(
            xce_object,
            map_cif_file_label,
            xce_object.map_cif_file_tasks,
            "XChemToolTips.map_cif_file_" "task_tip()",
            map_cif_file_buttons,
            "background: rgb(140, 255, " "150); ",
        )

        # define the combobox and buttons in dictionary key to determine behaviour
        xce_object.workflow_widget_dict["Maps"] = [
            xce_object.map_cif_file_tasks_combobox,
            xce_object.map_cif_file_task_run_button,
            xce_object.map_cif_file_task_status_button,
        ]

        ################################################################################
        #                                                                              #
        #                            HIT IDENTIFICATION BOX                            #
        #                                                                              #
        ################################################################################
        # settings for the run button

        # setup the run button with push button function
        xce_object.panddas_file_task_run_button = layout_functions.setup_push_button(
            xce_object, xce_object.panddas_file_task_run_button_dict
        )

        # setup the task button with push button function
        xce_object.panddas_file_task_status_button = layout_functions.setup_push_button(
            xce_object, xce_object.panddas_file_task_status_button_dict
        )

        # array of both button xce_objects to apply to bottom box layout
        panddas_file_buttons = [
            xce_object.panddas_file_task_run_button,
            xce_object.panddas_file_task_status_button,
        ]

        # label for the bottom box layout
        panddas_file_label = str(
            "<html><img src='"
            + str(icons_directory)
            + "/002-chinese-panda-bear.png'></html> Hit Identification"
        )

        # return the frame and combobox from the bottom box setup function
        (
            frame_panddas_file_task,
            xce_object.panddas_file_tasks_combobox,
        ) = layout_functions.bottom_box_setup(
            xce_object,
            panddas_file_label,
            xce_object.panddas_file_tasks,
            "XChemToolTips.panddas_file_" "task_tip()",
            panddas_file_buttons,
            "background: rgb(140,200,255)" "; ",
        )

        # define the combobox and buttons in dictionary key to determine behaviour
        xce_object.workflow_widget_dict["PANDDAs"] = [
            xce_object.panddas_file_tasks_combobox,
            xce_object.panddas_file_task_run_button,
            xce_object.panddas_file_task_status_button,
        ]

        ################################################################################
        #                                                                              #
        #                                REFINEMENT BOX                                #
        #                                                                              #
        ################################################################################
        # settings for the run button

        # setup the run button with push button function
        xce_object.refine_file_task_run_button = layout_functions.setup_push_button(
            xce_object, xce_object.refine_file_task_run_button_dict
        )

        # setup the task button with push button function
        xce_object.refine_file_task_status_button = layout_functions.setup_push_button(
            xce_object, xce_object.refine_file_task_status_button_dict
        )

        # array of both button xce_objects to apply to bottom box layout
        refine_file_buttons = [
            xce_object.refine_file_task_run_button,
            xce_object.refine_file_task_status_button,
        ]

        # label for the bottom box layout
        refine_file_label = str(
            "<html><img src='"
            + str(icons_directory)
            + "/001-ducky.png'></html> Refinement"
        )

        # return the frame and combobox from the bottom box setup function
        (
            frame_refine_file_task,
            xce_object.refine_file_tasks_combobox,
        ) = layout_functions.bottom_box_setup(
            xce_object,
            refine_file_label,
            xce_object.refine_file_tasks,
            "XChemToolTips.refine_file_task" "_tip()",
            refine_file_buttons,
            "background: rgb(245, 190, 255)" ";",
        )

        # define the combobox and buttons in dictionary key to determine behaviour
        xce_object.workflow_widget_dict["Refinement"] = [
            xce_object.refine_file_tasks_combobox,
            xce_object.refine_file_task_run_button,
            xce_object.refine_file_task_status_button,
        ]

        return (
            update_from_datasource_button,
            frame_dataset_task,
            frame_map_cif_file_task,
            frame_panddas_file_task,
            frame_refine_file_task,
        )

    def main_layout(self, xce_object):
        # initialise menu bar
        menu_bar = self.initialise_menu_bar(xce_object)

        # initialise bottom boxes
        (
            update_from_datasource_button,
            frame_dataset_task,
            frame_map_cif_file_task,
            frame_panddas_file_task,
            frame_refine_file_task,
        ) = self.initialise_bottom_boxes(xce_object)

        # Tab layout & content
        # --------------------
        #
        # Overview
        # |- datasource - TABLE
        # |- summary - GRAPH
        #
        # Datasets
        # |- summary - TABLE
        #
        # Maps - TABLE
        #
        # PANDDAS
        # |- pandda.analyse - TABLE
        # |- Dataset Summary ------------------
        # |- Processing Output                  |   HTML
        # |- pandda.inspect                     |
        # |- Statistical Map Summaries --------
        #
        # Refinement - TABLE
        #
        # Deposition
        #
        # Settings

        # Setup tabs
        OverviewTab().setup(xce_object)
        DatasetsTab().setup(xce_object)
        MapsTab().setup(xce_object)
        PanddaTab().setup(xce_object)
        RefinementTab().setup(xce_object)
        DepositionTab().setup(xce_object)
        SettingsTab().setup(xce_object)

        ################################################################################
        #                                                                              #
        #                                  STATUS BAR                                  #
        #                                                                              #
        ################################################################################
        xce_object.status_bar = QtGui.QStatusBar()
        xce_object.progress_bar = QtGui.QProgressBar()
        xce_object.progress_bar.setMaximum(100)
        xce_object.status_bar.setMaximumWidth(xce_object.screen.width())
        xce_object.progress_bar.setMaximumWidth(xce_object.screen.width())
        hbox_status = QtGui.QHBoxLayout()
        hbox_status.addWidget(xce_object.status_bar)
        hbox_status.addWidget(xce_object.progress_bar)

        vbox_main = QtGui.QVBoxLayout()
        menu_bar.setMaximumWidth(xce_object.screen.width())
        vbox_main.addWidget(menu_bar)
        xce_object.main_tab_widget.setMaximumSize(
            xce_object.screen.width(), xce_object.screen.height() - 245
        )
        vbox_main.addWidget(xce_object.main_tab_widget)

        hboxTaskFrames = QtGui.QHBoxLayout()

        hboxTaskFrames.addWidget(update_from_datasource_button)
        hboxTaskFrames.addWidget(frame_dataset_task)
        hboxTaskFrames.addWidget(frame_map_cif_file_task)
        hboxTaskFrames.addWidget(frame_panddas_file_task)
        hboxTaskFrames.addWidget(frame_refine_file_task)

        vbox_main.addLayout(hboxTaskFrames)

        vbox_main.addLayout(hbox_status)

        xce_object.window.setLayout(vbox_main)

        xce_object.status_bar.showMessage("Ready")
        xce_object.window.show()

        if xce_object.data_source_file != "":
            write_enabled = xce_object.check_write_permissions_of_data_source()
            if not write_enabled:
                xce_object.data_source_set = False

    def workflow(self, xce_object):
        ################################################################################
        #                                                                              #
        # ========================== WORKFLOW TASK CONTAINER ========================= #
        #                                                                              #
        ################################################################################

        # workflow task container - order of tabs as they appear for the main window
        xce_object.workflow = [
            "Overview",  # 0
            "Datasets",  # 1
            "Maps",  # 2
            "PANDDAs",  # 3
            "Refinement",  # 4
            "Deposition",  # 6
            "Settings",
        ]  # 5

        # dictionary with keys corresponding to each stage in the workflow
        xce_object.workflow_dict = {
            xce_object.workflow[0]: "Overview",
            xce_object.workflow[1]: "Datasets",
            xce_object.workflow[2]: "Maps",
            xce_object.workflow[3]: "PANDDAs",
            xce_object.workflow[4]: "Refinement",
            xce_object.workflow[6]: "Settings",
            xce_object.workflow[5]: "Deposition",
        }

        xce_object.workflow_widget_dict = {}

        # tab widget
        xce_object.main_tab_widget = QtGui.QTabWidget()
        xce_object.tab_dict = {}
        layout_functions.make_tab_dict(
            xce_object.workflow, xce_object.main_tab_widget, xce_object.tab_dict
        )
