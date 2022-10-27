import XChemMain
import layout
import sys
import os
from PyQt4 import QtGui, QtCore, QtWebKit

sys.path.append(os.path.join(os.getenv("XChemExplorer_DIR"), "gui_scripts"))
sys.path.append(os.path.join(os.getenv("XChemExplorer_DIR"), "lib"))


class DatasetsTab:
    def __init__(self):
        self.layout_funcs = layout.LayoutFuncs()

    def setup(self, xce_object):
        ################################################################################################################
        #                                                                                                              #
        #                                                 DATASETS TAB                                                 #
        #                                                                                                              #
        ################################################################################################################
        # define subtab list, widget and dict
        datasets_tab_list = ["Summary", "Reprocess"]
        xce_object.datasets_tab_widget = QtGui.QTabWidget()
        xce_object.datasets_tab_dict = {}

        # make subtabs
        self.layout_funcs.make_tab_dict(
            datasets_tab_list,
            xce_object.datasets_tab_widget,
            xce_object.datasets_tab_dict,
        )

        # main body - things that are always displayed
        # add a container to hold everythting and add to main tab layout
        xce_object.datasets_data_collection_vbox = QtGui.QVBoxLayout()

        # add a horizontal box to hold option to autocheck for new data
        xce_object.autocheck_hbox = QtGui.QHBoxLayout()

        # checkbox for autocollect
        xce_object.check_for_new_data_collection = QtGui.QCheckBox(
            "Check for new data collection every two minutes"
        )
        xce_object.check_for_new_data_collection = QtGui.QCheckBox(
            "Check for new data collection every two minutes"
        )
        self.layout_funcs.add_checkbox(
            xce_object,
            xce_object.check_for_new_data_collection,
            "xce_object.continously_check_for_new_data_collection",
        )

        # select target dropdown
        select_target_label = QtGui.QLabel("<b>Select Target: </b>")
        select_target_label.setAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)
        xce_object.target_selection_combobox = QtGui.QComboBox()
        xce_object.populate_target_selection_combobox(
            xce_object.target_selection_combobox
        )
        xce_object.target_selection_combobox.activated[str].connect(
            xce_object.target_selection_combobox_activated
        )
        xce_object.target = str(xce_object.target_selection_combobox.currentText())

        # array defining order of xce_objects to add
        xce_object.autocheck_hbox_widgets = [
            xce_object.check_for_new_data_collection,
            select_target_label,
            xce_object.target_selection_combobox,
        ]

        self.layout_funcs.add_to_box(
            xce_object.autocheck_hbox, xce_object.autocheck_hbox_widgets
        )  # add xce_objects in order

        # add target dropdown to top bar
        xce_object.datasets_data_collection_vbox.addLayout(xce_object.autocheck_hbox)

        # summary sub-tab
        # table
        xce_object.datasets_summary_table = QtGui.QTableWidget()

        xce_object.datasets_summary_table.resizeRowsToContents()
        xce_object.datasets_summary_table.resizeColumnsToContents()
        xce_object.datasets_summary_table.setSelectionBehavior(
            QtGui.QAbstractItemView.SelectRows
        )
        xce_object.datasets_summary_table.cellClicked.connect(
            xce_object.show_results_from_all_pipelines
        )

        self.layout_funcs.table_setup(
            xce_object.datasets_summary_table, xce_object.datasets_summary_table_columns
        )
        # setup layout to hold table
        xce_object.datasets_summarys_vbox_for_table = QtGui.QVBoxLayout()
        xce_object.datasets_summarys_vbox_for_table.addWidget(
            xce_object.datasets_summary_table
        )  # add table to layout
        xce_object.datasets_summarys_vbox_for_details = (
            QtGui.QVBoxLayout()
        )  # vbox for details
        # switch for displaying/updating table
        xce_object.data_collection_details_currently_on_display = None

        xce_object.datasets_data_collection_vbox.addWidget(
            xce_object.datasets_tab_widget
        )  # add subtab to main tab

        # reprocessing sub-tab
        # top options
        # data collection label
        dc_label = QtGui.QLabel("<b>Data collection directory: </b>")
        # align left and centre of container
        dc_label.setAlignment(QtCore.Qt.AlignLeft | QtCore.Qt.AlignVCenter)
        xce_object.diffraction_data_dir_label = QtGui.QLabel(
            xce_object.diffraction_data_directory
        )  # add directory as text
        xce_object.diffraction_data_dir_label.setAlignment(
            QtCore.Qt.AlignLeft | QtCore.Qt.AlignVCenter
        )  # align as above

        # select label
        select_button = QtGui.QPushButton("Select")
        # attach file open dialogue
        select_button.clicked.connect(xce_object.select_diffraction_data_directory)

        # search button
        search_button = QtGui.QPushButton("Search Datasets")
        # search for datasets in the selected directory
        search_button.clicked.connect(xce_object.search_for_datasets)

        # search info
        xce_object.diffraction_data_search_label = QtGui.QLabel(
            xce_object.diffraction_data_search_info
        )

        # translate label
        translate_label = QtGui.QLabel("<b>Translate:</b> datasetID -> sampleID")
        # align in centre of container
        translate_label.setAlignment(QtCore.Qt.AlignCenter)

        # CSV button
        csv_button = QtGui.QPushButton("Open CSV")
        csv_button.setStyleSheet("QPushButton { padding: 1px; margin: 1px }")
        # open the relevant csv file
        csv_button.clicked.connect(xce_object.translate_datasetID_to_sampleID)

        # create hbox to hold everything and add widgets to it
        xce_object.hbox_select = QtGui.QHBoxLayout()  # top options box
        xce_object.hbox_select_widgets = [
            dc_label,
            xce_object.diffraction_data_dir_label,
            select_button,
            search_button,
            xce_object.diffraction_data_search_label,
            translate_label,
            csv_button,
        ]  # array defining order of xce_objects to be added
        self.layout_funcs.add_to_box(
            xce_object.hbox_select, xce_object.hbox_select_widgets
        )  # add xce_objects in order

        # frame to hold everything
        frame_select = QtGui.QFrame()
        # apply to containing frame
        frame_select.setLayout(xce_object.hbox_select)

        # table - main body
        xce_object.datasets_reprocess_table = QtGui.QTableWidget()
        self.layout_funcs.table_setup(
            xce_object.datasets_reprocess_table, xce_object.datasets_reprocess_columns
        )  # setup

        # create context menu - no idea where this lives...
        xce_object.popMenu_for_datasets_reprocess_table = QtGui.QMenu()
        run_xia2_on_selected = QtGui.QAction(
            "mark selected for reprocessing", xce_object.window
        )
        run_xia2_on_selected.triggered.connect(xce_object.select_sample_for_xia2)
        xce_object.popMenu_for_datasets_reprocess_table.addAction(run_xia2_on_selected)
        xce_object.datasets_reprocess_table.setContextMenuPolicy(
            QtCore.Qt.CustomContextMenu
        )
        xce_object.datasets_reprocess_table.customContextMenuRequested.connect(
            xce_object.on_context_menu_reprocess_data
        )

        # options at bottom of tab
        # data processing label
        label = QtGui.QLabel("<b>Data processing protocol: </b>")
        label.setAlignment(QtCore.Qt.AlignLeft | QtCore.Qt.AlignVCenter)

        # option checkboxes
        xce_object.xia2_3d_checkbox = QtGui.QCheckBox("xia2 3d")
        xce_object.xia2_3dii_checkbox = QtGui.QCheckBox("xia2 3dii")
        xce_object.xia2_dials_checkbox = QtGui.QCheckBox("Dials")

        # spacegroup label
        sg_label = QtGui.QLabel("<b>Space Group:</b>")

        # spacegroup dropdown menu
        xce_object.reprocess_space_group_comboxbox = QtGui.QComboBox()
        xce_object.reprocess_space_group_comboxbox.addItem("ignore")
        for sg in XChemMain.space_group_list():
            xce_object.reprocess_space_group_comboxbox.addItem(sg)

        # mtz label
        mtz_label = QtGui.QLabel("<b>Reference MTZ:</b>")

        # file label
        xce_object.reprocess_reference_mtz_file_label = QtGui.QLabel(
            xce_object.diffraction_data_reference_mtz
        )

        # select button
        select_button = QtGui.QPushButton("Select")
        select_button.clicked.connect(xce_object.select_reprocess_reference_mtz)

        # define order of widgets to be added to options hbox
        hbox_options_widgets = [
            label,
            xce_object.xia2_3d_checkbox,
            xce_object.xia2_3dii_checkbox,
            xce_object.xia2_dials_checkbox,
            sg_label,
            xce_object.reprocess_space_group_comboxbox,
            mtz_label,
            xce_object.reprocess_reference_mtz_file_label,
            select_button,
        ]

        # create hbox, add everything to it and then put it in a frame
        hbox_options = QtGui.QHBoxLayout()
        self.layout_funcs.add_to_box(hbox_options, hbox_options_widgets)

        frame_options = QtGui.QFrame()
        frame_options.setLayout(hbox_options)

        # following are contained in vboxes
        # res limit isig label
        label = QtGui.QLabel("<b>Res.\nLimit:</b>\nMn<I/sig(I)>")
        label.setAlignment(QtCore.Qt.AlignCenter)

        # res limit isig dropdown menu
        xce_object.reprocess_isigma_combobox = QtGui.QComboBox()
        misigma = ["default", "4", "3", "2.5", "2", "1.5", "1", "0.5"]
        self.layout_funcs.populate_combobox(
            misigma, xce_object.reprocess_isigma_combobox
        )
        xce_object.reprocess_isigma_combobox.setCurrentIndex(0)
        xce_object.reprocess_isigma_combobox.setStyleSheet(
            " QComboBox { padding: 1px; margin: 1px }"
        )

        # create vertical box to add labels and dropdowns to, create box and put in frame
        vbox_isigma = QtGui.QVBoxLayout()
        vbox_isigma_widgets = [label, xce_object.reprocess_isigma_combobox]
        self.layout_funcs.add_to_box(vbox_isigma, vbox_isigma_widgets)
        frame_isigma = QtGui.QFrame()
        frame_isigma.setLayout(vbox_isigma)

        # res limit cc half label
        res_cc_label = QtGui.QLabel("<b>Res.\nLimit:</b>\nCC 1/2")
        res_cc_label.setAlignment(QtCore.Qt.AlignCenter)

        # res limit cc half dropdown
        xce_object.reprocess_cc_half_combobox = QtGui.QComboBox()
        cc_half = [
            "default",
            "0.9",
            "0.8",
            "0.7",
            "0.6",
            "0.5",
            "0.4",
            "0.3",
            "0.2",
            "0.1",
        ]
        self.layout_funcs.populate_combobox(
            cc_half, xce_object.reprocess_cc_half_combobox
        )
        xce_object.reprocess_cc_half_combobox.setCurrentIndex(0)
        xce_object.reprocess_cc_half_combobox.setStyleSheet(
            " QComboBox { padding: 1px; margin: 1px }"
        )

        # create a vbox for label and dropdown, and add items to it
        vbox_cc_half = QtGui.QVBoxLayout()
        vbox_cc_half_widgets = [res_cc_label, xce_object.reprocess_cc_half_combobox]
        self.layout_funcs.add_to_box(vbox_cc_half, vbox_cc_half_widgets)

        # create frame to hold everything and add vbox
        frame_cc_half = QtGui.QFrame()
        frame_cc_half.setLayout(vbox_cc_half)

        # create a hbox to hold the bottom frames and add everything
        data_protocol_hbox = QtGui.QHBoxLayout()
        data_protocol_hbox_widgets = [frame_options, frame_isigma, frame_cc_half]
        self.layout_funcs.add_to_box(data_protocol_hbox, data_protocol_hbox_widgets)

        # create frame to hold everything (horizontal)
        bottom_options_frame = QtGui.QFrame()
        bottom_options_frame.setLayout(data_protocol_hbox)

        # code below sets final layout for whole subtab
        # box to hold reprocessing subtab content
        xce_object.reprocess_vbox = QtGui.QVBoxLayout()
        xce_object.reprocess_hbox_widgets = [
            frame_select,
            xce_object.datasets_reprocess_table,
            bottom_options_frame,
        ]
        self.layout_funcs.add_to_box(
            xce_object.reprocess_vbox, xce_object.reprocess_hbox_widgets
        )
