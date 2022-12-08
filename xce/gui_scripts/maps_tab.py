from PyQt4 import QtCore, QtGui

from xce.gui_scripts import layout


class MapsTab:
    def __init__(self):
        self.layout_funcs = layout.LayoutFuncs()

    def setup(self, xce_object):
        ################################################################################
        #                                                                              #
        #                                   MAPS TAB                                   #
        #                                                                              #
        ################################################################################
        # select box for dimple
        xce_object.select_sample_for_dimple_box = QtGui.QCheckBox(
            "(de-)select all samples for DIMPLE"
        )
        self.layout_funcs.add_checkbox(
            xce_object,
            xce_object.select_sample_for_dimple_box,
            "xce_object.set_run_dimple_flag",
        )

        # set new reference button
        set_new_reference_button = QtGui.QPushButton(
            "Set New Reference (if applicable)"
        )
        set_new_reference_button.clicked.connect(
            xce_object.set_new_reference_if_applicable
        )

        # refresh button
        refresh_reference_file_list_button = QtGui.QPushButton(
            "Refresh reference file list"
        )
        refresh_reference_file_list_button.clicked.connect(
            xce_object.refresh_reference_file_list
        )

        # list and populate reference files
        xce_object.reference_file_list = xce_object.get_reference_file_list(" ")
        xce_object.reference_file_selection_combobox = QtGui.QComboBox()
        xce_object.populate_reference_combobox(
            xce_object.reference_file_selection_combobox
        )

        # setup hbox to hold everything and add widgets
        xce_object.maps_checkbutton_hbox = QtGui.QHBoxLayout()
        maps_checkbutton_widgets = [
            xce_object.select_sample_for_dimple_box,
            set_new_reference_button,
            refresh_reference_file_list_button,
            xce_object.reference_file_selection_combobox,
        ]

        self.layout_funcs.add_to_box(
            xce_object.maps_checkbutton_hbox, maps_checkbutton_widgets
        )

        # table setup
        xce_object.maps_table = QtGui.QTableWidget()
        self.layout_funcs.table_setup(
            xce_object.maps_table, xce_object.maps_table_columns
        )

        # box for table, add to box, add to tab
        xce_object.initial_model_vbox_for_table = QtGui.QVBoxLayout()
        xce_object.initial_model_vbox_for_table.addWidget(xce_object.maps_table)

        # create context menu... no idea where this lives again.
        xce_object.popMenu_for_maps_table = QtGui.QMenu()
        run_dimple = QtGui.QAction("mark selected for dimple run", xce_object.window)
        run_dimple.triggered.connect(xce_object.select_sample_for_dimple)
        xce_object.popMenu_for_maps_table.addAction(run_dimple)
        xce_object.maps_table.setContextMenuPolicy(QtCore.Qt.CustomContextMenu)
        xce_object.maps_table.customContextMenuRequested.connect(
            xce_object.on_context_menu_initial_model
        )
