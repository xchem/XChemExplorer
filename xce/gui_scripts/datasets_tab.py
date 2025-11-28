from PyQt4 import QtCore, QtGui

from xce.gui_scripts import layout_functions


class DatasetsTab:
    def setup(self, xce_object):
        ################################################################################
        #                                                                              #
        #                                 DATASETS TAB                                 #
        #                                                                              #
        ################################################################################

        # main body - things that are always displayed
        # add a container to hold everythting and add to main tab layout
        xce_object.datasets_data_collection_vbox = QtGui.QVBoxLayout()

        # add a horizontal box to hold option to autocheck for new data
        xce_object.autocheck_hbox = QtGui.QHBoxLayout()

        # LEGACY: This feature is disabled because it relies on unsupported functionality - removed from GUI to prevent issues(2025-11-28)
        # checkbox for autocollect
        # xce_object.check_for_new_data_collection = QtGui.QCheckBox(
        #     "Check for new data collection every two minutes"
        # )
        # layout_functions.add_checkbox(
        #     xce_object,
        #     xce_object.check_for_new_data_collection,
        #     "xce_object.continously_check_for_new_data_collection",
        # )

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
            # LEGACY: check_for_new_data_collection checkbox removed (2025-11-28)
            select_target_label,
            xce_object.target_selection_combobox,
        ]

        layout_functions.add_to_box(
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

        layout_functions.table_setup(
            xce_object.datasets_summary_table, xce_object.datasets_summary_table_columns
        )

        xce_object.datasets_data_collection_vbox.addWidget(
            xce_object.datasets_summary_table
        )  # add subtab to main tab
