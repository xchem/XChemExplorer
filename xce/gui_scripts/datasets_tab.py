import layout
from PyQt4 import QtGui, QtCore


class DatasetsTab:
    def __init__(self):
        self.layout_funcs = layout.LayoutFuncs()

    def setup(self, xce_object):
        ################################################################################
        #                                                                              #
        #                                 DATASETS TAB                                 #
        #                                                                              #
        ################################################################################
        # define subtab list, widget and dict
        datasets_tab_list = ["Summary"]
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

        xce_object.datasets_data_collection_vbox.addWidget(
            xce_object.datasets_tab_widget
        )  # add subtab to main tab
