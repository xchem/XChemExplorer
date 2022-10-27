import layout
from PyQt4 import QtGui


class RefinementTab:
    def __init__(self):
        self.layout_funcs = layout.LayoutFuncs()

    def setup(self, xce_object):
        ################################################################################################################
        #                                                                                                              #
        #                                                REFINEMENT TAB                                                #
        #                                                                                                              #
        ################################################################################################################
        xce_object.summary_vbox_for_table = QtGui.QVBoxLayout()

        # table
        xce_object.refinement_table = QtGui.QTableWidget()
        self.layout_funcs.table_setup(
            xce_object.refinement_table, xce_object.refinement_table_columns
        )
        xce_object.summary_vbox_for_table.addWidget(xce_object.refinement_table)
