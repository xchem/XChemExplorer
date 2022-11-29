import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt4agg import \
    FigureCanvasQTAgg as FigureCanvas
from PyQt4 import QtGui

import layout


class OverviewTab:
    def __init__(self):
        self.layout_funcs = layout.LayoutFuncs()

    def setup(self, xce_object):
        ################################################################################
        #                                                                              #
        #                                 OVERVIEW TAB                                 #
        #                                                                              #
        ################################################################################
        # define subtab list, widget and dict
        overview_tab_list = ["Data Source", "Summary"]
        xce_object.overview_tab_widget = QtGui.QTabWidget()
        xce_object.overview_tab_dict = {}

        # make subtabs
        self.layout_funcs.make_tab_dict(
            overview_tab_list,
            xce_object.overview_tab_widget,
            xce_object.overview_tab_dict,
        )

        # initiate the table in overview/datasource
        xce_object.overview_datasource_table = QtGui.QTableWidget()
        xce_object.overview_datasource_table.setSortingEnabled(True)
        xce_object.overview_datasource_table.resizeColumnsToContents()

        # initiate the graph in overview/summary
        xce_object.overview_figure, xce_object.overview_axes = plt.subplots()
        xce_object.overview_canvas = FigureCanvas(xce_object.overview_figure)
        xce_object.update_summary_plot()
