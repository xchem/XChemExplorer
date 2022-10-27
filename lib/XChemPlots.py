# last edited: 29/07/2016, 15:00

import numpy as np
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
import matplotlib.pyplot as plt
import XChemMain
import os
import sys

sys.path.append(os.path.join(os.getenv('XChemExplorer_DIR'), 'lib'))


class summary_plot(object):
    def __init__(self, datasource, overview_axes):
        self.datasource = datasource
        self.overview_axes = overview_axes

    def update_overview(self):

        db_dict = XChemMain.get_datasource_summary(self.datasource)

        category = ('Crystals mounted',
                    'Data Collection',
                    'Maps',
                    'PANDDA',
                    'Refinement',
                    'Comp Chem')

        N = len(category)
        Success = np.array([db_dict['nr_samples'],
                            db_dict['nr_data_collection_success'],
                            db_dict['nr_initial_maps_available'],
                            db_dict['nr_pandda_hits'],
                            5,
                            6])

        Processed = np.array([0,
                              0,
                              0,
                              db_dict['nr_pandda_processed'],
                              0,
                              0])

        Pending = np.array([0,
                            db_dict['nr_data_collection_pending'],
                            db_dict['nr_initial_maps_pending'],
                            db_dict['nr_pandda_pending'],
                            2,
                            1])

        Failure = np.array([db_dict['nr_samples_failed_to_mount'],
                            db_dict['nr_data_collection_failed'],
                            db_dict['nr_initial_maps_fail'],
                            db_dict['nr_pandda_reject'],
                            2,
                            3])

        ind = np.arange(N)  # the x locations for the groups
        width = 0.35       # the width of the bars

        p0 = self.overview_axes.bar(ind, Success,   width, color='g')
        p1 = self.overview_axes.bar(
            ind, Processed, width, color='b', bottom=Success)
        p2 = self.overview_axes.bar(
            ind, Pending,   width, color='y', bottom=Success+Processed)
        p3 = self.overview_axes.bar(
            ind, Failure,   width, color='r', bottom=Pending+Success+Processed)

        # add some text for labels, title and axes ticks
        self.overview_axes.set_ylabel('N')
        self.overview_axes.set_title('Overview')
        self.overview_axes.set_xticks(ind + width)
        self.overview_axes.set_xticklabels(category)

        self.overview_axes.legend(
            (p0[0], p1[0], p2[0], p3[0]), ('Success', 'Processed', 'Pending', 'Failure'))
