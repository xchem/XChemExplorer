import sys, os
from PyQt4 import QtGui, QtCore, QtWebKit

sys.path.append(os.path.join(os.getenv('XChemExplorer_DIR'), 'gui_scripts'))
sys.path.append(os.path.join(os.getenv('XChemExplorer_DIR'), 'lib'))

import layout
import XChemToolTips


class DepositionTab():
    def __init__(self):
        self.layout_funcs = layout.LayoutFuncs()

    def setup(self, xce_object):
        ################################################################################################################
        #                                                                                                              #
        #                                                DEPOSITION TAB                                                #
        #                                                                                                              #
        ################################################################################################################
        xce_object.deposition_vbox = QtGui.QVBoxLayout()

        scroll = QtGui.QScrollArea()
        xce_object.deposition_vbox.addWidget(scroll)
        scrollContent = QtGui.QWidget(scroll)
        scrollLayout = QtGui.QVBoxLayout(scrollContent)
        scrollContent.setLayout(scrollLayout)

        # deposition page heading
        deposition_page_heading = self.layout_funcs.add_depo_heading('Group deposition of bound-state structures & ground-state model')
        deposition_page_heading.setStyleSheet("font: bold 40pt Arial")

        deposition_page_introduction = QtGui.QLabel(XChemToolTips.deposition_introduction())
        deposition_page_introduction_link = QtGui.QLabel(XChemToolTips.deposition_introduction_link())
        deposition_page_introduction_link.setOpenExternalLinks(True)

        #
        # bound-state depostion
        #

        deposition_bound_state_heading = self.layout_funcs.add_depo_heading('Group deposition of bound-state structures')
        deposition_bound_state_heading.setStyleSheet("font: bold 20pt Arial")

        deposition_bound_state_prerequisites = self.layout_funcs.add_depo_heading('Prerequisites')
        deposition_bound_state_prerequisites.setStyleSheet("font: italic bold 17pt Arial")

        deposition_bound_state_prerequisites_text = self.layout_funcs.add_depo_text(XChemToolTips.deposition_bound_state_prerequisites())

        deposition_bound_state_preparation = self.layout_funcs.add_depo_heading('Procedure')
        deposition_bound_state_preparation.setStyleSheet("font: italic bold 17pt Arial ")

        xce_object.deposition_bound_state_preparation_step_one_text = self.layout_funcs.add_depo_text(XChemToolTips.deposition_bound_state_preparation_step_one_text())

        deposition_bounnd_state_preparation_ignore_event_map = QtGui.QCheckBox(XChemToolTips.deposition_bounnd_state_preparation_ignore_event_map())

        prepare_mmcif_button = QtGui.QPushButton('prepare mmcif')
        prepare_mmcif_button.clicked.connect(xce_object.prepare_models_for_deposition_ligand_bound)
        prepare_mmcif_button.setMaximumWidth(200)

        deposition_bound_state_preparation_step_two_text = self.layout_funcs.add_depo_text(XChemToolTips.deposition_bound_state_preparation_step_two_text())

        copy_mmcif_button = QtGui.QPushButton('copy mmcif')
        copy_mmcif_button.clicked.connect(xce_object.prepare_for_group_deposition_upload_ligand_bound)
        copy_mmcif_button.setMaximumWidth(200)

        pdb_group_deposition_instruction_one = self.layout_funcs.add_depo_text(XChemToolTips.pdb_group_deposition_instruction_one())

        pdb_group_deposition_link = QtGui.QLabel(XChemToolTips.pdb_group_deposition_link())
        pdb_group_deposition_link.setOpenExternalLinks(True)

        pdb_group_deposition_instruction_two = self.layout_funcs.add_depo_text(XChemToolTips.pdb_group_deposition_instruction_two())


        #
        # ground-state depostion
        #

        deposition_ground_state_heading = self.layout_funcs.add_depo_heading('Group deposition of ground-state model')
        deposition_ground_state_heading.setStyleSheet("font: bold 20pt Arial")



        deposition_ground_state_prerequisites = self.layout_funcs.add_depo_heading('Prerequisites')
        deposition_ground_state_prerequisites.setStyleSheet("font: italic bold 17pt Arial")

        deposition_ground_state_prerequisites_text = self.layout_funcs.add_depo_text(XChemToolTips.deposition_ground_state_prerequisites())

        deposition_ground_state_preparation = self.layout_funcs.add_depo_heading('Procedure')
        deposition_ground_state_preparation.setStyleSheet("font: italic bold 17pt Arial ")

        deposition_ground_state_preparation_step_one_text = self.layout_funcs.add_depo_text(XChemToolTips.deposition_ground_state_preparation_step_one_text())

        ground_state_pdb_button = QtGui.QPushButton('Select PDB file')
        ground_state_pdb_button.clicked.connect(xce_object.select_ground_state_pdb)
        ground_state_pdb_button.setMaximumWidth(200)
        xce_object.ground_state_pdb_button_label = QtGui.QLabel('')

        deposition_ground_state_log_info = self.layout_funcs.add_depo_text(XChemToolTips.deposition_ground_state_log_info())

        deposition_ground_state_preparation_step_two_text = self.layout_funcs.add_depo_text(XChemToolTips.deposition_ground_state_preparation_step_two_text())

        ground_state_mtz_button = QtGui.QPushButton('Select MTZ file')
        ground_state_mtz_button.clicked.connect(xce_object.select_ground_state_mtz)
        ground_state_mtz_button.setMaximumWidth(200)
        xce_object.ground_state_mtz_button_label = QtGui.QLabel('')

        deposition_ground_state_preparation_step_three_text = self.layout_funcs.add_depo_text(XChemToolTips.deposition_ground_state_preparation_step_three_text())

        deposition_ground_state_preparation_step_four_text = self.layout_funcs.add_depo_text(XChemToolTips.deposition_ground_state_preparation_step_four_text())

        add_ground_state_db_button = QtGui.QPushButton('Add to database')
        add_ground_state_db_button.clicked.connect(xce_object.add_ground_state_db)
        add_ground_state_db_button.setMaximumWidth(200)

        deposition_ground_state_preparation_step_five_text = self.layout_funcs.add_depo_text(XChemToolTips.deposition_ground_state_preparation_step_five_text())

        deposition_ground_state_preparation_step_six_text = self.layout_funcs.add_depo_text(XChemToolTips.deposition_ground_state_preparation_step_six_text())

        prepare_ground_state_mmcif_button = QtGui.QPushButton('Prepare mmcif')
        prepare_ground_state_mmcif_button.clicked.connect(xce_object.prepare_ground_state_mmcif)
        prepare_ground_state_mmcif_button.setMaximumWidth(200)

        deposition_ground_state_preparation_step_seven_text = self.layout_funcs.add_depo_text(XChemToolTips.deposition_ground_state_preparation_step_seven_text())

        copy_apo_mmcif_button = QtGui.QPushButton('copy mmcif')
        copy_apo_mmcif_button.clicked.connect(xce_object.prepare_for_group_deposition_upload_ground_state)
        copy_apo_mmcif_button.setMaximumWidth(200)

        deposition_ground_state_preparation_step_eight_text = self.layout_funcs.add_depo_text(XChemToolTips.deposition_ground_state_preparation_step_eight_text())


        #
        # after ligand_bound depostion
        #

        after_deposition_heading = self.layout_funcs.add_depo_heading('After deposition of ligand-bound structures')
        after_deposition_heading.setStyleSheet("font: bold 20pt Arial")



        after_deposition_preparation = self.layout_funcs.add_depo_heading('Procedure')
        after_deposition_preparation.setStyleSheet("font: italic bold 17pt Arial ")
        after_deposition_preparation_text = self.layout_funcs.add_depo_text(XChemToolTips.after_deposition_step_one_text())



###############################################

#        deposition_html_heading = self.layout_funcs.add_depo_heading('HTML export')
#        deposition_html_heading.setStyleSheet("font: 20pt Arial Bold")
#
#        introduction_text = self.layout_funcs.add_depo_text(XChemToolTips.html_summary_introduction())
#
#        html_export_button = QtGui.QPushButton('Export to HTML')
#        html_export_button.clicked.connect(xce_object.export_to_html)
#        html_export_button.setMaximumWidth(200)



        deposition_widget_list = [deposition_page_heading,
                                  QtGui.QLabel(' \n '),
                                  deposition_page_introduction,deposition_page_introduction_link,
                                  QtGui.QLabel(' \n '),

                                  deposition_bound_state_heading, QtGui.QLabel(' \n '),
                                  deposition_bound_state_prerequisites,
                                  deposition_bound_state_prerequisites_text, QtGui.QLabel(' \n '),
                                  deposition_bound_state_preparation,
                                  xce_object.deposition_bound_state_preparation_step_one_text,
                                  deposition_bounnd_state_preparation_ignore_event_map,
                                  prepare_mmcif_button,
                                  deposition_bound_state_preparation_step_two_text,
                                  copy_mmcif_button,
                                  pdb_group_deposition_instruction_one,
                                  pdb_group_deposition_link,
                                  pdb_group_deposition_instruction_two,

                                  QtGui.QLabel(' \n\n\n '),

                                  deposition_ground_state_heading, QtGui.QLabel(' \n '),
                                  deposition_ground_state_prerequisites,
                                  deposition_ground_state_prerequisites_text, QtGui.QLabel(' \n '),
                                  deposition_ground_state_preparation,
                                  deposition_ground_state_preparation_step_one_text,
                                  ground_state_pdb_button,
                                  deposition_ground_state_log_info,
                                  deposition_ground_state_preparation_step_two_text,
                                  ground_state_mtz_button,
                                  deposition_ground_state_preparation_step_three_text,
                                  deposition_ground_state_preparation_step_four_text,
                                  add_ground_state_db_button,
                                  deposition_ground_state_preparation_step_five_text,
                                  deposition_ground_state_preparation_step_six_text,
                                  prepare_ground_state_mmcif_button,
                                  deposition_ground_state_preparation_step_seven_text,
                                  copy_apo_mmcif_button,
                                  deposition_ground_state_preparation_step_eight_text,
                                  pdb_group_deposition_link,
                                  pdb_group_deposition_instruction_two,

                                  QtGui.QLabel(' \n\n\n '),

                                  after_deposition_heading, QtGui.QLabel(' \n '),
                                  after_deposition_preparation,
                                  after_deposition_preparation_text

                                  ]

#                                  QtGui.QLabel(' \n\n\n '),
#
#                                  deposition_html_heading, QtGui.QLabel(' \n '),
#                                  introduction_text,
#                                  html_export_button, QtGui.QLabel('  \n\n '),
#s                                  QtGui.QLabel('  \n  ')]
#
#        deposition_widget_list2 = [update_html_button, QtGui.QLabel('  '),
#                                   upload_html_heading, upload_html_text, QtGui.QLabel('  ')]

        self.layout_funcs.add_to_box(scrollLayout, deposition_widget_list)
#        scrollLayout.addLayout(hbox_zenodo_upload_id)
#        self.layout_funcs.add_to_box(scrollLayout, deposition_widget_list2)

        # container settings
        scrollLayout.addStretch(1)
        scroll.setWidget(scrollContent)

