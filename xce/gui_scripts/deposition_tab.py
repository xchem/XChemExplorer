from PyQt4 import QtGui

from xce.gui_scripts import layout_functions

from xce.lib import XChemToolTips


class DepositionTab:
    def setup(self, xce_object):
        ################################################################################
        #                                                                              #
        #                                DEPOSITION TAB                                #
        #                                                                              #
        ################################################################################
        xce_object.deposition_vbox = QtGui.QVBoxLayout()

        scroll = QtGui.QScrollArea()
        xce_object.deposition_vbox.addWidget(scroll)
        scrollContent = QtGui.QWidget(scroll)
        scrollLayout = QtGui.QVBoxLayout(scrollContent)
        scrollContent.setLayout(scrollLayout)

        # deposition page heading
        deposition_page_heading = layout_functions.add_depo_heading(
            "Group deposition of bound-state structures & ground-state model"
        )
        deposition_page_heading.setStyleSheet("font: bold 40pt Arial")

        deposition_page_introduction = QtGui.QLabel(
            XChemToolTips.deposition_introduction()
        )
        deposition_page_introduction.setWordWrap(True)
        deposition_page_introduction.setStyleSheet("font: 16pt Arial")
        deposition_page_introduction.setOpenExternalLinks(True)

        # bound-state depostion

        deposition_bound_state_heading = layout_functions.add_depo_heading(
            "Group deposition of bound-state structures"
        )
        deposition_bound_state_heading.setStyleSheet("font: bold 20pt Arial")

        deposition_bound_state_prerequisites = layout_functions.add_depo_heading(
            "Prerequisites"
        )
        deposition_bound_state_prerequisites.setStyleSheet(
            "font: italic bold 17pt Arial"
        )

        deposition_bound_state_prerequisites_text = layout_functions.add_depo_text(
            XChemToolTips.deposition_bound_state_prerequisites()
        )

        deposition_bound_state_preparation = layout_functions.add_depo_heading(
            "Procedure"
        )
        deposition_bound_state_preparation.setStyleSheet(
            "font: italic bold 17pt Arial "
        )

        deposition_bound_state_preparation_step_one_text = (
            layout_functions.add_depo_text(
                XChemToolTips.deposition_bound_state_preparation_step_one_text()
            )
        )

        xce_object.deposition_bounnd_state_preparation_ignore_event_map = (
            QtGui.QCheckBox(
                XChemToolTips.deposition_bounnd_state_preparation_ignore_event_map()
            )
        )

        prepare_mmcif_button = QtGui.QPushButton("prepare mmcif")
        prepare_mmcif_button.clicked.connect(
            xce_object.prepare_models_for_deposition_ligand_bound
        )
        prepare_mmcif_button.setMaximumWidth(200)

        deposition_bound_state_preparation_step_two_text = (
            layout_functions.add_depo_text(
                XChemToolTips.deposition_bound_state_preparation_step_two_text()
            )
        )

        copy_mmcif_button = QtGui.QPushButton("copy mmcif")
        copy_mmcif_button.clicked.connect(
            xce_object.prepare_for_group_deposition_upload_ligand_bound
        )
        copy_mmcif_button.setMaximumWidth(200)

        pdb_group_deposition_instruction_one = layout_functions.add_depo_text(
            XChemToolTips.pdb_group_deposition_instruction_one()
        )

        pdb_group_deposition_link = QtGui.QLabel(
            XChemToolTips.pdb_group_deposition_link()
        )
        pdb_group_deposition_link.setOpenExternalLinks(True)

        pdb_group_deposition_link_two = QtGui.QLabel(
            XChemToolTips.pdb_group_deposition_link()
        )
        pdb_group_deposition_link_two.setOpenExternalLinks(True)

        pdb_group_deposition_instruction_two = layout_functions.add_depo_text(
            XChemToolTips.pdb_group_deposition_instruction_two()
        )
        pdb_group_deposition_instruction_two_two = layout_functions.add_depo_text(
            XChemToolTips.pdb_group_deposition_instruction_two()
        )

        # ground-state depostion

        deposition_ground_state_heading = layout_functions.add_depo_heading(
            "Group deposition of ground-state model"
        )
        deposition_ground_state_heading.setStyleSheet("font: bold 20pt Arial")

        deposition_ground_state_preparation_step_one_text = (
            layout_functions.add_depo_text(
                XChemToolTips.notification_about_changes_to_apo_deposition()
            )
        )

        deposition_ground_state_prerequisites = layout_functions.add_depo_heading(
            "Prerequisites"
        )
        deposition_ground_state_prerequisites.setStyleSheet(
            "font: italic bold 17pt Arial"
        )

        deposition_ground_state_prerequisites_text = layout_functions.add_depo_text(
            XChemToolTips.deposition_ground_state_prerequisites()
        )

        deposition_ground_state_preparation = layout_functions.add_depo_heading(
            "Procedure"
        )
        deposition_ground_state_preparation.setStyleSheet(
            "font: italic bold 17pt Arial "
        )

        deposition_ground_state_preparation_step_three_text = (
            layout_functions.add_depo_text(
                XChemToolTips.deposition_ground_state_preparation_step_three_text()
            )
        )
        xce_object.ground_state_pandda_directory_label = QtGui.QLabel(
            xce_object.panddas_directory
        )
        xce_object.ground_state_pandda_directory_label.setStyleSheet("color: blue")

        deposition_ground_state_preparation_step_four_text = (
            layout_functions.add_depo_text(
                XChemToolTips.deposition_ground_state_preparation_step_four_text()
            )
        )

        add_ground_state_db_button = QtGui.QPushButton("Add to database")
        add_ground_state_db_button.clicked.connect(xce_object.add_ground_state_db)
        add_ground_state_db_button.setMaximumWidth(200)

        deposition_ground_state_preparation_step_five_text = (
            layout_functions.add_depo_text(
                XChemToolTips.deposition_ground_state_preparation_step_five_text()
            )
        )

        deposition_ground_state_preparation_step_six_text = (
            layout_functions.add_depo_text(
                XChemToolTips.deposition_ground_state_preparation_step_six_text()
            )
        )

        prepare_ground_state_mmcif_button = QtGui.QPushButton("Prepare mmcif")
        prepare_ground_state_mmcif_button.clicked.connect(
            xce_object.prepare_ground_state_mmcif
        )
        prepare_ground_state_mmcif_button.setMaximumWidth(200)

        deposition_ground_state_preparation_step_seven_text = (
            layout_functions.add_depo_text(
                XChemToolTips.deposition_ground_state_preparation_step_seven_text()
            )
        )

        copy_apo_mmcif_button = QtGui.QPushButton("copy mmcif")
        copy_apo_mmcif_button.clicked.connect(
            xce_object.prepare_for_group_deposition_upload_ground_state
        )
        copy_apo_mmcif_button.setMaximumWidth(200)

        deposition_ground_state_preparation_step_eight_text = (
            layout_functions.add_depo_text(
                XChemToolTips.deposition_ground_state_preparation_step_eight_text()
            )
        )

        # after ligand_bound depostion
        after_deposition_heading = layout_functions.add_depo_heading(
            "After deposition of ligand-bound structures"
        )
        after_deposition_heading.setStyleSheet("font: bold 20pt Arial")

        after_deposition_preparation = layout_functions.add_depo_heading("Procedure")
        after_deposition_preparation.setStyleSheet("font: italic bold 17pt Arial ")
        after_deposition_preparation_text = layout_functions.add_depo_text(
            XChemToolTips.after_deposition_step_one_text()
        )

        ###############################################

        deposition_widget_list = [
            deposition_page_heading,
            QtGui.QLabel(" \n "),
            deposition_page_introduction,
            QtGui.QLabel(" \n "),
            deposition_bound_state_heading,
            QtGui.QLabel(" \n "),
            deposition_bound_state_prerequisites,
            deposition_bound_state_prerequisites_text,
            QtGui.QLabel(" \n "),
            deposition_bound_state_preparation,
            deposition_bound_state_preparation_step_one_text,
            xce_object.deposition_bounnd_state_preparation_ignore_event_map,
            prepare_mmcif_button,
            deposition_bound_state_preparation_step_two_text,
            copy_mmcif_button,
            pdb_group_deposition_instruction_one,
            pdb_group_deposition_link,
            pdb_group_deposition_instruction_two_two,
            QtGui.QLabel(" \n\n\n "),
            after_deposition_heading,
            QtGui.QLabel(" \n "),
            after_deposition_preparation,
            after_deposition_preparation_text,
            QtGui.QLabel(" \n\n\n "),
            deposition_ground_state_heading,
            QtGui.QLabel(" \n "),
            deposition_ground_state_preparation_step_one_text,
            QtGui.QLabel(" \n "),
            deposition_ground_state_prerequisites,
            deposition_ground_state_prerequisites_text,
            QtGui.QLabel(" \n "),
            deposition_ground_state_preparation,
            deposition_ground_state_preparation_step_three_text,
            xce_object.ground_state_pandda_directory_label,
            deposition_ground_state_preparation_step_four_text,
            add_ground_state_db_button,
            deposition_ground_state_preparation_step_five_text,
            deposition_ground_state_preparation_step_six_text,
            prepare_ground_state_mmcif_button,
            deposition_ground_state_preparation_step_seven_text,
            copy_apo_mmcif_button,
            deposition_ground_state_preparation_step_eight_text,
            pdb_group_deposition_link_two,
            pdb_group_deposition_instruction_two,
            QtGui.QLabel(" \n\n\n "),
        ]

        layout_functions.add_to_box(scrollLayout, deposition_widget_list)

        # container settings
        scrollLayout.addStretch(1)
        scroll.setWidget(scrollContent)
