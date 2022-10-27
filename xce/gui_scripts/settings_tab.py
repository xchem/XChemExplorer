import layout
import os
from PyQt4 import QtGui


class SettingsTab:
    def __init__(self):
        self.layout_funcs = layout.LayoutFuncs()

    def setup(self, xce_object):
        ################################################################################################################
        #                                                                                                              #
        #                                                 SETTINGS TAB                                                 #
        #                                                                                                              #
        ################################################################################################################
        xce_object.settings_container = QtGui.QWidget()
        xce_object.buttons_etc = QtGui.QWidget()
        xce_object.settings_vbox = QtGui.QVBoxLayout()

        xce_object.scroll = QtGui.QScrollArea(xce_object.settings_container)
        xce_object.settings_vbox.addWidget(xce_object.scroll)
        scrollContent_settings = QtGui.QWidget(xce_object.scroll)

        scrollLayout_settings = QtGui.QVBoxLayout(scrollContent_settings)
        scrollContent_settings.setLayout(scrollLayout_settings)

        # Settings Tab
        xce_object.data_collection_vbox_for_settings = QtGui.QVBoxLayout()

        xce_object.buttons_etc.setLayout(xce_object.data_collection_vbox_for_settings)
        xce_object.scroll.setWidget(xce_object.buttons_etc)

        xce_object.initial_model_directory_label = (
            self.layout_funcs.settings_section_setup(
                xce_object.data_collection_vbox_for_settings,
                "\n\n<b>Project Directory: - REQUIRED -</b>",
                xce_object.initial_model_directory,
                "Select Project Directory",
                xce_object.settings_button_clicked,
            )
        )

        xce_object.reference_directory_label = self.layout_funcs.settings_section_setup(
            xce_object.data_collection_vbox_for_settings,
            "\n\n<b>Reference Structure Directory: - OPTIONAL -</b>",
            xce_object.reference_directory,
            "Select Reference Structure Directory",
            xce_object.settings_button_clicked,
        )

        if xce_object.data_source_file != "":
            xce_object.data_source_file_label_text = os.path.join(
                xce_object.database_directory, xce_object.data_source_file
            )
            xce_object.data_source_file_label = (
                self.layout_funcs.settings_section_setup(
                    xce_object.data_collection_vbox_for_settings,
                    "\n\n<b>Data Source: - REQUIRED -</b>",
                    xce_object.data_source_file_label_text,
                    "Select Data Source File",
                    xce_object.settings_button_clicked,
                )
            )
        else:
            xce_object.data_source_file_label_text = ""

            xce_object.data_source_file_label = (
                self.layout_funcs.settings_section_setup(
                    xce_object.data_collection_vbox_for_settings,
                    "\n\n<b>Data Source: - REQUIRED -</b>",
                    xce_object.data_source_file_label_text,
                    "Select Data Source File",
                    xce_object.settings_button_clicked,
                )
            )

        xce_object.data_collection_vbox_for_settings.addWidget(
            QtGui.QLabel(
                "\n\n<b>Data Collection Directory: (e.g. /dls/i04-1/data/2017/lb18145-70) -</b>"
            )
        )

        settings_beamline_frame = QtGui.QFrame()
        settings_beamline_frame.setFrameShape(QtGui.QFrame.StyledPanel)
        settings_beamline_vbox = QtGui.QVBoxLayout()

        settings_hbox_beamline_directory = QtGui.QHBoxLayout()
        xce_object.beamline_directory_label = QtGui.QLabel(
            xce_object.beamline_directory
        )
        settings_hbox_beamline_directory.addWidget(xce_object.beamline_directory_label)
        settings_button_beamline_directory = QtGui.QPushButton(
            "Select Data Collection Directory"
        )
        settings_button_beamline_directory.setMaximumWidth(500)

        settings_button_beamline_directory.clicked.connect(
            xce_object.settings_button_clicked
        )

        settings_hbox_beamline_directory.addWidget(settings_button_beamline_directory)
        settings_beamline_vbox.addLayout(settings_hbox_beamline_directory)
        xce_object.read_agamemnon = QtGui.QCheckBox("Read Agamemnon data structure")
        settings_beamline_vbox.addWidget(xce_object.read_agamemnon)

        #        settings_hbox_datasets_summary_file = QtGui.QHBoxLayout()
        #        xce_object.datasets_summary_file_label = QtGui.QLabel(xce_object.datasets_summary_file)
        #        settings_hbox_datasets_summary_file.addWidget(xce_object.datasets_summary_file_label)
        #        settings_button_datasets_summary_file = QtGui.QPushButton('Select Existing\nCollection Summary File')
        #        settings_button_datasets_summary_file.setMaximumWidth(247)
        #        settings_button_datasets_summary_file.clicked.connect(xce_object.settings_button_clicked)
        #        settings_hbox_datasets_summary_file.addWidget(settings_button_datasets_summary_file)
        #
        #        settings_button_new_datasets_summary_file = QtGui.QPushButton('Assign New\nCollection Summary File')
        #        settings_button_new_datasets_summary_file.clicked.connect(xce_object.settings_button_clicked)
        #        settings_button_new_datasets_summary_file.setMaximumWidth(247)
        #        settings_hbox_datasets_summary_file.addWidget(settings_button_new_datasets_summary_file)
        #
        #        settings_beamline_vbox.addLayout(settings_hbox_datasets_summary_file)

        settings_beamline_frame.setLayout(settings_beamline_vbox)
        xce_object.data_collection_vbox_for_settings.addWidget(settings_beamline_frame)

        xce_object.ccp4_scratch_directory_label = (
            self.layout_funcs.settings_section_setup(
                xce_object.data_collection_vbox_for_settings,
                "\n\n<b>CCP4_SCR Directory: - OPTIONAL -</b>",
                xce_object.ccp4_scratch_directory,
                "Select CCP4_SCR Directory",
                xce_object.settings_button_clicked,
            )
        )

        xce_object.panddas_directory_label = self.layout_funcs.settings_section_setup(
            xce_object.data_collection_vbox_for_settings,
            "\n\n<b>PANDDAs directory: - OPTIONAL -</b>",
            xce_object.panddas_directory,
            "Select PanDDA Directory",
            xce_object.settings_button_clicked,
        )

        xce_object.html_export_directory_label = (
            self.layout_funcs.settings_section_setup(
                xce_object.data_collection_vbox_for_settings,
                "\n\n<b>HTML export directory: - OPTIONAL -</b>",
                xce_object.html_export_directory,
                "Select HTML Export Directory",
                xce_object.settings_button_clicked,
            )
        )

        xce_object.group_deposition_directory_label = (
            self.layout_funcs.settings_section_setup(
                xce_object.data_collection_vbox_for_settings,
                "\n\n<b>Group deposition directory: - OPTIONAL -</b>",
                xce_object.group_deposit_directory,
                "Select Group deposition Directory",
                xce_object.settings_button_clicked,
            )
        )

        # xce_object.data_collection_vbox_for_settings.setSpacing(0)
        xce_object.data_collection_vbox_for_settings.setContentsMargins(30, 30, 30, 30)

        xce_object.buttons_etc.resize(
            xce_object.buttons_etc.sizeHint().width() + 100,
            xce_object.buttons_etc.sizeHint().height(),
        )
