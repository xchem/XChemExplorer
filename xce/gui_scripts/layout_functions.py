from PyQt4 import QtCore, QtGui
import os
from xce.lib import XChemToolTips  # noqa: F401


def make_tab_dict(tab_list, tab_widget, tab_dict):
    for page in tab_list:
        tab = QtGui.QWidget()
        vbox = QtGui.QVBoxLayout(tab)
        tab_widget.addTab(tab, page)
        tab_dict[page] = [tab, vbox]


def add_checkbox(xce_object, checkbox, function, checkopt=False):
    checkbox.toggle()
    checkbox.setChecked(checkopt)
    eval(str("checkbox.stateChanged.connect(" + function + ")"))


def table_setup(table, table_columns, sortingopt=True):
    table.setColumnCount(len(table_columns))
    table.setSortingEnabled(sortingopt)
    table.setHorizontalHeaderLabels(table_columns)
    table.resizeRowsToContents()
    table.resizeColumnsToContents()


def pandda_html(xce_object):
    if os.path.exists(str(xce_object.panddas_directory + "/interesting_datasets")):
        print(
            "WARNING: USING RESULTS FROM OLD PANDDA ANALYSE!"
            " THIS IS NOT FULLY SUPPORTED IN XCE2"
        )
        print(
            "PLEASE CHANGE YOUR PANDDA DIRECTORY TO A NEW RUN,"
            " OR USE THE OLD VERSION OF XCE!"
        )
        xce_object.pandda_initial_html_file = str(
            xce_object.panddas_directory + "/results_summareis/pandda_initial.html"
        )
        xce_object.pandda_analyse_html_file = str(
            xce_object.panddas_directory + "/results_summaries/pandda_analyse.html"
        )
    xce_object.pandda_initial_html_file = str(
        xce_object.panddas_directory
        + "/analyses/html_summaries/"
        + "pandda_initial.html"
    )
    xce_object.pandda_analyse_html_file = str(
        xce_object.panddas_directory
        + "/analyses/html_summaries/"
        + "pandda_analyse.html"
    )
    xce_object.pandda_inspect_html_file = str(
        xce_object.panddas_directory
        + "/analyses/html_summaries/"
        + "pandda_inspect.html"
    )


# function for datasource, run and status button setup
def setup_push_button(xce_object, button_dict):
    # use iterkeys to determine order of key by letter
    for name in sorted(button_dict.keys()):
        # add current item to menu bar
        button = eval('QtGui.QPushButton("' + str(button_dict[name][0]) + '")')
        # for each configuration item
        for button_config in button_dict[name][1]:
            eval(str("button.setToolTip(" + str(button_config[0]) + ")"))
            eval(str('button.setStyleSheet("' + str(button_config[1] + '")')))
            if len(button_config[2]) > 1:
                eval(str("button.setFont(" + str(button_config[2]) + ")"))
            eval(str("button.clicked.connect(" + str(button_config[3]) + ")"))

    return button


# function to setup one of the bottom boxes
def bottom_box_setup(
    xce_object, label, dropdown_options, dropdown_tooltip, buttons, colour
):
    frame = QtGui.QFrame()
    frame.setFrameShape(QtGui.QFrame.StyledPanel)
    frame.setStyleSheet(
        "QFrame {  "
        "border-radius: 1px; padding: 0px; margin: 0px;"
        " background-color: rgb(255, 255, 255); }"
    )

    vbox = QtGui.QVBoxLayout()
    label = QtGui.QLabel(label)
    label.setAlignment(QtCore.Qt.AlignCenter | QtCore.Qt.AlignVCenter)
    label.setStyleSheet(
        str(
            " QLabel { border: 1px solid rgb(184, 192, 210); border-radius: 1px;"
            + str(colour)
            + "padding: 3px; margin: 0px; font: bold 14pt}"
        )
    )
    vbox.addWidget(label)

    hboxAction = QtGui.QHBoxLayout()
    combobox = QtGui.QComboBox()
    for task in dropdown_options:
        combobox.addItem(task)
    eval("combobox.setToolTip(" + str(dropdown_tooltip) + ")")
    combobox.setStyleSheet(" QComboBox { padding: 1px; margin: 1px }")
    hboxAction.addWidget(combobox)

    vboxButton = QtGui.QVBoxLayout()
    for button in buttons:
        vboxButton.addWidget(button)
    hboxAction.addLayout(vboxButton)
    vbox.addLayout(hboxAction)
    vbox.setSpacing(0)
    vbox.setMargin(0)
    frame.setLayout(vbox)
    frame.setMaximumWidth((xce_object.screen.width() - 20) / 5)

    return frame, combobox


# function to add items to top menu bar
def setup_menubar(xce_object, menu_bar, menu_items_dict):
    # use iterkeys to determine order of key by letter
    for config in sorted(menu_items_dict.keys()):
        # add current item to menu bar
        menu = eval('menu_bar.addMenu("' + str(menu_items_dict[config][0]) + '")')
        # for each configuration item
        for menu_item in menu_items_dict[config][1]:
            # add the drop down option
            action = eval(
                str('QtGui.QAction("' + str(menu_item[0]) + '", xce_object.window)')
            )
            # add a shortcut if defined
            if len(menu_item[1]) > 1:
                eval(str('action.setShortcut("' + str(menu_item[1]) + '")'))
            # connect the relevant function and add as an action
            try:
                action.triggered.connect(menu_item[2])
                menu.addAction(action)
            except Exception as exception:
                print((menu_item[2]))
                raise exception

    return menu_bar


def add_to_box(frame, widgets_list):
    for widget in widgets_list:
        frame.addWidget(widget)


def populate_combobox(combobox_list, combobox):
    for item in combobox_list:
        combobox.addItem(item)


def add_depo_heading(heading_text):
    heading = QtGui.QLabel(str(heading_text))
    heading.setStyleSheet("font: bold 20pt Arial")

    return heading


def add_depo_text(text):
    out_text = QtGui.QLabel(text)
    out_text.setStyleSheet("font: 17pt Arial")

    return out_text


def settings_section_setup(vbox, label_text, directory, button_text, button_function):
    vbox.addWidget(QtGui.QLabel(label_text))

    hbox = QtGui.QHBoxLayout()
    directory_label = QtGui.QLabel(directory)
    hbox.addWidget(directory_label)
    button = QtGui.QPushButton(button_text)
    button.setMaximumWidth(500)
    button.clicked.connect(button_function)
    hbox.addWidget(button)

    vbox.addLayout(hbox)
    vbox.addWidget(QtGui.QLabel(" "))
    vbox.addWidget(QtGui.QLabel(" "))

    return directory_label


def add_widgets_layouts(xce_object):
    tab_add_widget = [
        [
            xce_object.tab_dict[xce_object.workflow_dict["Overview"]][1],
            xce_object.overview_tab_widget,
        ],
        [
            xce_object.overview_tab_dict["Data Source"][1],
            xce_object.overview_datasource_table,
        ],
        [
            xce_object.overview_tab_dict["Summary"][1],
            xce_object.overview_canvas,
        ],
        [
            xce_object.pandda_tab_dict["Dataset Summary"][1],
            xce_object.pandda_initial_html,
        ],
        [
            xce_object.pandda_tab_dict["Processing Output"][1],
            xce_object.pandda_analyse_html,
        ],
        [
            xce_object.pandda_tab_dict["pandda.inspect"][1],
            xce_object.pandda_inspect_html,
        ],
    ]

    tab_add_layout = [
        [
            xce_object.tab_dict[xce_object.workflow_dict["Datasets"]][1],
            xce_object.datasets_data_collection_vbox,
        ],
        [
            xce_object.tab_dict[xce_object.workflow_dict["Maps"]][1],
            xce_object.maps_checkbutton_hbox,
        ],
        [
            xce_object.tab_dict[xce_object.workflow_dict["Maps"]][1],
            xce_object.initial_model_vbox_for_table,
        ],
        [
            xce_object.pandda_tab_dict["Statistical Map Summaries"][1],
            xce_object.pandda_map_layout,
        ],
        [
            xce_object.pandda_tab_dict["pandda.analyse"][1],
            xce_object.pandda_analyse_hbox,
        ],
        [
            xce_object.tab_dict[xce_object.workflow_dict["PANDDAs"]][1],
            xce_object.panddas_results_vbox,
        ],
        [
            xce_object.tab_dict[xce_object.workflow_dict["Refinement"]][1],
            xce_object.summary_vbox_for_table,
        ],
        [
            xce_object.tab_dict[xce_object.workflow_dict["Deposition"]][1],
            xce_object.deposition_vbox,
        ],
        [
            xce_object.tab_dict[xce_object.workflow_dict["Settings"]][1],
            xce_object.settings_vbox,
        ],
    ]

    for item in tab_add_widget:
        item[0].addWidget(item[1])

    for item in tab_add_layout:
        item[0].addLayout(item[1])
