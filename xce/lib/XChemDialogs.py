from PyQt4 import QtGui, QtCore

import XChemDB


class select_columns_to_show(QtGui.QDialog):
    def __init__(self, data_source_file, parent=None):
        super(select_columns_to_show, self).__init__(parent)
        self.columns_in_data_source = XChemDB.data_source(
            data_source_file
        ).return_column_list()

        self.column_dict = {}

        layout = QtGui.QVBoxLayout(self)
        number_of_entries = len(self.columns_in_data_source)
        columns_shown_in_dialog_column = 25
        grid = QtGui.QGridLayout()
        x = 0
        y = 0
        columns_to_ignore = ["Sample ID", "ID"]
        for entries_added in range(number_of_entries):
            if not self.columns_in_data_source[entries_added][1] in columns_to_ignore:
                data_source_column = QtGui.QCheckBox(
                    self.columns_in_data_source[entries_added][1]
                )
                self.column_dict[entries_added] = data_source_column
                #            data_source_column.toggle()
                grid.addWidget(data_source_column, y, x)
                y += 1
            if y == columns_shown_in_dialog_column:
                y = 0
                x += 1
        layout.addLayout(grid)

        # OK and Cancel buttons
        buttons = QtGui.QDialogButtonBox(
            QtGui.QDialogButtonBox.Ok | QtGui.QDialogButtonBox.Cancel,
            QtCore.Qt.Horizontal,
            self,
        )
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)

    def data_source_column_active(self):
        columns_to_show = ["Sample ID"]
        for key in self.column_dict:
            if self.column_dict[key].isChecked():
                columns_to_show.append(self.columns_in_data_source[key][1])
        return columns_to_show

    # static method to create the dialog and return (date, time, accepted)
    @staticmethod
    def return_selected_columns(parent=None):
        dialog = select_columns_to_show(parent)
        result = dialog.exec_()
        columns = dialog.data_source_column_active()
        return (columns, result == QtGui.QDialog.Accepted)
