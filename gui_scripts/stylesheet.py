from PyQt4 import QtGui
import os


def set_stylesheet(xce_object):
    palette = QtGui.QPalette()

    palette.setColor(QtGui.QPalette.Background, QtGui.QColor("#ececec"))
    xce_object.setPalette(palette)

    icons_directory = os.path.join((os.getenv("XChemExplorer_DIR")), "icons")

    xce_object.setStyleSheet(
        """
    QApplication {
    padding: 0px;
    }
    QMenuBar {
    background-color: rgb(236, 236, 236);
    } 
    QMenuBar::item {
    background-color: rgb(236, 236, 236);
    }
    QMenu {
    background-color: rgb(236, 236, 236);
    border: 1px solid rgb(184, 192, 210);
    }
    QMenu::item::selected {
    background-color: rgb(28, 105, 241);
    color: rgb(255, 255, 255);
    }
    QComboBox {
    background-color: rgb(255, 255, 255);
    selection-background-color: rgb(28, 105, 241);
    border: 1px solid rgb(184, 192, 210);
    }
    QComboBox QAbstractItemView {
    background: rgb(255, 255, 255);
    border: 1px solid rgb(184, 192, 210);
    }
    QComboBox::down-arrow {
    image: url("""
        + icons_directory
        + """/drop-down.png);
    background-color: rgb(255, 255, 255);
    }
    QComboBox::drop-down {
    background-color: rgb(255, 255, 255);
    }
    QPushButton {
    background-color: rgb(214, 230, 244);
    border: 1px solid rgb(184, 192, 210);
    padding: 3px;
    }
    QFrame {
    background-color: rgb(236, 236, 236);
    }
    QTableWidget {
    background-color: rgb(255, 255, 255);
    }
    QHeaderView::section {
    background-color: rgb(236, 236, 236);
    border: 1px solid rgb(184, 192, 210);
    }
    QTabWidget::pane {
    border-top: 1px solid rgb(184, 192, 210);
    border-bottom: 1px solid rgb(184, 192, 210);
    }
    QTabBar::tab {
    background-color: rgb(197,197,197);
    border: 1px solid rgb(184, 192, 210);
    padding: 3px;
    }
    QTabBar::tab::selected {
    background-color: rgb(214, 230, 244);
    border: 1px solid rgb(184, 192, 210);
    padding: 3px;
    }
    QScrollBar {
    background:  rgb(236, 236, 236);
    }
    """
    )

    QtGui.qApp.setStyle("Cleanlooks")
