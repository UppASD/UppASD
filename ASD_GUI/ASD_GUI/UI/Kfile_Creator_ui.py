# -*- coding: utf-8 -*-

################################################################################
## Form generated from reading UI file 'Kfile_Creator.ui'
##
## Created by: Qt User Interface Compiler version 6.8.0
##
## WARNING! All changes made in this file will be lost when recompiling UI file!
################################################################################

from PySide6.QtCore import (QCoreApplication, QDate, QDateTime, QLocale,
    QMetaObject, QObject, QPoint, QRect,
    QSize, QTime, QUrl, Qt)
from PySide6.QtGui import (QBrush, QColor, QConicalGradient, QCursor,
    QFont, QFontDatabase, QGradient, QIcon,
    QImage, QKeySequence, QLinearGradient, QPainter,
    QPalette, QPixmap, QRadialGradient, QTransform)
from PySide6.QtWidgets import (QApplication, QDialog, QGridLayout, QGroupBox,
    QHBoxLayout, QHeaderView, QPushButton, QScrollArea,
    QSizePolicy, QSpacerItem, QTableWidget, QTableWidgetItem,
    QWidget)

class Ui_InpKfileCreate(object):
    def setupUi(self, InpKfileCreate):
        if not InpKfileCreate.objectName():
            InpKfileCreate.setObjectName(u"InpKfileCreate")
        InpKfileCreate.resize(885, 500)
        self.gridLayout_3 = QGridLayout(InpKfileCreate)
        self.gridLayout_3.setObjectName(u"gridLayout_3")
        self.InpKfileScrollArea = QScrollArea(InpKfileCreate)
        self.InpKfileScrollArea.setObjectName(u"InpKfileScrollArea")
        self.InpKfileScrollArea.setWidgetResizable(True)
        self.scrollAreaWidgetContents = QWidget()
        self.scrollAreaWidgetContents.setObjectName(u"scrollAreaWidgetContents")
        self.scrollAreaWidgetContents.setGeometry(QRect(0, 0, 860, 500))
        self.gridLayout_2 = QGridLayout(self.scrollAreaWidgetContents)
        self.gridLayout_2.setObjectName(u"gridLayout_2")
        self.InKfileButtonBox = QGroupBox(self.scrollAreaWidgetContents)
        self.InKfileButtonBox.setObjectName(u"InKfileButtonBox")
        self.horizontalLayout = QHBoxLayout(self.InKfileButtonBox)
        self.horizontalLayout.setObjectName(u"horizontalLayout")
        self.InpKfileDone = QPushButton(self.InKfileButtonBox)
        self.InpKfileDone.setObjectName(u"InpKfileDone")

        self.horizontalLayout.addWidget(self.InpKfileDone)

        self.InpKfileCancel = QPushButton(self.InKfileButtonBox)
        self.InpKfileCancel.setObjectName(u"InpKfileCancel")

        self.horizontalLayout.addWidget(self.InpKfileCancel)


        self.gridLayout_2.addWidget(self.InKfileButtonBox, 3, 0, 1, 2)

        self.InKfileBox = QGroupBox(self.scrollAreaWidgetContents)
        self.InKfileBox.setObjectName(u"InKfileBox")
        self.gridLayout = QGridLayout(self.InKfileBox)
        self.gridLayout.setObjectName(u"gridLayout")
        self.InKfileTable = QTableWidget(self.InKfileBox)
        if (self.InKfileTable.columnCount() < 8):
            self.InKfileTable.setColumnCount(8)
        __qtablewidgetitem = QTableWidgetItem()
        self.InKfileTable.setHorizontalHeaderItem(0, __qtablewidgetitem)
        __qtablewidgetitem1 = QTableWidgetItem()
        self.InKfileTable.setHorizontalHeaderItem(1, __qtablewidgetitem1)
        __qtablewidgetitem2 = QTableWidgetItem()
        self.InKfileTable.setHorizontalHeaderItem(2, __qtablewidgetitem2)
        __qtablewidgetitem3 = QTableWidgetItem()
        self.InKfileTable.setHorizontalHeaderItem(3, __qtablewidgetitem3)
        __qtablewidgetitem4 = QTableWidgetItem()
        self.InKfileTable.setHorizontalHeaderItem(4, __qtablewidgetitem4)
        __qtablewidgetitem5 = QTableWidgetItem()
        self.InKfileTable.setHorizontalHeaderItem(5, __qtablewidgetitem5)
        __qtablewidgetitem6 = QTableWidgetItem()
        self.InKfileTable.setHorizontalHeaderItem(6, __qtablewidgetitem6)
        __qtablewidgetitem7 = QTableWidgetItem()
        self.InKfileTable.setHorizontalHeaderItem(7, __qtablewidgetitem7)
        if (self.InKfileTable.rowCount() < 1):
            self.InKfileTable.setRowCount(1)
        __qtablewidgetitem8 = QTableWidgetItem()
        self.InKfileTable.setItem(0, 0, __qtablewidgetitem8)
        __qtablewidgetitem9 = QTableWidgetItem()
        self.InKfileTable.setItem(0, 1, __qtablewidgetitem9)
        __qtablewidgetitem10 = QTableWidgetItem()
        self.InKfileTable.setItem(0, 2, __qtablewidgetitem10)
        __qtablewidgetitem11 = QTableWidgetItem()
        self.InKfileTable.setItem(0, 3, __qtablewidgetitem11)
        __qtablewidgetitem12 = QTableWidgetItem()
        self.InKfileTable.setItem(0, 4, __qtablewidgetitem12)
        __qtablewidgetitem13 = QTableWidgetItem()
        self.InKfileTable.setItem(0, 5, __qtablewidgetitem13)
        __qtablewidgetitem14 = QTableWidgetItem()
        self.InKfileTable.setItem(0, 6, __qtablewidgetitem14)
        __qtablewidgetitem15 = QTableWidgetItem()
        self.InKfileTable.setItem(0, 7, __qtablewidgetitem15)
        self.InKfileTable.setObjectName(u"InKfileTable")
        self.InKfileTable.setLineWidth(2)
        self.InKfileTable.setMidLineWidth(1)
        self.InKfileTable.setRowCount(1)

        self.gridLayout.addWidget(self.InKfileTable, 0, 0, 1, 2)

        self.InKfileAddRow = QPushButton(self.InKfileBox)
        self.InKfileAddRow.setObjectName(u"InKfileAddRow")

        self.gridLayout.addWidget(self.InKfileAddRow, 1, 0, 1, 1)

        self.InKfileDelRow = QPushButton(self.InKfileBox)
        self.InKfileDelRow.setObjectName(u"InKfileDelRow")

        self.gridLayout.addWidget(self.InKfileDelRow, 1, 1, 1, 1)


        self.gridLayout_2.addWidget(self.InKfileBox, 0, 0, 1, 2)

        self.verticalSpacer = QSpacerItem(20, 29, QSizePolicy.Policy.Minimum, QSizePolicy.Policy.Expanding)

        self.gridLayout_2.addItem(self.verticalSpacer, 2, 0, 1, 1)

        self.InpKfileScrollArea.setWidget(self.scrollAreaWidgetContents)

        self.gridLayout_3.addWidget(self.InpKfileScrollArea, 0, 0, 1, 1)


        self.retranslateUi(InpKfileCreate)

        QMetaObject.connectSlotsByName(InpKfileCreate)
    # setupUi

    def retranslateUi(self, InpKfileCreate):
        InpKfileCreate.setWindowTitle(QCoreApplication.translate("InpKfileCreate", u"Kfile Creation", None))
        self.InKfileButtonBox.setTitle("")
        self.InpKfileDone.setText(QCoreApplication.translate("InpKfileCreate", u"Done!", None))
        self.InpKfileCancel.setText(QCoreApplication.translate("InpKfileCreate", u"Cancel", None))
        self.InKfileBox.setTitle("")
        ___qtablewidgetitem = self.InKfileTable.horizontalHeaderItem(0)
        ___qtablewidgetitem.setText(QCoreApplication.translate("InpKfileCreate", u"Site i", None));
        ___qtablewidgetitem1 = self.InKfileTable.horizontalHeaderItem(1)
        ___qtablewidgetitem1.setText(QCoreApplication.translate("InpKfileCreate", u"Site j", None));
        ___qtablewidgetitem2 = self.InKfileTable.horizontalHeaderItem(2)
        ___qtablewidgetitem2.setText(QCoreApplication.translate("InpKfileCreate", u"X Dir.", None));
        ___qtablewidgetitem3 = self.InKfileTable.horizontalHeaderItem(3)
        ___qtablewidgetitem3.setText(QCoreApplication.translate("InpKfileCreate", u"Y Dir.", None));
        ___qtablewidgetitem4 = self.InKfileTable.horizontalHeaderItem(4)
        ___qtablewidgetitem4.setText(QCoreApplication.translate("InpKfileCreate", u"Z Dir.", None));
        ___qtablewidgetitem5 = self.InKfileTable.horizontalHeaderItem(5)
        ___qtablewidgetitem5.setText(QCoreApplication.translate("InpKfileCreate", u"K x", None));
        ___qtablewidgetitem6 = self.InKfileTable.horizontalHeaderItem(6)
        ___qtablewidgetitem6.setText(QCoreApplication.translate("InpKfileCreate", u"K y", None));
        ___qtablewidgetitem7 = self.InKfileTable.horizontalHeaderItem(7)
        ___qtablewidgetitem7.setText(QCoreApplication.translate("InpKfileCreate", u"K z", None));

        __sortingEnabled = self.InKfileTable.isSortingEnabled()
        self.InKfileTable.setSortingEnabled(False)
        ___qtablewidgetitem8 = self.InKfileTable.item(0, 0)
        ___qtablewidgetitem8.setText(QCoreApplication.translate("InpKfileCreate", u"1", None));
        ___qtablewidgetitem9 = self.InKfileTable.item(0, 1)
        ___qtablewidgetitem9.setText(QCoreApplication.translate("InpKfileCreate", u"1", None));
        ___qtablewidgetitem10 = self.InKfileTable.item(0, 2)
        ___qtablewidgetitem10.setText(QCoreApplication.translate("InpKfileCreate", u"1.0", None));
        ___qtablewidgetitem11 = self.InKfileTable.item(0, 3)
        ___qtablewidgetitem11.setText(QCoreApplication.translate("InpKfileCreate", u"0.0", None));
        ___qtablewidgetitem12 = self.InKfileTable.item(0, 4)
        ___qtablewidgetitem12.setText(QCoreApplication.translate("InpKfileCreate", u"0.0", None));
        ___qtablewidgetitem13 = self.InKfileTable.item(0, 5)
        ___qtablewidgetitem13.setText(QCoreApplication.translate("InpKfileCreate", u"1.0", None));
        ___qtablewidgetitem14 = self.InKfileTable.item(0, 6)
        ___qtablewidgetitem14.setText(QCoreApplication.translate("InpKfileCreate", u"0.0", None));
        ___qtablewidgetitem15 = self.InKfileTable.item(0, 7)
        ___qtablewidgetitem15.setText(QCoreApplication.translate("InpKfileCreate", u"0.0", None));
        self.InKfileTable.setSortingEnabled(__sortingEnabled)

        self.InKfileAddRow.setText(QCoreApplication.translate("InpKfileCreate", u"Add Row", None))
        self.InKfileDelRow.setText(QCoreApplication.translate("InpKfileCreate", u"Delete Row", None))
    # retranslateUi

