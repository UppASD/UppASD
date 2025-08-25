# -*- coding: utf-8 -*-

################################################################################
## Form generated from reading UI file 'DMfile_Creator.ui'
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
    QHBoxLayout, QHeaderView, QLabel, QPushButton,
    QScrollArea, QSizePolicy, QSpacerItem, QSpinBox,
    QTableWidget, QTableWidgetItem, QWidget)

class Ui_InpDMfileCreate(object):
    def setupUi(self, InpDMfileCreate):
        if not InpDMfileCreate.objectName():
            InpDMfileCreate.setObjectName(u"InpDMfileCreate")
        InpDMfileCreate.resize(885, 500)
        self.gridLayout_3 = QGridLayout(InpDMfileCreate)
        self.gridLayout_3.setObjectName(u"gridLayout_3")
        self.InpDMfileScrollArea = QScrollArea(InpDMfileCreate)
        self.InpDMfileScrollArea.setObjectName(u"InpDMfileScrollArea")
        self.InpDMfileScrollArea.setWidgetResizable(True)
        self.scrollAreaWidgetContents = QWidget()
        self.scrollAreaWidgetContents.setObjectName(u"scrollAreaWidgetContents")
        self.scrollAreaWidgetContents.setGeometry(QRect(0, 0, 860, 500))
        self.gridLayout_2 = QGridLayout(self.scrollAreaWidgetContents)
        self.gridLayout_2.setObjectName(u"gridLayout_2")
        self.InDMfileButtonBox = QGroupBox(self.scrollAreaWidgetContents)
        self.InDMfileButtonBox.setObjectName(u"InDMfileButtonBox")
        self.horizontalLayout = QHBoxLayout(self.InDMfileButtonBox)
        self.horizontalLayout.setObjectName(u"horizontalLayout")
        self.InpDMfileDone = QPushButton(self.InDMfileButtonBox)
        self.InpDMfileDone.setObjectName(u"InpDMfileDone")

        self.horizontalLayout.addWidget(self.InpDMfileDone)

        self.InpDMfileCancel = QPushButton(self.InDMfileButtonBox)
        self.InpDMfileCancel.setObjectName(u"InpDMfileCancel")

        self.horizontalLayout.addWidget(self.InpDMfileCancel)


        self.gridLayout_2.addWidget(self.InDMfileButtonBox, 3, 0, 1, 2)

        self.InDMfileBox = QGroupBox(self.scrollAreaWidgetContents)
        self.InDMfileBox.setObjectName(u"InDMfileBox")
        self.gridLayout = QGridLayout(self.InDMfileBox)
        self.gridLayout.setObjectName(u"gridLayout")
        self.InDMfileTable = QTableWidget(self.InDMfileBox)
        if (self.InDMfileTable.columnCount() < 8):
            self.InDMfileTable.setColumnCount(8)
        __qtablewidgetitem = QTableWidgetItem()
        self.InDMfileTable.setHorizontalHeaderItem(0, __qtablewidgetitem)
        __qtablewidgetitem1 = QTableWidgetItem()
        self.InDMfileTable.setHorizontalHeaderItem(1, __qtablewidgetitem1)
        __qtablewidgetitem2 = QTableWidgetItem()
        self.InDMfileTable.setHorizontalHeaderItem(2, __qtablewidgetitem2)
        __qtablewidgetitem3 = QTableWidgetItem()
        self.InDMfileTable.setHorizontalHeaderItem(3, __qtablewidgetitem3)
        __qtablewidgetitem4 = QTableWidgetItem()
        self.InDMfileTable.setHorizontalHeaderItem(4, __qtablewidgetitem4)
        __qtablewidgetitem5 = QTableWidgetItem()
        self.InDMfileTable.setHorizontalHeaderItem(5, __qtablewidgetitem5)
        __qtablewidgetitem6 = QTableWidgetItem()
        self.InDMfileTable.setHorizontalHeaderItem(6, __qtablewidgetitem6)
        __qtablewidgetitem7 = QTableWidgetItem()
        self.InDMfileTable.setHorizontalHeaderItem(7, __qtablewidgetitem7)
        if (self.InDMfileTable.rowCount() < 1):
            self.InDMfileTable.setRowCount(1)
        __qtablewidgetitem8 = QTableWidgetItem()
        self.InDMfileTable.setItem(0, 0, __qtablewidgetitem8)
        __qtablewidgetitem9 = QTableWidgetItem()
        self.InDMfileTable.setItem(0, 1, __qtablewidgetitem9)
        __qtablewidgetitem10 = QTableWidgetItem()
        self.InDMfileTable.setItem(0, 2, __qtablewidgetitem10)
        __qtablewidgetitem11 = QTableWidgetItem()
        self.InDMfileTable.setItem(0, 3, __qtablewidgetitem11)
        __qtablewidgetitem12 = QTableWidgetItem()
        self.InDMfileTable.setItem(0, 4, __qtablewidgetitem12)
        __qtablewidgetitem13 = QTableWidgetItem()
        self.InDMfileTable.setItem(0, 5, __qtablewidgetitem13)
        __qtablewidgetitem14 = QTableWidgetItem()
        self.InDMfileTable.setItem(0, 6, __qtablewidgetitem14)
        __qtablewidgetitem15 = QTableWidgetItem()
        self.InDMfileTable.setItem(0, 7, __qtablewidgetitem15)
        self.InDMfileTable.setObjectName(u"InDMfileTable")
        self.InDMfileTable.setLineWidth(2)
        self.InDMfileTable.setMidLineWidth(1)
        self.InDMfileTable.setRowCount(1)

        self.gridLayout.addWidget(self.InDMfileTable, 0, 0, 1, 2)

        self.InDMfileAddRow = QPushButton(self.InDMfileBox)
        self.InDMfileAddRow.setObjectName(u"InDMfileAddRow")

        self.gridLayout.addWidget(self.InDMfileAddRow, 1, 0, 1, 1)

        self.InDMfileDelRow = QPushButton(self.InDMfileBox)
        self.InDMfileDelRow.setObjectName(u"InDMfileDelRow")

        self.gridLayout.addWidget(self.InDMfileDelRow, 1, 1, 1, 1)

        self.InDMfileGenVectors = QPushButton(self.InDMfileBox)
        self.InDMfileGenVectors.setObjectName(u"InDMfileGenVectors")

        self.gridLayout.addWidget(self.InDMfileGenVectors, 2, 1, 1, 1)

        self.InDMfileCutoffLabel = QLabel(self.InDMfileBox)
        self.InDMfileCutoffLabel.setObjectName(u"InDMfileCutoffLabel")

        self.gridLayout.addWidget(self.InDMfileCutoffLabel, 2, 0, 1, 1)

        self.InDMfileNNCutoff = QSpinBox(self.InDMfileBox)
        self.InDMfileNNCutoff.setObjectName(u"InDMfileNNCutoff")
        self.InDMfileNNCutoff.setMaximum(5)

        self.gridLayout.addWidget(self.InDMfileNNCutoff, 3, 0, 1, 1)


        self.gridLayout_2.addWidget(self.InDMfileBox, 0, 0, 1, 2)

        self.verticalSpacer = QSpacerItem(20, 29, QSizePolicy.Policy.Minimum, QSizePolicy.Policy.Expanding)

        self.gridLayout_2.addItem(self.verticalSpacer, 2, 0, 1, 1)

        self.InpDMfileScrollArea.setWidget(self.scrollAreaWidgetContents)

        self.gridLayout_3.addWidget(self.InpDMfileScrollArea, 0, 0, 1, 1)


        self.retranslateUi(InpDMfileCreate)

        QMetaObject.connectSlotsByName(InpDMfileCreate)
    # setupUi

    def retranslateUi(self, InpDMfileCreate):
        InpDMfileCreate.setWindowTitle(QCoreApplication.translate("InpDMfileCreate", u"DMfile Creation", None))
        self.InDMfileButtonBox.setTitle("")
        self.InpDMfileDone.setText(QCoreApplication.translate("InpDMfileCreate", u"Done!", None))
        self.InpDMfileCancel.setText(QCoreApplication.translate("InpDMfileCreate", u"Cancel", None))
        self.InDMfileBox.setTitle("")
        ___qtablewidgetitem = self.InDMfileTable.horizontalHeaderItem(0)
        ___qtablewidgetitem.setText(QCoreApplication.translate("InpDMfileCreate", u"Site i", None));
        ___qtablewidgetitem1 = self.InDMfileTable.horizontalHeaderItem(1)
        ___qtablewidgetitem1.setText(QCoreApplication.translate("InpDMfileCreate", u"Site j", None));
        ___qtablewidgetitem2 = self.InDMfileTable.horizontalHeaderItem(2)
        ___qtablewidgetitem2.setText(QCoreApplication.translate("InpDMfileCreate", u"X Dir.", None));
        ___qtablewidgetitem3 = self.InDMfileTable.horizontalHeaderItem(3)
        ___qtablewidgetitem3.setText(QCoreApplication.translate("InpDMfileCreate", u"Y Dir.", None));
        ___qtablewidgetitem4 = self.InDMfileTable.horizontalHeaderItem(4)
        ___qtablewidgetitem4.setText(QCoreApplication.translate("InpDMfileCreate", u"Z Dir.", None));
        ___qtablewidgetitem5 = self.InDMfileTable.horizontalHeaderItem(5)
        ___qtablewidgetitem5.setText(QCoreApplication.translate("InpDMfileCreate", u"DM x", None));
        ___qtablewidgetitem6 = self.InDMfileTable.horizontalHeaderItem(6)
        ___qtablewidgetitem6.setText(QCoreApplication.translate("InpDMfileCreate", u"DM y", None));
        ___qtablewidgetitem7 = self.InDMfileTable.horizontalHeaderItem(7)
        ___qtablewidgetitem7.setText(QCoreApplication.translate("InpDMfileCreate", u"DM z", None));

        __sortingEnabled = self.InDMfileTable.isSortingEnabled()
        self.InDMfileTable.setSortingEnabled(False)
        ___qtablewidgetitem8 = self.InDMfileTable.item(0, 0)
        ___qtablewidgetitem8.setText(QCoreApplication.translate("InpDMfileCreate", u"1", None));
        ___qtablewidgetitem9 = self.InDMfileTable.item(0, 1)
        ___qtablewidgetitem9.setText(QCoreApplication.translate("InpDMfileCreate", u"1", None));
        ___qtablewidgetitem10 = self.InDMfileTable.item(0, 2)
        ___qtablewidgetitem10.setText(QCoreApplication.translate("InpDMfileCreate", u"1.0", None));
        ___qtablewidgetitem11 = self.InDMfileTable.item(0, 3)
        ___qtablewidgetitem11.setText(QCoreApplication.translate("InpDMfileCreate", u"0.0", None));
        ___qtablewidgetitem12 = self.InDMfileTable.item(0, 4)
        ___qtablewidgetitem12.setText(QCoreApplication.translate("InpDMfileCreate", u"0.0", None));
        ___qtablewidgetitem13 = self.InDMfileTable.item(0, 5)
        ___qtablewidgetitem13.setText(QCoreApplication.translate("InpDMfileCreate", u"1.0", None));
        ___qtablewidgetitem14 = self.InDMfileTable.item(0, 6)
        ___qtablewidgetitem14.setText(QCoreApplication.translate("InpDMfileCreate", u"0.0", None));
        ___qtablewidgetitem15 = self.InDMfileTable.item(0, 7)
        ___qtablewidgetitem15.setText(QCoreApplication.translate("InpDMfileCreate", u"0.0", None));
        self.InDMfileTable.setSortingEnabled(__sortingEnabled)

        self.InDMfileAddRow.setText(QCoreApplication.translate("InpDMfileCreate", u"Add Row", None))
        self.InDMfileDelRow.setText(QCoreApplication.translate("InpDMfileCreate", u"Delete Row", None))
        self.InDMfileGenVectors.setText(QCoreApplication.translate("InpDMfileCreate", u"Generate neighbor vectors", None))
        self.InDMfileCutoffLabel.setText(QCoreApplication.translate("InpDMfileCreate", u"Cutoff for neighbor search", None))
#if QT_CONFIG(tooltip)
        self.InDMfileNNCutoff.setToolTip(QCoreApplication.translate("InpDMfileCreate", u"Set cutoff of nearest neighbor search", None))
#endif // QT_CONFIG(tooltip)
#if QT_CONFIG(statustip)
        self.InDMfileNNCutoff.setStatusTip(QCoreApplication.translate("InpDMfileCreate", u"Set cutoff of nearest neighbor search", None))
#endif // QT_CONFIG(statustip)
#if QT_CONFIG(whatsthis)
        self.InDMfileNNCutoff.setWhatsThis(QCoreApplication.translate("InpDMfileCreate", u"Set cutoff of nearest neighbor search", None))
#endif // QT_CONFIG(whatsthis)
    # retranslateUi

