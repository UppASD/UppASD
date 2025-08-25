# -*- coding: utf-8 -*-

################################################################################
## Form generated from reading UI file 'Jfile_Creator.ui'
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
from PySide6.QtWidgets import (QApplication, QCheckBox, QDialog, QGridLayout,
    QGroupBox, QHBoxLayout, QHeaderView, QPushButton,
    QScrollArea, QSizePolicy, QSpacerItem, QSpinBox,
    QTableWidget, QTableWidgetItem, QWidget)

class Ui_InpJfileCreate(object):
    def setupUi(self, InpJfileCreate):
        if not InpJfileCreate.objectName():
            InpJfileCreate.setObjectName(u"InpJfileCreate")
        InpJfileCreate.resize(696, 500)
        self.gridLayout_3 = QGridLayout(InpJfileCreate)
        self.gridLayout_3.setObjectName(u"gridLayout_3")
        self.InpJfileScrollArea = QScrollArea(InpJfileCreate)
        self.InpJfileScrollArea.setObjectName(u"InpJfileScrollArea")
        self.InpJfileScrollArea.setWidgetResizable(True)
        self.scrollAreaWidgetContents = QWidget()
        self.scrollAreaWidgetContents.setObjectName(u"scrollAreaWidgetContents")
        self.scrollAreaWidgetContents.setGeometry(QRect(0, 0, 670, 500))
        self.gridLayout_2 = QGridLayout(self.scrollAreaWidgetContents)
        self.gridLayout_2.setObjectName(u"gridLayout_2")
        self.InJfileButtonBox = QGroupBox(self.scrollAreaWidgetContents)
        self.InJfileButtonBox.setObjectName(u"InJfileButtonBox")
        self.horizontalLayout = QHBoxLayout(self.InJfileButtonBox)
        self.horizontalLayout.setObjectName(u"horizontalLayout")
        self.InpJfileDone = QPushButton(self.InJfileButtonBox)
        self.InpJfileDone.setObjectName(u"InpJfileDone")

        self.horizontalLayout.addWidget(self.InpJfileDone)

        self.InpJfileCancel = QPushButton(self.InJfileButtonBox)
        self.InpJfileCancel.setObjectName(u"InpJfileCancel")

        self.horizontalLayout.addWidget(self.InpJfileCancel)


        self.gridLayout_2.addWidget(self.InJfileButtonBox, 3, 0, 1, 2)

        self.InJfileBox = QGroupBox(self.scrollAreaWidgetContents)
        self.InJfileBox.setObjectName(u"InJfileBox")
        self.gridLayout = QGridLayout(self.InJfileBox)
        self.gridLayout.setObjectName(u"gridLayout")
        self.InJfileTable = QTableWidget(self.InJfileBox)
        if (self.InJfileTable.columnCount() < 6):
            self.InJfileTable.setColumnCount(6)
        __qtablewidgetitem = QTableWidgetItem()
        self.InJfileTable.setHorizontalHeaderItem(0, __qtablewidgetitem)
        __qtablewidgetitem1 = QTableWidgetItem()
        self.InJfileTable.setHorizontalHeaderItem(1, __qtablewidgetitem1)
        __qtablewidgetitem2 = QTableWidgetItem()
        self.InJfileTable.setHorizontalHeaderItem(2, __qtablewidgetitem2)
        __qtablewidgetitem3 = QTableWidgetItem()
        self.InJfileTable.setHorizontalHeaderItem(3, __qtablewidgetitem3)
        __qtablewidgetitem4 = QTableWidgetItem()
        self.InJfileTable.setHorizontalHeaderItem(4, __qtablewidgetitem4)
        __qtablewidgetitem5 = QTableWidgetItem()
        self.InJfileTable.setHorizontalHeaderItem(5, __qtablewidgetitem5)
        if (self.InJfileTable.rowCount() < 1):
            self.InJfileTable.setRowCount(1)
        __qtablewidgetitem6 = QTableWidgetItem()
        self.InJfileTable.setItem(0, 0, __qtablewidgetitem6)
        __qtablewidgetitem7 = QTableWidgetItem()
        self.InJfileTable.setItem(0, 1, __qtablewidgetitem7)
        __qtablewidgetitem8 = QTableWidgetItem()
        self.InJfileTable.setItem(0, 2, __qtablewidgetitem8)
        __qtablewidgetitem9 = QTableWidgetItem()
        self.InJfileTable.setItem(0, 3, __qtablewidgetitem9)
        __qtablewidgetitem10 = QTableWidgetItem()
        self.InJfileTable.setItem(0, 4, __qtablewidgetitem10)
        __qtablewidgetitem11 = QTableWidgetItem()
        self.InJfileTable.setItem(0, 5, __qtablewidgetitem11)
        self.InJfileTable.setObjectName(u"InJfileTable")
        self.InJfileTable.setLineWidth(2)
        self.InJfileTable.setMidLineWidth(1)
        self.InJfileTable.setRowCount(1)

        self.gridLayout.addWidget(self.InJfileTable, 0, 0, 1, 2)

        self.InJfileAddRow = QPushButton(self.InJfileBox)
        self.InJfileAddRow.setObjectName(u"InJfileAddRow")

        self.gridLayout.addWidget(self.InJfileAddRow, 1, 0, 1, 1)

        self.InJfileDelRow = QPushButton(self.InJfileBox)
        self.InJfileDelRow.setObjectName(u"InJfileDelRow")

        self.gridLayout.addWidget(self.InJfileDelRow, 1, 1, 1, 1)

        self.InJfileGenVectors = QPushButton(self.InJfileBox)
        self.InJfileGenVectors.setObjectName(u"InJfileGenVectors")

        self.gridLayout.addWidget(self.InJfileGenVectors, 2, 1, 1, 1)

        self.InJfileSymCheck = QCheckBox(self.InJfileBox)
        self.InJfileSymCheck.setObjectName(u"InJfileSymCheck")

        self.gridLayout.addWidget(self.InJfileSymCheck, 2, 0, 1, 1)

        self.InJfileNNCutoff = QSpinBox(self.InJfileBox)
        self.InJfileNNCutoff.setObjectName(u"InJfileNNCutoff")
        self.InJfileNNCutoff.setMaximum(5)

        self.gridLayout.addWidget(self.InJfileNNCutoff, 3, 0, 1, 1)


        self.gridLayout_2.addWidget(self.InJfileBox, 0, 0, 1, 2)

        self.verticalSpacer = QSpacerItem(20, 29, QSizePolicy.Policy.Minimum, QSizePolicy.Policy.Expanding)

        self.gridLayout_2.addItem(self.verticalSpacer, 2, 0, 1, 1)

        self.InpJfileScrollArea.setWidget(self.scrollAreaWidgetContents)

        self.gridLayout_3.addWidget(self.InpJfileScrollArea, 0, 0, 1, 1)


        self.retranslateUi(InpJfileCreate)

        QMetaObject.connectSlotsByName(InpJfileCreate)
    # setupUi

    def retranslateUi(self, InpJfileCreate):
        InpJfileCreate.setWindowTitle(QCoreApplication.translate("InpJfileCreate", u"Jfile Creation", None))
        self.InJfileButtonBox.setTitle("")
        self.InpJfileDone.setText(QCoreApplication.translate("InpJfileCreate", u"Done!", None))
        self.InpJfileCancel.setText(QCoreApplication.translate("InpJfileCreate", u"Cancel", None))
        self.InJfileBox.setTitle("")
        ___qtablewidgetitem = self.InJfileTable.horizontalHeaderItem(0)
        ___qtablewidgetitem.setText(QCoreApplication.translate("InpJfileCreate", u"Site i", None));
        ___qtablewidgetitem1 = self.InJfileTable.horizontalHeaderItem(1)
        ___qtablewidgetitem1.setText(QCoreApplication.translate("InpJfileCreate", u"Site j", None));
        ___qtablewidgetitem2 = self.InJfileTable.horizontalHeaderItem(2)
        ___qtablewidgetitem2.setText(QCoreApplication.translate("InpJfileCreate", u"X Dir.", None));
        ___qtablewidgetitem3 = self.InJfileTable.horizontalHeaderItem(3)
        ___qtablewidgetitem3.setText(QCoreApplication.translate("InpJfileCreate", u"Y Dir.", None));
        ___qtablewidgetitem4 = self.InJfileTable.horizontalHeaderItem(4)
        ___qtablewidgetitem4.setText(QCoreApplication.translate("InpJfileCreate", u"Z Dir.", None));
        ___qtablewidgetitem5 = self.InJfileTable.horizontalHeaderItem(5)
        ___qtablewidgetitem5.setText(QCoreApplication.translate("InpJfileCreate", u"J Mag", None));

        __sortingEnabled = self.InJfileTable.isSortingEnabled()
        self.InJfileTable.setSortingEnabled(False)
        ___qtablewidgetitem6 = self.InJfileTable.item(0, 0)
        ___qtablewidgetitem6.setText(QCoreApplication.translate("InpJfileCreate", u"1", None));
        ___qtablewidgetitem7 = self.InJfileTable.item(0, 1)
        ___qtablewidgetitem7.setText(QCoreApplication.translate("InpJfileCreate", u"1", None));
        ___qtablewidgetitem8 = self.InJfileTable.item(0, 2)
        ___qtablewidgetitem8.setText(QCoreApplication.translate("InpJfileCreate", u"1.0", None));
        ___qtablewidgetitem9 = self.InJfileTable.item(0, 3)
        ___qtablewidgetitem9.setText(QCoreApplication.translate("InpJfileCreate", u"0.0", None));
        ___qtablewidgetitem10 = self.InJfileTable.item(0, 4)
        ___qtablewidgetitem10.setText(QCoreApplication.translate("InpJfileCreate", u"0.0", None));
        ___qtablewidgetitem11 = self.InJfileTable.item(0, 5)
        ___qtablewidgetitem11.setText(QCoreApplication.translate("InpJfileCreate", u"1.0", None));
        self.InJfileTable.setSortingEnabled(__sortingEnabled)

        self.InJfileAddRow.setText(QCoreApplication.translate("InpJfileCreate", u"Add Row", None))
        self.InJfileDelRow.setText(QCoreApplication.translate("InpJfileCreate", u"Delete Row", None))
        self.InJfileGenVectors.setText(QCoreApplication.translate("InpJfileCreate", u"Generate neighbor vectors", None))
        self.InJfileSymCheck.setText(QCoreApplication.translate("InpJfileCreate", u"Use symmetry", None))
#if QT_CONFIG(tooltip)
        self.InJfileNNCutoff.setToolTip(QCoreApplication.translate("InpJfileCreate", u"Set cutoff of nearest neighbor search", None))
#endif // QT_CONFIG(tooltip)
#if QT_CONFIG(statustip)
        self.InJfileNNCutoff.setStatusTip(QCoreApplication.translate("InpJfileCreate", u"Set cutoff of nearest neighbor search", None))
#endif // QT_CONFIG(statustip)
#if QT_CONFIG(whatsthis)
        self.InJfileNNCutoff.setWhatsThis(QCoreApplication.translate("InpJfileCreate", u"Set cutoff of nearest neighbor search", None))
#endif // QT_CONFIG(whatsthis)
    # retranslateUi

