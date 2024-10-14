# -*- coding: utf-8 -*-

################################################################################
## Form generated from reading UI file 'PyNCBIminer_main10.ui'
##
## Created by: Qt User Interface Compiler version 5.15.2
##
## WARNING! All changes made in this file will be lost when recompiling UI file!
################################################################################

from PySide2.QtCore import *
from PySide2.QtGui import *
from PySide2.QtWidgets import *


class Ui_MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setupUi(self)

    def setupUi(self, MainWindow):
        if not MainWindow.objectName():
            MainWindow.setObjectName(u"MainWindow")
        MainWindow.resize(2500, 1284)
        font = QFont()
        font.setFamily(u"Arial")
        font.setPointSize(10)
        MainWindow.setFont(font)
        self.actionChange_font = QAction(MainWindow)
        self.actionChange_font.setObjectName(u"actionChange_font")
        self.actionInstall_MAFFT = QAction(MainWindow)
        self.actionInstall_MAFFT.setObjectName(u"actionInstall_MAFFT")
        self.actionInstall_trimAl = QAction(MainWindow)
        self.actionInstall_trimAl.setObjectName(u"actionInstall_trimAl")
        self.centralwidget = QWidget(MainWindow)
        self.centralwidget.setObjectName(u"centralwidget")
        font1 = QFont()
        font1.setFamily(u"Arial")
        font1.setPointSize(10)
        font1.setBold(False)
        font1.setWeight(50)
        self.centralwidget.setFont(font1)
        self.horizontalLayout = QHBoxLayout(self.centralwidget)
        self.horizontalLayout.setSpacing(6)
        self.horizontalLayout.setObjectName(u"horizontalLayout")
        self.horizontalLayout.setSizeConstraint(QLayout.SetNoConstraint)
        self.horizontalLayout.setContentsMargins(9, 9, 9, 9)
        self.tabWidget = QTabWidget(self.centralwidget)
        self.tabWidget.setObjectName(u"tabWidget")
        sizePolicy = QSizePolicy(QSizePolicy.Preferred, QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.tabWidget.sizePolicy().hasHeightForWidth())
        self.tabWidget.setSizePolicy(sizePolicy)
        self.tabWidget.setFont(font1)
        self.tab = QWidget()
        self.tab.setObjectName(u"tab")
        font2 = QFont()
        font2.setFamily(u"Century")
        font2.setPointSize(9)
        font2.setBold(False)
        font2.setWeight(50)
        self.tab.setFont(font2)
        self.verticalLayout_10 = QVBoxLayout(self.tab)
        self.verticalLayout_10.setObjectName(u"verticalLayout_10")
        self.groupBox_7 = QGroupBox(self.tab)
        self.groupBox_7.setObjectName(u"groupBox_7")
        self.groupBox_7.setFont(font1)
        self.verticalLayout = QVBoxLayout(self.groupBox_7)
        self.verticalLayout.setObjectName(u"verticalLayout")
        self.horizontalLayout_5 = QHBoxLayout()
        self.horizontalLayout_5.setObjectName(u"horizontalLayout_5")
        self.label_32 = QLabel(self.groupBox_7)
        self.label_32.setObjectName(u"label_32")

        self.horizontalLayout_5.addWidget(self.label_32)

        self.wd = QLineEdit(self.groupBox_7)
        self.wd.setObjectName(u"wd")
        font3 = QFont()
        font3.setFamily(u"Arial")
        font3.setPointSize(9)
        font3.setBold(False)
        font3.setWeight(50)
        self.wd.setFont(font3)

        self.horizontalLayout_5.addWidget(self.wd)

        self.view = QPushButton(self.groupBox_7)
        self.view.setObjectName(u"view")

        self.horizontalLayout_5.addWidget(self.view)


        self.verticalLayout.addLayout(self.horizontalLayout_5)


        self.verticalLayout_10.addWidget(self.groupBox_7)

        self.horizontalLayout_23 = QHBoxLayout()
        self.horizontalLayout_23.setObjectName(u"horizontalLayout_23")
        self.groupBox_8 = QGroupBox(self.tab)
        self.groupBox_8.setObjectName(u"groupBox_8")
        sizePolicy1 = QSizePolicy(QSizePolicy.Preferred, QSizePolicy.Minimum)
        sizePolicy1.setHorizontalStretch(0)
        sizePolicy1.setVerticalStretch(0)
        sizePolicy1.setHeightForWidth(self.groupBox_8.sizePolicy().hasHeightForWidth())
        self.groupBox_8.setSizePolicy(sizePolicy1)
        self.groupBox_8.setFont(font1)
        self.verticalLayout_2 = QVBoxLayout(self.groupBox_8)
        self.verticalLayout_2.setObjectName(u"verticalLayout_2")
        self.label_24 = QLabel(self.groupBox_8)
        self.label_24.setObjectName(u"label_24")

        self.verticalLayout_2.addWidget(self.label_24)

        self.taxonomy = QPlainTextEdit(self.groupBox_8)
        self.taxonomy.setObjectName(u"taxonomy")
        sizePolicy2 = QSizePolicy(QSizePolicy.Preferred, QSizePolicy.Expanding)
        sizePolicy2.setHorizontalStretch(0)
        sizePolicy2.setVerticalStretch(0)
        sizePolicy2.setHeightForWidth(self.taxonomy.sizePolicy().hasHeightForWidth())
        self.taxonomy.setSizePolicy(sizePolicy2)
        self.taxonomy.setFont(font3)

        self.verticalLayout_2.addWidget(self.taxonomy)

        self.horizontalLayout_4 = QHBoxLayout()
        self.horizontalLayout_4.setObjectName(u"horizontalLayout_4")
        self.label_4 = QLabel(self.groupBox_8)
        self.label_4.setObjectName(u"label_4")

        self.horizontalLayout_4.addWidget(self.label_4)

        self.target_region = QComboBox(self.groupBox_8)
        self.target_region.addItem("")
        self.target_region.addItem("")
        self.target_region.addItem("")
        self.target_region.addItem("")
        self.target_region.addItem("")
        self.target_region.addItem("")
        self.target_region.addItem("")
        self.target_region.addItem("")
        self.target_region.setObjectName(u"target_region")
        font4 = QFont()
        font4.setFamily(u"Times New Roman")
        font4.setPointSize(9)
        font4.setBold(False)
        font4.setWeight(50)
        self.target_region.setFont(font4)
        self.target_region.setEditable(True)
        self.target_region.setMaxVisibleItems(5)
        self.target_region.setMaxCount(2147483645)

        self.horizontalLayout_4.addWidget(self.target_region)


        self.verticalLayout_2.addLayout(self.horizontalLayout_4)

        self.horizontalLayout_2 = QHBoxLayout()
        self.horizontalLayout_2.setObjectName(u"horizontalLayout_2")
        self.set_target_region = QPushButton(self.groupBox_8)
        self.set_target_region.setObjectName(u"set_target_region")

        self.horizontalLayout_2.addWidget(self.set_target_region)

        self.save_settings = QPushButton(self.groupBox_8)
        self.save_settings.setObjectName(u"save_settings")

        self.horizontalLayout_2.addWidget(self.save_settings)


        self.verticalLayout_2.addLayout(self.horizontalLayout_2)

        self.label_25 = QLabel(self.groupBox_8)
        self.label_25.setObjectName(u"label_25")

        self.verticalLayout_2.addWidget(self.label_25)

        self.entrez_qualifier = QPlainTextEdit(self.groupBox_8)
        self.entrez_qualifier.setObjectName(u"entrez_qualifier")
        self.entrez_qualifier.setFont(font3)

        self.verticalLayout_2.addWidget(self.entrez_qualifier)

        self.horizontalLayout_33 = QHBoxLayout()
        self.horizontalLayout_33.setObjectName(u"horizontalLayout_33")
        self.label_31 = QLabel(self.groupBox_8)
        self.label_31.setObjectName(u"label_31")

        self.horizontalLayout_33.addWidget(self.label_31)

        self.horizontalLayout_11 = QHBoxLayout()
        self.horizontalLayout_11.setObjectName(u"horizontalLayout_11")
        self.label_43 = QLabel(self.groupBox_8)
        self.label_43.setObjectName(u"label_43")

        self.horizontalLayout_11.addWidget(self.label_43)

        self.date_from = QLineEdit(self.groupBox_8)
        self.date_from.setObjectName(u"date_from")

        self.horizontalLayout_11.addWidget(self.date_from)

        self.label_42 = QLabel(self.groupBox_8)
        self.label_42.setObjectName(u"label_42")

        self.horizontalLayout_11.addWidget(self.label_42)

        self.date_to = QLineEdit(self.groupBox_8)
        self.date_to.setObjectName(u"date_to")

        self.horizontalLayout_11.addWidget(self.date_to)


        self.horizontalLayout_33.addLayout(self.horizontalLayout_11)


        self.verticalLayout_2.addLayout(self.horizontalLayout_33)

        self.horizontalLayout_10 = QHBoxLayout()
        self.horizontalLayout_10.setObjectName(u"horizontalLayout_10")
        self.label_27 = QLabel(self.groupBox_8)
        self.label_27.setObjectName(u"label_27")

        self.horizontalLayout_10.addWidget(self.label_27)

        self.entrez_email = QLineEdit(self.groupBox_8)
        self.entrez_email.setObjectName(u"entrez_email")
        self.entrez_email.setFont(font3)

        self.horizontalLayout_10.addWidget(self.entrez_email)


        self.verticalLayout_2.addLayout(self.horizontalLayout_10)

        self.horizontalLayout_32 = QHBoxLayout()
        self.horizontalLayout_32.setObjectName(u"horizontalLayout_32")
        self.marker_summary = QCheckBox(self.groupBox_8)
        self.marker_summary.setObjectName(u"marker_summary")

        self.horizontalLayout_32.addWidget(self.marker_summary)

        self.esearch = QPushButton(self.groupBox_8)
        self.esearch.setObjectName(u"esearch")

        self.horizontalLayout_32.addWidget(self.esearch)


        self.verticalLayout_2.addLayout(self.horizontalLayout_32)


        self.horizontalLayout_23.addWidget(self.groupBox_8)

        self.groupBox_5 = QGroupBox(self.tab)
        self.groupBox_5.setObjectName(u"groupBox_5")
        sizePolicy3 = QSizePolicy(QSizePolicy.Preferred, QSizePolicy.Maximum)
        sizePolicy3.setHorizontalStretch(0)
        sizePolicy3.setVerticalStretch(0)
        sizePolicy3.setHeightForWidth(self.groupBox_5.sizePolicy().hasHeightForWidth())
        self.groupBox_5.setSizePolicy(sizePolicy3)
        self.groupBox_5.setFont(font1)
        self.verticalLayout_3 = QVBoxLayout(self.groupBox_5)
        self.verticalLayout_3.setObjectName(u"verticalLayout_3")
        self.formLayout_5 = QFormLayout()
        self.formLayout_5.setObjectName(u"formLayout_5")
        self.label_33 = QLabel(self.groupBox_5)
        self.label_33.setObjectName(u"label_33")

        self.formLayout_5.setWidget(0, QFormLayout.LabelRole, self.label_33)

        self.initial_queries = QPlainTextEdit(self.groupBox_5)
        self.initial_queries.setObjectName(u"initial_queries")
        self.initial_queries.setFont(font3)

        self.formLayout_5.setWidget(0, QFormLayout.FieldRole, self.initial_queries)

        self.label_30 = QLabel(self.groupBox_5)
        self.label_30.setObjectName(u"label_30")

        self.formLayout_5.setWidget(1, QFormLayout.LabelRole, self.label_30)

        self.key_annotations = QPlainTextEdit(self.groupBox_5)
        self.key_annotations.setObjectName(u"key_annotations")
        self.key_annotations.setFont(font3)

        self.formLayout_5.setWidget(1, QFormLayout.FieldRole, self.key_annotations)

        self.label_29 = QLabel(self.groupBox_5)
        self.label_29.setObjectName(u"label_29")

        self.formLayout_5.setWidget(2, QFormLayout.LabelRole, self.label_29)

        self.exclude_sources = QPlainTextEdit(self.groupBox_5)
        self.exclude_sources.setObjectName(u"exclude_sources")
        self.exclude_sources.setFont(font3)

        self.formLayout_5.setWidget(2, QFormLayout.FieldRole, self.exclude_sources)

        self.label_18 = QLabel(self.groupBox_5)
        self.label_18.setObjectName(u"label_18")

        self.formLayout_5.setWidget(3, QFormLayout.LabelRole, self.label_18)

        self.max_length = QLineEdit(self.groupBox_5)
        self.max_length.setObjectName(u"max_length")
        self.max_length.setFont(font3)

        self.formLayout_5.setWidget(3, QFormLayout.FieldRole, self.max_length)

        self.label_7 = QLabel(self.groupBox_5)
        self.label_7.setObjectName(u"label_7")

        self.formLayout_5.setWidget(6, QFormLayout.LabelRole, self.label_7)

        self.expect_value = QLineEdit(self.groupBox_5)
        self.expect_value.setObjectName(u"expect_value")
        self.expect_value.setFont(font3)

        self.formLayout_5.setWidget(6, QFormLayout.FieldRole, self.expect_value)

        self.label_15 = QLabel(self.groupBox_5)
        self.label_15.setObjectName(u"label_15")

        self.formLayout_5.setWidget(7, QFormLayout.LabelRole, self.label_15)

        self.word_size = QLineEdit(self.groupBox_5)
        self.word_size.setObjectName(u"word_size")
        self.word_size.setFont(font3)

        self.formLayout_5.setWidget(7, QFormLayout.FieldRole, self.word_size)

        self.label_21 = QLabel(self.groupBox_5)
        self.label_21.setObjectName(u"label_21")

        self.formLayout_5.setWidget(8, QFormLayout.LabelRole, self.label_21)

        self.gap_costs = QLineEdit(self.groupBox_5)
        self.gap_costs.setObjectName(u"gap_costs")
        self.gap_costs.setFont(font3)

        self.formLayout_5.setWidget(8, QFormLayout.FieldRole, self.gap_costs)

        self.label_23 = QLabel(self.groupBox_5)
        self.label_23.setObjectName(u"label_23")

        self.formLayout_5.setWidget(9, QFormLayout.LabelRole, self.label_23)

        self.nucl_reward = QLineEdit(self.groupBox_5)
        self.nucl_reward.setObjectName(u"nucl_reward")
        self.nucl_reward.setFont(font3)

        self.formLayout_5.setWidget(9, QFormLayout.FieldRole, self.nucl_reward)

        self.label_22 = QLabel(self.groupBox_5)
        self.label_22.setObjectName(u"label_22")

        self.formLayout_5.setWidget(10, QFormLayout.LabelRole, self.label_22)

        self.nucl_penalty = QLineEdit(self.groupBox_5)
        self.nucl_penalty.setObjectName(u"nucl_penalty")
        self.nucl_penalty.setFont(font3)

        self.formLayout_5.setWidget(10, QFormLayout.FieldRole, self.nucl_penalty)


        self.verticalLayout_3.addLayout(self.formLayout_5)


        self.horizontalLayout_23.addWidget(self.groupBox_5)

        self.horizontalLayout_23.setStretch(0, 1)
        self.horizontalLayout_23.setStretch(1, 1)

        self.verticalLayout_10.addLayout(self.horizontalLayout_23)

        self.horizontalLayout_12 = QHBoxLayout()
        self.horizontalLayout_12.setObjectName(u"horizontalLayout_12")
        self.submit_new_blast = QPushButton(self.tab)
        self.submit_new_blast.setObjectName(u"submit_new_blast")
        self.submit_new_blast.setFont(font1)

        self.horizontalLayout_12.addWidget(self.submit_new_blast)

        self.load_previous_job = QPushButton(self.tab)
        self.load_previous_job.setObjectName(u"load_previous_job")
        self.load_previous_job.setFont(font1)

        self.horizontalLayout_12.addWidget(self.load_previous_job)


        self.verticalLayout_10.addLayout(self.horizontalLayout_12)

        self.tabWidget.addTab(self.tab, "")
        self.tab_2 = QWidget()
        self.tab_2.setObjectName(u"tab_2")
        self.verticalLayout_5 = QVBoxLayout(self.tab_2)
        self.verticalLayout_5.setObjectName(u"verticalLayout_5")
        self.groupBox_2 = QGroupBox(self.tab_2)
        self.groupBox_2.setObjectName(u"groupBox_2")
        self.verticalLayout_9 = QVBoxLayout(self.groupBox_2)
        self.verticalLayout_9.setObjectName(u"verticalLayout_9")
        self.control_extension = QCheckBox(self.groupBox_2)
        self.buttonGroup = QButtonGroup(MainWindow)
        self.buttonGroup.setObjectName(u"buttonGroup")
        self.buttonGroup.addButton(self.control_extension)
        self.control_extension.setObjectName(u"control_extension")

        self.verticalLayout_9.addWidget(self.control_extension)

        self.horizontalLayout_31 = QHBoxLayout()
        self.horizontalLayout_31.setObjectName(u"horizontalLayout_31")
        self.reduce_dataset = QCheckBox(self.groupBox_2)
        self.buttonGroup.addButton(self.reduce_dataset)
        self.reduce_dataset.setObjectName(u"reduce_dataset")

        self.horizontalLayout_31.addWidget(self.reduce_dataset)

        self.horizontalSpacer_6 = QSpacerItem(40, 20, QSizePolicy.Expanding, QSizePolicy.Minimum)

        self.horizontalLayout_31.addItem(self.horizontalSpacer_6)

        self.horizontalLayout_30 = QHBoxLayout()
        self.horizontalLayout_30.setObjectName(u"horizontalLayout_30")
        self.label_9 = QLabel(self.groupBox_2)
        self.label_9.setObjectName(u"label_9")

        self.horizontalLayout_30.addWidget(self.label_9)

        self.consensus_value = QComboBox(self.groupBox_2)
        self.consensus_value.addItem("")
        self.consensus_value.addItem("")
        self.consensus_value.setObjectName(u"consensus_value")
        sizePolicy4 = QSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        sizePolicy4.setHorizontalStretch(0)
        sizePolicy4.setVerticalStretch(0)
        sizePolicy4.setHeightForWidth(self.consensus_value.sizePolicy().hasHeightForWidth())
        self.consensus_value.setSizePolicy(sizePolicy4)
        self.consensus_value.setFont(font3)
        self.consensus_value.setEditable(True)
        self.consensus_value.setMaxVisibleItems(5)
        self.consensus_value.setMaxCount(2147483645)

        self.horizontalLayout_30.addWidget(self.consensus_value)


        self.horizontalLayout_31.addLayout(self.horizontalLayout_30)

        self.horizontalSpacer_5 = QSpacerItem(40, 20, QSizePolicy.Expanding, QSizePolicy.Minimum)

        self.horizontalLayout_31.addItem(self.horizontalSpacer_5)

        self.horizontalLayout_29 = QHBoxLayout()
        self.horizontalLayout_29.setObjectName(u"horizontalLayout_29")
        self.label_5 = QLabel(self.groupBox_2)
        self.label_5.setObjectName(u"label_5")

        self.horizontalLayout_29.addWidget(self.label_5)

        self.len_threshold = QLineEdit(self.groupBox_2)
        self.len_threshold.setObjectName(u"len_threshold")

        self.horizontalLayout_29.addWidget(self.len_threshold)


        self.horizontalLayout_31.addLayout(self.horizontalLayout_29)

        self.horizontalLayout_31.setStretch(0, 2)
        self.horizontalLayout_31.setStretch(1, 1)
        self.horizontalLayout_31.setStretch(2, 2)
        self.horizontalLayout_31.setStretch(3, 1)
        self.horizontalLayout_31.setStretch(4, 2)

        self.verticalLayout_9.addLayout(self.horizontalLayout_31)

        self.horizontalLayout_21 = QHBoxLayout()
        self.horizontalLayout_21.setObjectName(u"horizontalLayout_21")
        self.verticalLayout_4 = QVBoxLayout()
        self.verticalLayout_4.setObjectName(u"verticalLayout_4")
        self.horizontalLayout_13 = QHBoxLayout()
        self.horizontalLayout_13.setObjectName(u"horizontalLayout_13")
        self.label_8 = QLabel(self.groupBox_2)
        self.label_8.setObjectName(u"label_8")

        self.horizontalLayout_13.addWidget(self.label_8)

        self.in_path1 = QLineEdit(self.groupBox_2)
        self.in_path1.setObjectName(u"in_path1")

        self.horizontalLayout_13.addWidget(self.in_path1)

        self.view_in_path1 = QPushButton(self.groupBox_2)
        self.view_in_path1.setObjectName(u"view_in_path1")

        self.horizontalLayout_13.addWidget(self.view_in_path1)

        self.horizontalLayout_13.setStretch(0, 1)
        self.horizontalLayout_13.setStretch(1, 7)
        self.horizontalLayout_13.setStretch(2, 1)

        self.verticalLayout_4.addLayout(self.horizontalLayout_13)

        self.horizontalLayout_14 = QHBoxLayout()
        self.horizontalLayout_14.setObjectName(u"horizontalLayout_14")
        self.label_14 = QLabel(self.groupBox_2)
        self.label_14.setObjectName(u"label_14")

        self.horizontalLayout_14.addWidget(self.label_14)

        self.out_path1 = QLineEdit(self.groupBox_2)
        self.out_path1.setObjectName(u"out_path1")

        self.horizontalLayout_14.addWidget(self.out_path1)

        self.view_out_path1 = QPushButton(self.groupBox_2)
        self.view_out_path1.setObjectName(u"view_out_path1")

        self.horizontalLayout_14.addWidget(self.view_out_path1)

        self.horizontalLayout_14.setStretch(0, 1)
        self.horizontalLayout_14.setStretch(1, 7)
        self.horizontalLayout_14.setStretch(2, 1)

        self.verticalLayout_4.addLayout(self.horizontalLayout_14)


        self.horizontalLayout_21.addLayout(self.verticalLayout_4)

        self.run_filtering = QPushButton(self.groupBox_2)
        self.run_filtering.setObjectName(u"run_filtering")
        sizePolicy5 = QSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        sizePolicy5.setHorizontalStretch(0)
        sizePolicy5.setVerticalStretch(0)
        sizePolicy5.setHeightForWidth(self.run_filtering.sizePolicy().hasHeightForWidth())
        self.run_filtering.setSizePolicy(sizePolicy5)
        self.run_filtering.setAutoDefault(False)
        self.run_filtering.setFlat(False)

        self.horizontalLayout_21.addWidget(self.run_filtering)

        self.horizontalLayout_21.setStretch(0, 4)
        self.horizontalLayout_21.setStretch(1, 1)

        self.verticalLayout_9.addLayout(self.horizontalLayout_21)


        self.verticalLayout_5.addWidget(self.groupBox_2)

        self.groupBox = QGroupBox(self.tab_2)
        self.groupBox.setObjectName(u"groupBox")
        self.verticalLayout_11 = QVBoxLayout(self.groupBox)
        self.verticalLayout_11.setObjectName(u"verticalLayout_11")
        self.horizontalLayout_43 = QHBoxLayout()
        self.horizontalLayout_43.setObjectName(u"horizontalLayout_43")
        self.radioButton_3 = QRadioButton(self.groupBox)
        self.buttonGroup_2 = QButtonGroup(MainWindow)
        self.buttonGroup_2.setObjectName(u"buttonGroup_2")
        self.buttonGroup_2.addButton(self.radioButton_3)
        self.radioButton_3.setObjectName(u"radioButton_3")
        self.radioButton_3.setChecked(True)

        self.horizontalLayout_43.addWidget(self.radioButton_3)

        self.radioButton_4 = QRadioButton(self.groupBox)
        self.buttonGroup_2.addButton(self.radioButton_4)
        self.radioButton_4.setObjectName(u"radioButton_4")

        self.horizontalLayout_43.addWidget(self.radioButton_4)


        self.verticalLayout_11.addLayout(self.horizontalLayout_43)

        self.horizontalLayout_26 = QHBoxLayout()
        self.horizontalLayout_26.setObjectName(u"horizontalLayout_26")
        self.verticalLayout_6 = QVBoxLayout()
        self.verticalLayout_6.setObjectName(u"verticalLayout_6")
        self.horizontalLayout_8 = QHBoxLayout()
        self.horizontalLayout_8.setObjectName(u"horizontalLayout_8")
        self.label_13 = QLabel(self.groupBox)
        self.label_13.setObjectName(u"label_13")

        self.horizontalLayout_8.addWidget(self.label_13)

        self.in_path2 = QLineEdit(self.groupBox)
        self.in_path2.setObjectName(u"in_path2")

        self.horizontalLayout_8.addWidget(self.in_path2)

        self.view_in_path2 = QPushButton(self.groupBox)
        self.view_in_path2.setObjectName(u"view_in_path2")

        self.horizontalLayout_8.addWidget(self.view_in_path2)

        self.horizontalLayout_8.setStretch(0, 1)
        self.horizontalLayout_8.setStretch(1, 7)
        self.horizontalLayout_8.setStretch(2, 1)

        self.verticalLayout_6.addLayout(self.horizontalLayout_8)

        self.horizontalLayout_17 = QHBoxLayout()
        self.horizontalLayout_17.setObjectName(u"horizontalLayout_17")
        self.label_16 = QLabel(self.groupBox)
        self.label_16.setObjectName(u"label_16")

        self.horizontalLayout_17.addWidget(self.label_16)

        self.out_path2 = QLineEdit(self.groupBox)
        self.out_path2.setObjectName(u"out_path2")

        self.horizontalLayout_17.addWidget(self.out_path2)

        self.view_out_path2 = QPushButton(self.groupBox)
        self.view_out_path2.setObjectName(u"view_out_path2")

        self.horizontalLayout_17.addWidget(self.view_out_path2)

        self.horizontalLayout_17.setStretch(0, 1)
        self.horizontalLayout_17.setStretch(1, 7)
        self.horizontalLayout_17.setStretch(2, 1)

        self.verticalLayout_6.addLayout(self.horizontalLayout_17)

        self.horizontalLayout_25 = QHBoxLayout()
        self.horizontalLayout_25.setObjectName(u"horizontalLayout_25")
        self.horizontalLayout_9 = QHBoxLayout()
        self.horizontalLayout_9.setObjectName(u"horizontalLayout_9")
        self.label_2 = QLabel(self.groupBox)
        self.label_2.setObjectName(u"label_2")

        self.horizontalLayout_9.addWidget(self.label_2)

        self.ali_thr = QLineEdit(self.groupBox)
        self.ali_thr.setObjectName(u"ali_thr")
        sizePolicy4.setHeightForWidth(self.ali_thr.sizePolicy().hasHeightForWidth())
        self.ali_thr.setSizePolicy(sizePolicy4)

        self.horizontalLayout_9.addWidget(self.ali_thr)


        self.horizontalLayout_25.addLayout(self.horizontalLayout_9)

        self.horizontalLayout_16 = QHBoxLayout()
        self.horizontalLayout_16.setObjectName(u"horizontalLayout_16")
        self.label_3 = QLabel(self.groupBox)
        self.label_3.setObjectName(u"label_3")

        self.horizontalLayout_16.addWidget(self.label_3)

        self.ali_reo = QComboBox(self.groupBox)
        self.ali_reo.addItem("")
        self.ali_reo.addItem("")
        self.ali_reo.setObjectName(u"ali_reo")
        sizePolicy4.setHeightForWidth(self.ali_reo.sizePolicy().hasHeightForWidth())
        self.ali_reo.setSizePolicy(sizePolicy4)
        self.ali_reo.setFont(font3)
        self.ali_reo.setEditable(True)
        self.ali_reo.setMaxVisibleItems(5)
        self.ali_reo.setMaxCount(2147483645)

        self.horizontalLayout_16.addWidget(self.ali_reo)


        self.horizontalLayout_25.addLayout(self.horizontalLayout_16)


        self.verticalLayout_6.addLayout(self.horizontalLayout_25)

        self.horizontalLayout_20 = QHBoxLayout()
        self.horizontalLayout_20.setObjectName(u"horizontalLayout_20")
        self.label = QLabel(self.groupBox)
        self.label.setObjectName(u"label")

        self.horizontalLayout_20.addWidget(self.label)

        self.ali_alg = QComboBox(self.groupBox)
        self.ali_alg.addItem("")
        self.ali_alg.addItem("")
        self.ali_alg.setObjectName(u"ali_alg")
        sizePolicy4.setHeightForWidth(self.ali_alg.sizePolicy().hasHeightForWidth())
        self.ali_alg.setSizePolicy(sizePolicy4)
        self.ali_alg.setFont(font3)
        self.ali_alg.setEditable(True)
        self.ali_alg.setMaxVisibleItems(5)
        self.ali_alg.setMaxCount(2147483645)

        self.horizontalLayout_20.addWidget(self.ali_alg)


        self.verticalLayout_6.addLayout(self.horizontalLayout_20)


        self.horizontalLayout_26.addLayout(self.verticalLayout_6)

        self.run_alignment = QPushButton(self.groupBox)
        self.run_alignment.setObjectName(u"run_alignment")
        sizePolicy5.setHeightForWidth(self.run_alignment.sizePolicy().hasHeightForWidth())
        self.run_alignment.setSizePolicy(sizePolicy5)

        self.horizontalLayout_26.addWidget(self.run_alignment)

        self.horizontalLayout_26.setStretch(0, 4)
        self.horizontalLayout_26.setStretch(1, 1)

        self.verticalLayout_11.addLayout(self.horizontalLayout_26)


        self.verticalLayout_5.addWidget(self.groupBox)

        self.groupBox_4 = QGroupBox(self.tab_2)
        self.groupBox_4.setObjectName(u"groupBox_4")
        self.verticalLayout_12 = QVBoxLayout(self.groupBox_4)
        self.verticalLayout_12.setObjectName(u"verticalLayout_12")
        self.horizontalLayout_27 = QHBoxLayout()
        self.horizontalLayout_27.setObjectName(u"horizontalLayout_27")
        self.radioButton_5 = QRadioButton(self.groupBox_4)
        self.buttonGroup_3 = QButtonGroup(MainWindow)
        self.buttonGroup_3.setObjectName(u"buttonGroup_3")
        self.buttonGroup_3.addButton(self.radioButton_5)
        self.radioButton_5.setObjectName(u"radioButton_5")
        self.radioButton_5.setChecked(True)

        self.horizontalLayout_27.addWidget(self.radioButton_5)

        self.radioButton_6 = QRadioButton(self.groupBox_4)
        self.buttonGroup_3.addButton(self.radioButton_6)
        self.radioButton_6.setObjectName(u"radioButton_6")

        self.horizontalLayout_27.addWidget(self.radioButton_6)

        self.horizontalLayout_27.setStretch(0, 1)
        self.horizontalLayout_27.setStretch(1, 1)

        self.verticalLayout_12.addLayout(self.horizontalLayout_27)

        self.horizontalLayout_18 = QHBoxLayout()
        self.horizontalLayout_18.setObjectName(u"horizontalLayout_18")
        self.verticalLayout_7 = QVBoxLayout()
        self.verticalLayout_7.setObjectName(u"verticalLayout_7")
        self.horizontalLayout_22 = QHBoxLayout()
        self.horizontalLayout_22.setObjectName(u"horizontalLayout_22")
        self.label_17 = QLabel(self.groupBox_4)
        self.label_17.setObjectName(u"label_17")

        self.horizontalLayout_22.addWidget(self.label_17)

        self.in_path3 = QLineEdit(self.groupBox_4)
        self.in_path3.setObjectName(u"in_path3")

        self.horizontalLayout_22.addWidget(self.in_path3)

        self.view_in_path3 = QPushButton(self.groupBox_4)
        self.view_in_path3.setObjectName(u"view_in_path3")

        self.horizontalLayout_22.addWidget(self.view_in_path3)

        self.horizontalLayout_22.setStretch(0, 1)
        self.horizontalLayout_22.setStretch(1, 7)
        self.horizontalLayout_22.setStretch(2, 1)

        self.verticalLayout_7.addLayout(self.horizontalLayout_22)

        self.horizontalLayout_6 = QHBoxLayout()
        self.horizontalLayout_6.setObjectName(u"horizontalLayout_6")
        self.label_26 = QLabel(self.groupBox_4)
        self.label_26.setObjectName(u"label_26")

        self.horizontalLayout_6.addWidget(self.label_26)

        self.out_path3 = QLineEdit(self.groupBox_4)
        self.out_path3.setObjectName(u"out_path3")

        self.horizontalLayout_6.addWidget(self.out_path3)

        self.view_out_path3 = QPushButton(self.groupBox_4)
        self.view_out_path3.setObjectName(u"view_out_path3")

        self.horizontalLayout_6.addWidget(self.view_out_path3)

        self.horizontalLayout_6.setStretch(0, 1)
        self.horizontalLayout_6.setStretch(1, 7)
        self.horizontalLayout_6.setStretch(2, 1)

        self.verticalLayout_7.addLayout(self.horizontalLayout_6)

        self.horizontalLayout_19 = QHBoxLayout()
        self.horizontalLayout_19.setObjectName(u"horizontalLayout_19")
        self.horizontalLayout_15 = QHBoxLayout()
        self.horizontalLayout_15.setObjectName(u"horizontalLayout_15")
        self.label_6 = QLabel(self.groupBox_4)
        self.label_6.setObjectName(u"label_6")
        sizePolicy6 = QSizePolicy(QSizePolicy.Minimum, QSizePolicy.Preferred)
        sizePolicy6.setHorizontalStretch(0)
        sizePolicy6.setVerticalStretch(0)
        sizePolicy6.setHeightForWidth(self.label_6.sizePolicy().hasHeightForWidth())
        self.label_6.setSizePolicy(sizePolicy6)

        self.horizontalLayout_15.addWidget(self.label_6)

        self.tri_met = QComboBox(self.groupBox_4)
        self.tri_met.addItem("")
        self.tri_met.addItem("")
        self.tri_met.addItem("")
        self.tri_met.addItem("")
        self.tri_met.addItem("")
        self.tri_met.setObjectName(u"tri_met")
        sizePolicy4.setHeightForWidth(self.tri_met.sizePolicy().hasHeightForWidth())
        self.tri_met.setSizePolicy(sizePolicy4)
        self.tri_met.setFont(font3)
        self.tri_met.setEditable(True)
        self.tri_met.setMaxVisibleItems(5)
        self.tri_met.setMaxCount(2147483645)

        self.horizontalLayout_15.addWidget(self.tri_met)


        self.horizontalLayout_19.addLayout(self.horizontalLayout_15)

        self.horizontalLayout_19.setStretch(0, 3)

        self.verticalLayout_7.addLayout(self.horizontalLayout_19)

        self.horizontalLayout_57 = QHBoxLayout()
        self.horizontalLayout_57.setObjectName(u"horizontalLayout_57")
        self.horizontalLayout_52 = QHBoxLayout()
        self.horizontalLayout_52.setObjectName(u"horizontalLayout_52")
        self.label_48 = QLabel(self.groupBox_4)
        self.label_48.setObjectName(u"label_48")

        self.horizontalLayout_52.addWidget(self.label_48)

        self.tri_gt = QLineEdit(self.groupBox_4)
        self.tri_gt.setObjectName(u"tri_gt")

        self.horizontalLayout_52.addWidget(self.tri_gt)


        self.horizontalLayout_57.addLayout(self.horizontalLayout_52)

        self.horizontalLayout_53 = QHBoxLayout()
        self.horizontalLayout_53.setObjectName(u"horizontalLayout_53")
        self.label_49 = QLabel(self.groupBox_4)
        self.label_49.setObjectName(u"label_49")

        self.horizontalLayout_53.addWidget(self.label_49)

        self.tri_st = QLineEdit(self.groupBox_4)
        self.tri_st.setObjectName(u"tri_st")

        self.horizontalLayout_53.addWidget(self.tri_st)


        self.horizontalLayout_57.addLayout(self.horizontalLayout_53)

        self.horizontalLayout_57.setStretch(0, 1)
        self.horizontalLayout_57.setStretch(1, 1)

        self.verticalLayout_7.addLayout(self.horizontalLayout_57)

        self.horizontalLayout_61 = QHBoxLayout()
        self.horizontalLayout_61.setObjectName(u"horizontalLayout_61")
        self.horizontalLayout_54 = QHBoxLayout()
        self.horizontalLayout_54.setObjectName(u"horizontalLayout_54")
        self.label_50 = QLabel(self.groupBox_4)
        self.label_50.setObjectName(u"label_50")

        self.horizontalLayout_54.addWidget(self.label_50)

        self.tri_ct = QLineEdit(self.groupBox_4)
        self.tri_ct.setObjectName(u"tri_ct")

        self.horizontalLayout_54.addWidget(self.tri_ct)


        self.horizontalLayout_61.addLayout(self.horizontalLayout_54)

        self.horizontalLayout_55 = QHBoxLayout()
        self.horizontalLayout_55.setObjectName(u"horizontalLayout_55")
        self.label_51 = QLabel(self.groupBox_4)
        self.label_51.setObjectName(u"label_51")

        self.horizontalLayout_55.addWidget(self.label_51)

        self.tri_con = QLineEdit(self.groupBox_4)
        self.tri_con.setObjectName(u"tri_con")

        self.horizontalLayout_55.addWidget(self.tri_con)


        self.horizontalLayout_61.addLayout(self.horizontalLayout_55)

        self.horizontalLayout_61.setStretch(0, 1)
        self.horizontalLayout_61.setStretch(1, 1)

        self.verticalLayout_7.addLayout(self.horizontalLayout_61)


        self.horizontalLayout_18.addLayout(self.verticalLayout_7)

        self.run_trimming = QPushButton(self.groupBox_4)
        self.run_trimming.setObjectName(u"run_trimming")
        sizePolicy5.setHeightForWidth(self.run_trimming.sizePolicy().hasHeightForWidth())
        self.run_trimming.setSizePolicy(sizePolicy5)

        self.horizontalLayout_18.addWidget(self.run_trimming)

        self.horizontalLayout_18.setStretch(0, 4)
        self.horizontalLayout_18.setStretch(1, 1)

        self.verticalLayout_12.addLayout(self.horizontalLayout_18)


        self.verticalLayout_5.addWidget(self.groupBox_4)

        self.groupBox_6 = QGroupBox(self.tab_2)
        self.groupBox_6.setObjectName(u"groupBox_6")
        self.horizontalLayout_28 = QHBoxLayout(self.groupBox_6)
        self.horizontalLayout_28.setObjectName(u"horizontalLayout_28")
        self.verticalLayout_8 = QVBoxLayout()
        self.verticalLayout_8.setObjectName(u"verticalLayout_8")
        self.horizontalLayout_64 = QHBoxLayout()
        self.horizontalLayout_64.setObjectName(u"horizontalLayout_64")
        self.label_55 = QLabel(self.groupBox_6)
        self.label_55.setObjectName(u"label_55")

        self.horizontalLayout_64.addWidget(self.label_55)

        self.in_path4 = QLineEdit(self.groupBox_6)
        self.in_path4.setObjectName(u"in_path4")

        self.horizontalLayout_64.addWidget(self.in_path4)

        self.view_in_path4 = QPushButton(self.groupBox_6)
        self.view_in_path4.setObjectName(u"view_in_path4")

        self.horizontalLayout_64.addWidget(self.view_in_path4)

        self.horizontalLayout_64.setStretch(0, 1)
        self.horizontalLayout_64.setStretch(1, 7)
        self.horizontalLayout_64.setStretch(2, 1)

        self.verticalLayout_8.addLayout(self.horizontalLayout_64)

        self.horizontalLayout_24 = QHBoxLayout()
        self.horizontalLayout_24.setObjectName(u"horizontalLayout_24")
        self.label_39 = QLabel(self.groupBox_6)
        self.label_39.setObjectName(u"label_39")

        self.horizontalLayout_24.addWidget(self.label_39)

        self.out_path4 = QLineEdit(self.groupBox_6)
        self.out_path4.setObjectName(u"out_path4")

        self.horizontalLayout_24.addWidget(self.out_path4)

        self.view_out_path4 = QPushButton(self.groupBox_6)
        self.view_out_path4.setObjectName(u"view_out_path4")

        self.horizontalLayout_24.addWidget(self.view_out_path4)

        self.horizontalLayout_24.setStretch(0, 1)
        self.horizontalLayout_24.setStretch(1, 7)
        self.horizontalLayout_24.setStretch(2, 1)

        self.verticalLayout_8.addLayout(self.horizontalLayout_24)


        self.horizontalLayout_28.addLayout(self.verticalLayout_8)

        self.horizontalLayout_7 = QHBoxLayout()
        self.horizontalLayout_7.setObjectName(u"horizontalLayout_7")
        self.run_concatenation = QPushButton(self.groupBox_6)
        self.run_concatenation.setObjectName(u"run_concatenation")
        sizePolicy7 = QSizePolicy(QSizePolicy.Expanding, QSizePolicy.Ignored)
        sizePolicy7.setHorizontalStretch(0)
        sizePolicy7.setVerticalStretch(0)
        sizePolicy7.setHeightForWidth(self.run_concatenation.sizePolicy().hasHeightForWidth())
        self.run_concatenation.setSizePolicy(sizePolicy7)

        self.horizontalLayout_7.addWidget(self.run_concatenation)


        self.horizontalLayout_28.addLayout(self.horizontalLayout_7)

        self.horizontalLayout_28.setStretch(0, 4)
        self.horizontalLayout_28.setStretch(1, 1)

        self.verticalLayout_5.addWidget(self.groupBox_6)

        self.verticalLayout_5.setStretch(1, 2)
        self.verticalLayout_5.setStretch(2, 3)
        self.verticalLayout_5.setStretch(3, 1)
        self.tabWidget.addTab(self.tab_2, "")

        self.horizontalLayout.addWidget(self.tabWidget)

        self.groupBox_3 = QGroupBox(self.centralwidget)
        self.groupBox_3.setObjectName(u"groupBox_3")
        self.horizontalLayout_3 = QHBoxLayout(self.groupBox_3)
        self.horizontalLayout_3.setSpacing(4)
        self.horizontalLayout_3.setObjectName(u"horizontalLayout_3")
        self.horizontalLayout_3.setSizeConstraint(QLayout.SetMaximumSize)
        self.message_box = QPlainTextEdit(self.groupBox_3)
        self.message_box.setObjectName(u"message_box")
        sizePolicy2.setHeightForWidth(self.message_box.sizePolicy().hasHeightForWidth())
        self.message_box.setSizePolicy(sizePolicy2)
        self.message_box.setFont(font3)

        self.horizontalLayout_3.addWidget(self.message_box)

        self.horizontalLayout_3.setStretch(0, 1)

        self.horizontalLayout.addWidget(self.groupBox_3)

        self.horizontalLayout.setStretch(0, 2)
        self.horizontalLayout.setStretch(1, 1)
        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QMenuBar(MainWindow)
        self.menubar.setObjectName(u"menubar")
        self.menubar.setGeometry(QRect(0, 0, 1499, 21))
        self.menuTools = QMenu(self.menubar)
        self.menuTools.setObjectName(u"menuTools")
        self.menuTools.setFont(font)
        self.menuInstall_dependencies = QMenu(self.menuTools)
        self.menuInstall_dependencies.setObjectName(u"menuInstall_dependencies")
        self.menuInstall_dependencies.setFont(font)
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QStatusBar(MainWindow)
        self.statusbar.setObjectName(u"statusbar")
        MainWindow.setStatusBar(self.statusbar)

        self.menubar.addAction(self.menuTools.menuAction())
        self.menuTools.addAction(self.menuInstall_dependencies.menuAction())
        self.menuInstall_dependencies.addAction(self.actionInstall_MAFFT)
        self.menuInstall_dependencies.addAction(self.actionInstall_trimAl)

        self.retranslateUi(MainWindow)

        self.tabWidget.setCurrentIndex(0)
        self.run_filtering.setDefault(False)


        QMetaObject.connectSlotsByName(MainWindow)
    # setupUi

    def retranslateUi(self, MainWindow):
        MainWindow.setWindowTitle(QCoreApplication.translate("MainWindow", u"PyNCBIminer", None))
        self.actionChange_font.setText(QCoreApplication.translate("MainWindow", u"Change font", None))
        self.actionInstall_MAFFT.setText(QCoreApplication.translate("MainWindow", u"Install MAFFT", None))
        self.actionInstall_trimAl.setText(QCoreApplication.translate("MainWindow", u"Install Trimal", None))
        self.groupBox_7.setTitle(QCoreApplication.translate("MainWindow", u"Working Directory", None))
        self.label_32.setText(QCoreApplication.translate("MainWindow", u"working directory:", None))
        self.wd.setPlaceholderText(QCoreApplication.translate("MainWindow", u"absolute path of working directory", None))
        self.view.setText(QCoreApplication.translate("MainWindow", u"view", None))
        self.groupBox_8.setTitle(QCoreApplication.translate("MainWindow", u"Basic Settings", None))
        self.label_24.setText(QCoreApplication.translate("MainWindow", u"target groups:", None))
        self.taxonomy.setPlaceholderText(QCoreApplication.translate("MainWindow", u"one taxon per line", None))
        self.label_4.setText(QCoreApplication.translate("MainWindow", u"select or input target regeion:", None))
        self.target_region.setItemText(0, "")
        self.target_region.setItemText(1, QCoreApplication.translate("MainWindow", u"ITS", None))
        self.target_region.setItemText(2, QCoreApplication.translate("MainWindow", u"rbcL", None))
        self.target_region.setItemText(3, QCoreApplication.translate("MainWindow", u"matK", None))
        self.target_region.setItemText(4, QCoreApplication.translate("MainWindow", u"trnL-trnF", None))
        self.target_region.setItemText(5, QCoreApplication.translate("MainWindow", u"psbA-trnH", None))
        self.target_region.setItemText(6, QCoreApplication.translate("MainWindow", u"ndhF", None))
        self.target_region.setItemText(7, QCoreApplication.translate("MainWindow", u"rpoB", None))

        self.target_region.setCurrentText("")
        self.target_region.setPlaceholderText("")
        self.set_target_region.setText(QCoreApplication.translate("MainWindow", u"set target region", None))
        self.save_settings.setText(QCoreApplication.translate("MainWindow", u"save settings", None))
        self.label_25.setText(QCoreApplication.translate("MainWindow", u"entrez qualifier:", None))
        self.entrez_qualifier.setPlaceholderText(QCoreApplication.translate("MainWindow", u"constraint on BLAST search", None))
        self.label_31.setText(QCoreApplication.translate("MainWindow", u"publication date:  ", None))
        self.label_43.setText(QCoreApplication.translate("MainWindow", u"from", None))
        self.date_from.setPlaceholderText(QCoreApplication.translate("MainWindow", u"YYYY/MM/DD", None))
        self.label_42.setText(QCoreApplication.translate("MainWindow", u"to", None))
        self.date_to.setPlaceholderText(QCoreApplication.translate("MainWindow", u"YYYY/MM/DD", None))
        self.label_27.setText(QCoreApplication.translate("MainWindow", u"entrez email:", None))
        self.entrez_email.setPlaceholderText(QCoreApplication.translate("MainWindow", u"user's email", None))
        self.marker_summary.setText(QCoreApplication.translate("MainWindow", u"summarize widely used marker", None))
        self.esearch.setText(QCoreApplication.translate("MainWindow", u"entrez search", None))
        self.groupBox_5.setTitle(QCoreApplication.translate("MainWindow", u"Advanced Settings", None))
        self.label_33.setText(QCoreApplication.translate("MainWindow", u"initial queries:", None))
        self.initial_queries.setPlaceholderText(QCoreApplication.translate("MainWindow", u"paste sequences in fasta format here", None))
        self.label_30.setText(QCoreApplication.translate("MainWindow", u"key annotations:", None))
        self.key_annotations.setPlaceholderText(QCoreApplication.translate("MainWindow", u"one keyword per line", None))
        self.label_29.setText(QCoreApplication.translate("MainWindow", u"exclude sources:", None))
        self.exclude_sources.setPlaceholderText(QCoreApplication.translate("MainWindow", u"one keyword per line", None))
        self.label_18.setText(QCoreApplication.translate("MainWindow", u"max length:", None))
        self.max_length.setPlaceholderText(QCoreApplication.translate("MainWindow", u"integer", None))
        self.label_7.setText(QCoreApplication.translate("MainWindow", u"expect value:", None))
        self.expect_value.setPlaceholderText(QCoreApplication.translate("MainWindow", u"nonnegative number", None))
        self.label_15.setText(QCoreApplication.translate("MainWindow", u"word size:", None))
        self.word_size.setPlaceholderText(QCoreApplication.translate("MainWindow", u"positive integer", None))
        self.label_21.setText(QCoreApplication.translate("MainWindow", u"gap costs:", None))
        self.gap_costs.setPlaceholderText(QCoreApplication.translate("MainWindow", u"two positive integers separated  such as \u201c11 1\u201d", None))
        self.label_23.setText(QCoreApplication.translate("MainWindow", u"nucleotide reward:", None))
        self.nucl_reward.setPlaceholderText(QCoreApplication.translate("MainWindow", u"nonnegative number", None))
        self.label_22.setText(QCoreApplication.translate("MainWindow", u"nucleotide penalty:", None))
        self.nucl_penalty.setText("")
        self.nucl_penalty.setPlaceholderText(QCoreApplication.translate("MainWindow", u"nonpositive integer", None))
        self.submit_new_blast.setText(QCoreApplication.translate("MainWindow", u"submit new BLAST", None))
        self.load_previous_job.setText(QCoreApplication.translate("MainWindow", u"load previous job", None))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab), QCoreApplication.translate("MainWindow", u"Sequence Retrieval", None))
        self.groupBox_2.setTitle(QCoreApplication.translate("MainWindow", u"Sequences Filtering", None))
        self.control_extension.setText(QCoreApplication.translate("MainWindow", u"Extended segments refinement", None))
        self.reduce_dataset.setText(QCoreApplication.translate("MainWindow", u"Species-level sequence selection", None))
        self.label_9.setText(QCoreApplication.translate("MainWindow", u"abnormal index", None))
        self.consensus_value.setItemText(0, QCoreApplication.translate("MainWindow", u"True", None))
        self.consensus_value.setItemText(1, QCoreApplication.translate("MainWindow", u"False", None))

        self.consensus_value.setCurrentText(QCoreApplication.translate("MainWindow", u"True", None))
        self.label_5.setText(QCoreApplication.translate("MainWindow", u"length threshold:", None))
        self.len_threshold.setText(QCoreApplication.translate("MainWindow", u"100", None))
        self.label_8.setText(QCoreApplication.translate("MainWindow", u"input path:", None))
        self.in_path1.setPlaceholderText(QCoreApplication.translate("MainWindow", u"one working directory or the parent directory of multiple working directories", None))
        self.view_in_path1.setText(QCoreApplication.translate("MainWindow", u"view", None))
        self.label_14.setText(QCoreApplication.translate("MainWindow", u"output path:", None))
        self.out_path1.setPlaceholderText(QCoreApplication.translate("MainWindow", u"the same as input path by default", None))
        self.view_out_path1.setText(QCoreApplication.translate("MainWindow", u"view", None))
        self.run_filtering.setText(QCoreApplication.translate("MainWindow", u"run", None))
        self.groupBox.setTitle(QCoreApplication.translate("MainWindow", u"Sequences Alignment", None))
        self.radioButton_3.setText(QCoreApplication.translate("MainWindow", u"input one file", None))
        self.radioButton_4.setText(QCoreApplication.translate("MainWindow", u"input multiple files in one folder", None))
        self.label_13.setText(QCoreApplication.translate("MainWindow", u"input path:", None))
        self.in_path2.setPlaceholderText(QCoreApplication.translate("MainWindow", u"the path of one fasta file or the folder path that contains multiple fasta files", None))
        self.view_in_path2.setText(QCoreApplication.translate("MainWindow", u"view", None))
        self.label_16.setText(QCoreApplication.translate("MainWindow", u"output path:", None))
        self.out_path2.setPlaceholderText(QCoreApplication.translate("MainWindow", u"one folder to save the trimmeded fasta files, create a new one if does not exists", None))
        self.view_out_path2.setText(QCoreApplication.translate("MainWindow", u"view", None))
        self.label_2.setText(QCoreApplication.translate("MainWindow", u"thread:", None))
        self.ali_thr.setText(QCoreApplication.translate("MainWindow", u"-1", None))
        self.label_3.setText(QCoreApplication.translate("MainWindow", u"reorder:", None))
        self.ali_reo.setItemText(0, QCoreApplication.translate("MainWindow", u"True", None))
        self.ali_reo.setItemText(1, QCoreApplication.translate("MainWindow", u"False", None))

        self.ali_reo.setCurrentText(QCoreApplication.translate("MainWindow", u"True", None))
        self.label.setText(QCoreApplication.translate("MainWindow", u"alignment strategies:", None))
        self.ali_alg.setItemText(0, QCoreApplication.translate("MainWindow", u"auto (depends on data size)", None))
        self.ali_alg.setItemText(1, QCoreApplication.translate("MainWindow", u"add (use long sequences as backbone to align fragment sequences)", None))

        self.ali_alg.setCurrentText(QCoreApplication.translate("MainWindow", u"auto (depends on data size)", None))
        self.run_alignment.setText(QCoreApplication.translate("MainWindow", u"run", None))
        self.groupBox_4.setTitle(QCoreApplication.translate("MainWindow", u"Alignments Trimming", None))
        self.radioButton_5.setText(QCoreApplication.translate("MainWindow", u"input one file", None))
        self.radioButton_6.setText(QCoreApplication.translate("MainWindow", u"input multiple files in one folder", None))
        self.label_17.setText(QCoreApplication.translate("MainWindow", u"input path:", None))
        self.in_path3.setText("")
        self.in_path3.setPlaceholderText(QCoreApplication.translate("MainWindow", u"the path of one fasta file or the folder path that contains multiple fasta files", None))
        self.view_in_path3.setText(QCoreApplication.translate("MainWindow", u"view", None))
        self.label_26.setText(QCoreApplication.translate("MainWindow", u"output path:", None))
        self.out_path3.setPlaceholderText(QCoreApplication.translate("MainWindow", u"one folder to save the aligned fasta files, create a new one if does not exists", None))
        self.view_out_path3.setText(QCoreApplication.translate("MainWindow", u"view", None))
        self.label_6.setText(QCoreApplication.translate("MainWindow", u"automatic methods", None))
        self.tri_met.setItemText(0, QCoreApplication.translate("MainWindow", u"automated1 (heuristic selection based on similarity statistics)", None))
        self.tri_met.setItemText(1, QCoreApplication.translate("MainWindow", u"gappyout (uses information based on gaps' distribution)", None))
        self.tri_met.setItemText(2, QCoreApplication.translate("MainWindow", u"strict (automatic selection on \"strict\" mode)", None))
        self.tri_met.setItemText(3, QCoreApplication.translate("MainWindow", u"strictplus (automatic selection on \"strictplus\" mode)", None))
        self.tri_met.setItemText(4, QCoreApplication.translate("MainWindow", u"user defined method (set thresholds of non gap, similarity, consistency...)", None))

        self.tri_met.setCurrentText(QCoreApplication.translate("MainWindow", u"automated1 (heuristic selection based on similarity statistics)", None))
        self.label_48.setText(QCoreApplication.translate("MainWindow", u"non gap threshold (0-1):", None))
        self.tri_gt.setText("")
        self.tri_gt.setPlaceholderText("")
        self.label_49.setText(QCoreApplication.translate("MainWindow", u"similarity threshold (0-1):", None))
        self.tri_st.setText("")
        self.label_50.setText(QCoreApplication.translate("MainWindow", u"consistency threshold (0-1):", None))
        self.tri_ct.setText("")
        self.label_51.setText(QCoreApplication.translate("MainWindow", u"min percentage to conserve (0-100):", None))
        self.tri_con.setText("")
        self.run_trimming.setText(QCoreApplication.translate("MainWindow", u"run", None))
        self.groupBox_6.setTitle(QCoreApplication.translate("MainWindow", u"Alignments concatenation", None))
        self.label_55.setText(QCoreApplication.translate("MainWindow", u"input path:", None))
        self.in_path4.setPlaceholderText(QCoreApplication.translate("MainWindow", u"the  folder path that contains multiple fasta files", None))
        self.view_in_path4.setText(QCoreApplication.translate("MainWindow", u"view", None))
        self.label_39.setText(QCoreApplication.translate("MainWindow", u"output path:", None))
        self.out_path4.setPlaceholderText(QCoreApplication.translate("MainWindow", u"one folder to save the concatenation results, create a new one if does not exists", None))
        self.view_out_path4.setText(QCoreApplication.translate("MainWindow", u"view", None))
        self.run_concatenation.setText(QCoreApplication.translate("MainWindow", u"run", None))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab_2), QCoreApplication.translate("MainWindow", u"Supermatrix Construction", None))
        self.groupBox_3.setTitle(QCoreApplication.translate("MainWindow", u"Message Box", None))
        self.message_box.setPlaceholderText(QCoreApplication.translate("MainWindow", u"Welcome to use PyNCBIminer!", None))
        self.menuTools.setTitle(QCoreApplication.translate("MainWindow", u"Tools", None))
        self.menuInstall_dependencies.setTitle(QCoreApplication.translate("MainWindow", u"Install dependencies", None))
    # retranslateUi

