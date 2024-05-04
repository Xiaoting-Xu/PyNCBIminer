# *-* coding:utf-8 *-*
# @Time:2022/5/6 15:14
# @Author: Ruijing Cheng
# @File:pyNCBIminer_00_main.py
# @Software:PyCharm


import os
import sys
from ui_main import Ui_MainWindow
from PySide2.QtWidgets import QApplication, QMessageBox, QFileDialog, QMainWindow
from PySide2.QtUiTools import QUiLoader
from PySide2.QtCore import QObject, Signal, QEventLoop, QTimer, Slot, SIGNAL
from PySide2.QtGui import QTextCursor
from pathlib import Path
import threading
from my_filter import call_miner_filter
from tools import print_line
from install_dependencies import install_mafft, install_trimal
from call_mafft2 import mafft
from call_trimal import trimal
from my_concatenation import my_concatenation
from iterated_blast import iterated_blast_main
from my_entrez import format_entrez_query, entrez_count
import warnings

root_path = os.path.dirname(os.path.realpath(sys.argv[0]))

import ssl
ssl._create_default_https_context = ssl._create_unverified_context

#from certificates import set_certifi
#set_certifi()

class EmittingStr(QObject):
    """
    define a signal of sending string
    """

    textWritten = Signal(str)

    def write(self, text):
        self.textWritten.emit(str(text))
        loop = QEventLoop()
        QTimer.singleShot(50, loop.quit)
        loop.exec_()

    def flush(self):
        pass


class MainWindow(QMainWindow):
    """
    define the main window of graphical user interface
    """
    @Slot()
    def outputWritten(self, text):
        """
        define the slot to receive the signal of sending string
        """

        cursor = self.ui.message_box.textCursor()
        cursor.movePosition(QTextCursor.End)
        cursor.insertText(text)
        self.ui.message_box.setTextCursor(cursor)
        self.ui.message_box.ensureCursorVisible()

    def __init__(self):
        # use uiloader
        # self.ui = QUiLoader().load('ui/PyNCBIminer_main9.ui')

        # use ui_main.py
        super(MainWindow, self).__init__()
        self.ui = Ui_MainWindow()
        self.ui.setupUi(self)

        # resize the window
        screen = QApplication.desktop()
        newW = int(screen.width() * 0.8)
        newH = int(newW * 0.5)
        newLeft = int(screen.width() * 0.05)
        newTop = int(screen.height() * 0.05)
        self.ui.resize(newW, newH)
        self.ui.move(newLeft, newTop)

        # connect menu with actions
        self.ui.actionInstall_MAFFT.triggered.connect(self.run_install_mafft)
        self.ui.actionInstall_trimAl.triggered.connect(self.run_install_trimal)

        # connect buttons and functions
        # buttons in Sequence Retrieving module
        self.ui.set_target_region.clicked.connect(self.set_target_region)
        self.ui.esearch.clicked.connect(self.my_esearch)
        self.ui.submit_new_blast.clicked.connect(self.submit_new_blast)
        self.ui.load_previous_job.clicked.connect(self.load_previous_job)
        self.ui.view.clicked.connect(self.view_wd)

        
        self.ui.target_region.setEditable(True)
        self.ui.set_target_region.setEnabled(False)
        target_region_list = os.listdir(Path(root_path)/Path(r"./blast_parameters"))
        target_region_list = [Path(x).stem for x in target_region_list]
        self.ui.target_region.clear()
        self.ui.target_region.addItems([""] + target_region_list)
        self.ui.target_region.currentIndexChanged.connect(self.select_target_region)

        self.ui.save_settings.clicked.connect(self.save_settings)

        self.ui.len_threshold.setEnabled(False)
        # self.ui.name_correction.setEnabled(False)
        self.ui.reduce_dataset.stateChanged.connect(self.set_reduce_threshold)
        self.ui.out_path1.setEnabled(False)
        self.ui.view_out_path1.setEnabled(False)

        self.ui.ali_alg.setEditable(False)

        # allow to do control extension, then reduce dataset
        self.ui.buttonGroup.setExclusive(False)

        self.ui.tri_met.setEditable(False)
        self.ui.tri_gt.setEnabled(False)
        self.ui.tri_st.setEnabled(False)
        self.ui.tri_ct.setEnabled(False)
        self.ui.tri_con.setEnabled(False)
        self.ui.tri_met.currentIndexChanged.connect(self.select_tri_method)

        # buttons in Supermatrix Construction module
        # filter sequences
        self.ui.view_in_path1.clicked.connect(self.view_in_path1)
        self.ui.view_out_path1.clicked.connect(self.view_out_path1)
        self.ui.run_filtering.clicked.connect(self.run_filtering2)
        # align sequences
        self.ui.view_in_path2.clicked.connect(self.view_in_path2)
        self.ui.view_out_path2.clicked.connect(self.view_out_path2)
        self.ui.run_alignment.clicked.connect(self.run_alignment)
        # trim alignments
        self.ui.view_in_path3.clicked.connect(self.view_in_path3)
        self.ui.view_out_path3.clicked.connect(self.view_out_path3)
        self.ui.run_trimming.clicked.connect(self.run_trimming)
        # concatenate alignments
        self.ui.view_in_path4.clicked.connect(self.view_in_path4)
        self.ui.view_out_path4.clicked.connect(self.view_out_path4)
        self.ui.run_concatenation.clicked.connect(self.run_concatenation)
        # # ML inference
        # self.ui.view_in_path5.clicked.connect(self.view_in_path5)
        # # self.ui.view_out_path5.clicked.connect(self.view_out_path5)
        # self.ui.view_partition_file.clicked.connect(self.view_partition_file)
        # self.ui.view_constrain_file.clicked.connect(self.view_constrain_file)
        #
        # self.ui.run_iqtree.clicked.connect(self.run_iqtree)

        # redirect the output messages to the message box
        sys.stdout = EmittingStr()
        self.ui.message_box.connect(sys.stdout, SIGNAL("textWritten(QString)"), self.outputWritten)
        sys.stderr = EmittingStr()
        self.ui.message_box.connect(sys.stderr, SIGNAL("textWritten(QString)"), self.outputWritten)

        # --------------------------------- Tools ---------------------------------
        # root_path = os.getcwd()
        # print(root_path)
        mafft_path = Path(root_path) / Path(r"./mafft/mafft-mac")
        print("Dependencies path: ")
        trimal_path = Path(root_path) / Path(r"./trimal/trimAl/bin")

        if os.path.exists(Path(root_path) / Path("./mafft/mafft-mac/mafft.bat")):
            os.environ["PATH"] = os.environ["PATH"] + ";" + str(mafft_path)
            print(mafft_path)
        elif os.system("mafft --help") == 0:
            print("mafft already in system environment.")
        else:
            print("Please check if MAFFT is installed.")

        if os.path.exists(Path(root_path) / Path("./trimal/trimAl/bin/trimal.exe")):
            os.environ["PATH"] = os.environ["PATH"] + ";" + str(trimal_path)
            print(trimal_path)
        elif os.system("trimal -h") == 0:
            print("trimal already in system environment.")
        else:
            print("Please check if trimAl is installed.")

    def run_install_mafft(self):
        self.ui.thread = threading.Thread(target=install_mafft)
        self.ui.thread.setDaemon(True)
        self.ui.thread.start()
        QMessageBox.about(self.ui, "Install MAFFT", "Installation started!")

    def run_install_trimal(self):
        self.ui.thread = threading.Thread(target=install_trimal)
        self.ui.thread.setDaemon(True)
        self.ui.thread.start()
        QMessageBox.about(self.ui, "Install trimAl", "Installation started!")

    # --------------------------------- Sequence Retrieving ---------------------------------
    def select_target_region(self):
        """

        :return:
        """
        target_region = self.ui.target_region.currentText()
        if target_region == "":
            self.ui.target_region.setEditable(True)
            self.ui.save_settings.setEnabled(True)
            self.ui.set_target_region.setEnabled(False)
        else:
            self.ui.target_region.setEditable(False)
            self.ui.save_settings.setEnabled(False)
            self.ui.set_target_region.setEnabled(True)

    def set_target_region(self):
        """
        set target region in the Sequence Retrieving module
        implemented markers includes: ITS, rbcL, matK, trnL-trnF, psbA-trnH, ndhF and rpoB
        automatically adds initial queries, BLAST parameters, key annotations and exclude sources
        :return: None
        """
        # todo: save parameters in separate files and allow users to save their settings.

        target_region = self.ui.target_region.currentText()
        if target_region == "":
            return

        parameters_dict = {}
        with open(Path(root_path)/Path("blast_parameters") / Path(target_region + ".txt"), "r") as fr:
            parameters = fr.read().splitlines()
            for parameter in parameters:
                if parameter.strip() != "":
                    parameters_dict[parameter.split("\t")[0]] = str(parameter.split("\t")[1])

        target_region = parameters_dict["target_region"]
        entrez_qualifier = parameters_dict["entrez_qualifier"]
        max_length = parameters_dict["max_length"]
        key_annotations = parameters_dict["key_annotations"]
        exclude_sources = parameters_dict["exclude_sources"]
        expect_value = parameters_dict["expect_value"]
        gap_costs = parameters_dict["gap_costs"]
        word_size = parameters_dict["word_size"]
        nucl_reward = parameters_dict["nucl_reward"]
        nucl_penalty = parameters_dict["nucl_penalty"]

        self.ui.initial_queries.clear()
        for file in os.listdir(Path(root_path)/Path("initial_queries") / Path(target_region)):
            with open(Path(root_path)/Path("initial_queries") / Path(target_region) / Path(file), "r") as fr:
                seq = fr.read()
                self.ui.initial_queries.appendPlainText(seq + "\n")

        self.ui.entrez_qualifier.setPlainText(entrez_qualifier)
        self.ui.max_length.setText(max_length)
        self.ui.expect_value.setText(expect_value)
        self.ui.gap_costs.setText(gap_costs)
        self.ui.word_size.setText(word_size)
        self.ui.nucl_reward.setText(nucl_reward)
        self.ui.nucl_penalty.setText(nucl_penalty)
        self.ui.key_annotations.setPlainText(key_annotations.replace("|", "\n"))
        self.ui.exclude_sources.setPlainText(exclude_sources.replace("|", "\n"))

        print("set target region: %s" % target_region)

    def save_settings(self):
        """

        :return:
        """
        target_region = self.ui.target_region.currentText()

        # this conflict will not happen
        # if os.path.exists(Path("parameters") / Path(target_region+".txt")):
        #     print("Setting for this target region already exists, do you want to replace it?")
        #     QMessageBox.about(self.ui, "Save settings", "Setting for this target region already exists, do you want to replace it?", QMessageBox.No)

        entrez_qualifier = self.ui.entrez_qualifier.toPlainText().strip()
        max_length = self.ui.max_length.text().strip()  # positive integer
        key_annotations = self.ui.key_annotations.toPlainText().strip()
        exclude_sources = self.ui.exclude_sources.toPlainText().strip()
        expect_value = self.ui.expect_value.text().strip()  # numeric
        gap_costs = self.ui.gap_costs.text().strip()  # two integer
        word_size = self.ui.word_size.text().strip()  # positive integer
        nucl_reward = self.ui.nucl_reward.text().strip()
        nucl_penalty = self.ui.nucl_penalty.text().strip()

        with open(Path(root_path)/Path("blast_parameters") / Path(target_region + ".txt"), "w") as fw:
            fw.write("target_region\t" + target_region + "\n")
            fw.write("entrez_qualifier\t" + entrez_qualifier + "\n")
            fw.write("max_length\t" + str(max_length) + "\n")
            fw.write("key_annotations\t" + key_annotations.replace("\n", "|") + "\n")
            fw.write("exclude_sources\t" + exclude_sources.replace("\n", "|") + "\n")
            fw.write("expect_value\t" + str(expect_value) + "\n")
            fw.write("gap_costs\t" + gap_costs + "\n")
            fw.write("word_size\t" + str(word_size) + "\n")
            fw.write("nucl_reward\t" + str(nucl_reward) + "\n")
            fw.write("nucl_penalty\t" + str(nucl_penalty) + "\n")

        initial_queries = self.ui.initial_queries.toPlainText().strip()
        os.makedirs(Path("initial_queries") / Path(target_region))
        with open(Path("initial_queries") / Path(target_region) / Path(target_region + ".fasta"), "w") as fw:
            fw.write(initial_queries)

        self.ui.target_region.addItems([target_region])
        self.ui.target_region.setEditable(False)
        self.ui.save_settings.setEnabled(False)
        self.ui.set_target_region.setEnabled(True)
        print("save settings: %s" % target_region)
        # todo: save initial queriess

    def view_wd(self):
        """
        view working directory in the Sequence Retrieving module
        :return: None
        """
        wd = QFileDialog.getExistingDirectory(self.ui, "select directory")
        self.ui.wd.setText(wd)

    def my_esearch(self):
        """
        connects the esearch button in the Sequence Retrieving module with my_esearch function
        :return: None
        """
        print_line()
        taxonomy = self.ui.taxonomy.toPlainText().strip()
        organisms = taxonomy.splitlines()
        organisms = [x for x in organisms if len(x) > 0]
        entrez_qualifier = self.ui.entrez_qualifier.toPlainText().strip()
        entrez_email = self.ui.entrez_email.text().strip()
        date_from = self.ui.date_from.text().strip()
        date_to = self.ui.date_to.text().strip()
        entrez_query = format_entrez_query(organisms=organisms, entrez_qualifier=entrez_qualifier, date_from=date_from,
                                           date_to=date_to)
        print("Entrez query: %s" % entrez_query)
        print("Entrez email: %s" % entrez_email)

        self.ui.thread = threading.Thread(target=entrez_count, args=(entrez_query, entrez_email))
        self.ui.thread.setDaemon(True)
        self.ui.thread.start()

    def submit_new_blast(self):
        """
        read the parameters in the Sequence Retrieving module
        connects the submit_new_blast button in the Sequence Retrieving module with multi_blast_main function
        :return:
        """
        print_line("*")
        print("Submitting New BlAST...")

        # read the current text of BLAST parameters in the Sequence Retrieving module
        target_region = self.ui.target_region.currentText()
        taxonomy = self.ui.taxonomy.toPlainText().strip()
        entrez_qualifier = self.ui.entrez_qualifier.toPlainText().strip()
        entrez_email = self.ui.entrez_email.text().strip()
        max_length = self.ui.max_length.text().strip()  # positive integer
        key_annotations = self.ui.key_annotations.toPlainText().strip()
        exclude_sources = self.ui.exclude_sources.toPlainText().strip()
        expect_value = self.ui.expect_value.text().strip()  # numeric
        gap_costs = self.ui.gap_costs.text().strip()  # two integer
        word_size = self.ui.word_size.text().strip()  # positive integer
        nucl_reward = self.ui.nucl_reward.text().strip()
        nucl_penalty = self.ui.nucl_penalty.text().strip()
        initial_queries = self.ui.initial_queries.toPlainText().strip()
        wd = self.ui.wd.text().strip()
        date_from = self.ui.date_from.text().strip()
        date_to = self.ui.date_to.text().strip()

        # check BLAST parameters
        # max_length, word_size: positive integer
        # expect_value: nonnegative number
        # nucl_reward: nonnegative integer
        # nucl_penalty: nonpositive integer
        # gap_costs: two positive integers separated by a space
        try:
            max_length = int(max_length)
        except ValueError:
            print("Invalid value for max_length: %s." % max_length)
            return
        else:
            if max_length < 0:
                print("max_length needs to be a positive integer")
                return

        try:
            word_size = int(word_size)
        except ValueError:
            print("Invalid value for word_size: %s." % word_size)
            print("Please visit https://ncbi.github.io/blast-cloud/dev/api.html for details about allowed values.")
            return
        else:
            if word_size < 0:
                print("word_size needs to be a positive integer")
                print("Please visit https://ncbi.github.io/blast-cloud/dev/api.html for details about allowed values.")
                return

        try:
            expect_value = float(expect_value)
        except ValueError:
            print("Invalid value for expect_value: %s." % expect_value)
            print("Please visit https://ncbi.github.io/blast-cloud/dev/api.html for details about allowed values.")
            return
        else:
            if expect_value < 0:
                print("expect_value needs to be a nonnegative number")
                print("Please visit https://ncbi.github.io/blast-cloud/dev/api.html for details about allowed values.")
                return

        try:
            nucl_reward = int(nucl_reward)
        except ValueError:
            print("Invalid value for nucl_reward: %s." % nucl_reward)
            print("Please visit https://ncbi.github.io/blast-cloud/dev/api.html for details about allowed values.")
            return
        else:
            if nucl_reward < 0:
                print("nucl_reward needs to be a nonnegative integer")
                print("Please visit https://ncbi.github.io/blast-cloud/dev/api.html for details about allowed values.")
                return

        try:
            nucl_penalty = int(nucl_penalty)
        except ValueError:
            print("Invalid value for nucl_penalty: %s." % nucl_penalty)
            print("Please visit https://ncbi.github.io/blast-cloud/dev/api.html for details about allowed values.")
            return
        else:
            if nucl_penalty > 0:
                print("nucl_penalty needs to be a nonpositive integer")
                print("Please visit https://ncbi.github.io/blast-cloud/dev/api.html for details about allowed values.")
                return

        try:
            gap_costs_list = gap_costs.split(" ")
            assert len(gap_costs_list) == 2
        except AssertionError:
            print("Invalid value for gap_costs: %s." % gap_costs)
            print("Please visit https://ncbi.github.io/blast-cloud/dev/api.html for details about allowed values.")
            return
        else:
            try:
                cost0 = int(gap_costs_list[0])
                cost1 = int(gap_costs_list[1])
            except ValueError:
                print("Invalid value for gap_costs: %s." % gap_costs)
                print("Please visit https://ncbi.github.io/blast-cloud/dev/api.html for details about allowed values.")
                return
            else:
                if cost0 < 0 or cost1 < 0:
                    print("gap_costs need to be two positive integers separated by a space")
                    print(
                        "Please visit https://ncbi.github.io/blast-cloud/dev/api.html for details about allowed values.")
                    return

        if not os.path.exists(Path(wd)):
            os.makedirs(Path(wd))
        if not os.path.exists(Path(wd) / Path("parameters")):
            os.makedirs(Path(wd) / Path("parameters"))
        if not os.path.exists(Path(wd) / Path("parameters") / Path("ref_seq")):
            os.makedirs(Path(wd) / Path("parameters") / Path("ref_seq"))
        if not os.path.exists(Path(wd) / Path("parameters") / Path("ref_msa")):
            os.makedirs(Path(wd) / Path("parameters") / Path("ref_msa"))
        if not os.path.exists(Path(wd) / Path("tmp_files")):
            os.makedirs(Path(wd) / Path("tmp_files"))
        if not os.path.exists(Path(wd) / Path("results")):
            os.makedirs(Path(wd) / Path("results"))

        organisms = taxonomy.splitlines()
        organisms = [x for x in organisms if len(x) > 0]
        entrez_query = format_entrez_query(organisms=organisms, entrez_qualifier=entrez_qualifier, date_from=date_from,
                                           date_to=date_to)
        print("Entrez query: %s" % entrez_query)
        print("Entrez email: %s" % entrez_email)
        count = entrez_count(entrez_query, entrez_email)

        # save BLAST parameters in the blast_parameters.txt file
        with open(Path(wd) / Path("parameters") / Path("blast_parameters.txt"), "w") as fw:
            fw.write("target_region\t" + target_region + "\n")
            fw.write("taxonomy\t" + taxonomy.replace("\n", "|") + "\n")
            fw.write("entrez_qualifier\t" + entrez_qualifier + "\n")
            fw.write("date_from\t" + date_from + "\n")
            fw.write("date_to\t" + date_to + "\n")
            fw.write("entrez_email\t" + entrez_email + "\n")
            fw.write("entrez_count\t" + str(count) + "\n")
            fw.write("max_length\t" + str(max_length) + "\n")
            fw.write("key_annotations\t" + key_annotations.replace("\n", "|") + "\n")
            fw.write("exclude_sources\t" + exclude_sources.replace("\n", "|") + "\n")
            fw.write("expect_value\t" + str(expect_value) + "\n")
            fw.write("gap_costs\t" + gap_costs + "\n")
            fw.write("word_size\t" + str(word_size) + "\n")
            fw.write("nucl_reward\t" + str(nucl_reward) + "\n")
            fw.write("nucl_penalty\t" + str(nucl_penalty) + "\n")

        print("BLAST parameters saved in parameters folder as blast_parameters.txt")
        if len(initial_queries) > 0:
            with open(Path(wd) / Path("parameters") / Path("initial_queries.fasta"), "w") as fw:
                fw.write(initial_queries)
        else:
            print("Please add initial queries!")
            return
        # todo: if only use one query, then no need to align

        key_annotations = key_annotations.splitlines()
        key_annotations = [x for x in key_annotations if len(x) > 0]
        exclude_sources = exclude_sources.splitlines()
        exclude_sources = [x for x in exclude_sources if len(x) > 0]

        ref_number = 5

        self.ui.thread = threading.Thread(target=iterated_blast_main, args=(wd, organisms, count,
                                                                            expect_value, gap_costs, word_size,
                                                                            nucl_reward, nucl_penalty, max_length,
                                                                            key_annotations, exclude_sources,
                                                                            ref_number,
                                                                            date_from, date_to, entrez_email))
        self.ui.thread.setDaemon(True)
        self.ui.thread.start()
        self.ui.submit_new_blast.setEnabled(False)
        self.ui.load_previous_job.setEnabled(False)
        QMessageBox.about(self.ui, "submit_new_blast", "New blast submitted!")

    def load_previous_job(self):
        """
        read BLAST parameters in the blast_parameters.txt file
        connects the load_previous_job button in the Sequence Retrieving module with multi_blast_main function
        :return:
        """
        print_line("*")
        print("Loading previous job...")

        # check if the working directory exists
        wd = self.ui.wd.text()
        if not os.path.exists(wd):
            print("Working directory does not exists, please submit new BLAST.")
            return
        if not os.path.exists(Path(wd) / Path("parameters")):
            print("Can't find parameters directory, please submit new BLAST.")
            return
        if not os.path.exists(Path(wd) / Path("tmp_files")):
            print("Can't find tmp_files directory, please submit new BLAST.")
            return
        if not os.path.exists(Path(wd) / Path("results")):
            print("Can't find results directory, please submit new BLAST.")
            return

        # read BLAST parameters in the blast_parameters.txt file
        parameters_files = os.listdir(Path(wd) / Path("parameters"))
        if "blast_parameters.txt" not in parameters_files:
            # todo: check BLAST parameters
            print("Can't find BLAST parameters, please submit new BLAST.")
            return
        if "initial_queries.fasta" not in parameters_files:
            print("Can't find initial queries, please submit new BLAST.")
            return
        else:
            self.ui.initial_queries.clear()
            with open(Path(wd) / Path("parameters") / Path("initial_queries.fasta"), "r") as fr:
                seq = fr.read()
                self.ui.initial_queries.appendPlainText(seq)

        parameters_dict = {}
        with open(Path(wd) / Path("parameters") / Path("blast_parameters.txt"), "r") as fr:
            parameters = fr.read().splitlines()
            for parameter in parameters:
                if parameter.strip() != "":
                    parameters_dict[parameter.split("\t")[0]] = str(parameter.split("\t")[1])

        target_region = parameters_dict["target_region"]
        taxonomy = parameters_dict["taxonomy"]
        entrez_qualifier = parameters_dict["entrez_qualifier"]
        date_from = parameters_dict["date_from"]
        date_to = parameters_dict["date_to"]
        entrez_email = parameters_dict["entrez_email"]
        count = parameters_dict["entrez_count"]
        max_length = int(parameters_dict["max_length"])
        key_annotations = parameters_dict["key_annotations"]
        exclude_sources = parameters_dict["exclude_sources"]
        expect_value = parameters_dict["expect_value"]
        gap_costs = parameters_dict["gap_costs"]
        word_size = parameters_dict["word_size"]
        nucl_reward = parameters_dict["nucl_reward"]
        nucl_penalty = parameters_dict["nucl_penalty"]

        # show parameters in the Sequence Retrieving panel
        # target_region_list = ["ITS", "rbcL", "matK", "trnL-trnF", "psbA-trnH", "ndhF", "rpoB"]
        target_region_list = os.listdir(Path(root_path)/Path("blast_parameters"))
        target_region_list = [Path(x).stem for x in target_region_list]
        if target_region in target_region_list:
            target_region_list.remove(target_region)
        target_region_list.insert(0, "")
        target_region_list.insert(0, target_region)
        self.ui.target_region.clear()
        self.ui.target_region.addItems(target_region_list)
        self.ui.taxonomy.setPlainText(taxonomy.replace("|", "\n"))
        self.ui.entrez_qualifier.setPlainText(entrez_qualifier)
        self.ui.date_from.setText(date_from)
        self.ui.date_to.setText(date_to)
        self.ui.entrez_email.setText(entrez_email)
        self.ui.max_length.setText(str(max_length))
        self.ui.expect_value.setText(str(expect_value))
        self.ui.gap_costs.setText(gap_costs)
        self.ui.word_size.setText(str(word_size))
        self.ui.nucl_reward.setText(str(nucl_reward))
        self.ui.nucl_penalty.setText(str(nucl_penalty))
        self.ui.key_annotations.setPlainText(key_annotations.replace("|", "\n"))
        self.ui.exclude_sources.setPlainText(exclude_sources.replace("|", "\n"))

        organisms = taxonomy.split("|")
        organisms = [x for x in organisms if len(x) > 0]
        key_annotations = key_annotations.split("|")
        key_annotations = [x for x in key_annotations if len(x) > 0]
        exclude_sources = exclude_sources.split("|")
        exclude_sources = [x for x in exclude_sources if len(x) > 0]

        # todo: show warnings when the size of queries file is zero
        queries_file_list = os.listdir(Path(wd) / Path("parameters") / Path("ref_seq"))
        if len(queries_file_list) > 0:
            round_list = []
            for queries_file in queries_file_list:
                round_list.append(int(Path(queries_file).stem.split("_")[-1]))
            round_list.sort()
            blast_round = round_list[-1]
        else:
            blast_round = 1
        ref_number = 5
        count = int(count)
        self.ui.thread = threading.Thread(target=iterated_blast_main, args=(wd, organisms, count,
                                                                            expect_value, gap_costs, word_size,
                                                                            nucl_reward, nucl_penalty, max_length,
                                                                            key_annotations, exclude_sources,
                                                                            ref_number,
                                                                            date_from, date_to, entrez_email,
                                                                            blast_round))
        self.ui.thread.setDaemon(True)
        self.ui.thread.start()
        self.ui.submit_new_blast.setEnabled(False)
        self.ui.load_previous_job.setEnabled(False)
        QMessageBox.about(self.ui, "load_previous_job", "Previous job loaded!")

    def stop_current_job(self):
        pass
        # todo: stop thread by raising exception
        # todo: release resources

    # -------------------------------- Supermatrix Construction ----------------------------------
    # view input path and output path of Sequence Filtering
    def view_in_path1(self):
        # if self.ui.buttonGroup.checkedButton().text() == "Remove exceptional records":
        #     path = QFileDialog.getOpenFileName(self.ui, "select file path", r"D:\\", "file type (*.fasta *.fas *.fa)")
        #     self.ui.in_path1.setText(path[0])
        # elif self.ui.buttonGroup.checkedButton().text() == "Combine species":
        #     path = QFileDialog.getOpenFileName(self.ui, "select file path", r"D:\\", "file type (*.fasta *.fas *.fa)")
        #     self.ui.in_path1.setText(path[0])
        # else:
        #     path = QFileDialog.getExistingDirectory(self.ui, "select file path", r"D:\\")
        #     self.ui.in_path1.setText(path)
        path = QFileDialog.getExistingDirectory(self.ui, "select file path", r"D:\\")
        self.ui.in_path1.setText(path)


    def view_out_path1(self):
        path = QFileDialog.getExistingDirectory(self.ui, "select file path", r"D:\\")
        self.ui.out_path1.setText(path)

    # view input path and output path of Sequence Alignment
    def view_in_path2(self):
        if self.ui.buttonGroup_2.checkedButton().text() == "input one file":
            path = QFileDialog.getOpenFileName(self.ui, "select file path", r"D:\\", "file type (*.fasta *.fas *.fa)")
            self.ui.in_path2.setText(path[0])
        else:
            path = QFileDialog.getExistingDirectory(self.ui, "select file path", r"D:\\")
            self.ui.in_path2.setText(path)

    def view_out_path2(self):
        path = QFileDialog.getExistingDirectory(self.ui, "select file path", r"D:\\")
        self.ui.out_path2.setText(path)

    # view input path and output path of Alignments Trimming
    def view_in_path3(self):
        if self.ui.buttonGroup_3.checkedButton().text() == "input one file":
            path = QFileDialog.getOpenFileName(self.ui, "select file path", r"D:\\", "file type (*.fasta *.fas *.fa)")
            self.ui.in_path3.setText(path[0])
        else:
            path = QFileDialog.getExistingDirectory(self.ui, "select file path", r"D:\\")
            self.ui.in_path3.setText(path)

    def view_out_path3(self):
        path = QFileDialog.getExistingDirectory(self.ui, "select file path", r"D:\\")
        self.ui.out_path3.setText(path)

    # view input path and output path of Alignments Concatenation
    def view_in_path4(self):
        path = QFileDialog.getExistingDirectory(self.ui, "select file path", r"D:\\")
        self.ui.in_path4.setText(path)

    def view_out_path4(self):
        path = QFileDialog.getExistingDirectory(self.ui, "select file path", r"D:\\")
        self.ui.out_path4.setText(path)

    def view_in_path5(self):
        path = QFileDialog.getOpenFileName(self.ui, "select file path", r"D:\\",
                                           "file type (*.fasta *.fas *.fa *.phylip *.phy)")
        self.ui.in_path5.setText(path[0])

    def set_reduce_threshold(self):
        """

        :return:
        """
        if self.ui.reduce_dataset.isChecked():
            self.ui.len_threshold.setEnabled(True)
            # self.ui.name_correction.setEnabled(True)
        else:
            self.ui.len_threshold.setEnabled(False)
            # self.ui.name_correction.setEnabled(False)

    def run_filtering2(self):
        in_path1 = self.ui.in_path1.text().strip()
        # out_path1 = self.ui.out_path1.text().strip()
        out_path1 = in_path1
        len_threshold = int(self.ui.len_threshold.text().strip())
        # name_correction = eval(self.ui.name_correction.currentText().strip())
        name_correction = False

        if self.ui.control_extension.isChecked() and self.ui.reduce_dataset.isChecked():
            print("Control extension and reduce dataset...")
            action = 3
            self.ui.thread = threading.Thread(target=call_miner_filter,
                                              args=(in_path1, out_path1, action, len_threshold, name_correction))
            self.ui.thread.setDaemon(True)
            self.ui.thread.start()
        elif self.ui.control_extension.isChecked():
            print("Control extension...")
            action = 1  # call_miner_filter(in_path, out_path, action, len_shresh, max_num)
            self.ui.thread = threading.Thread(target=call_miner_filter,
                                              args=(in_path1, out_path1, action, len_threshold, name_correction))
            self.ui.thread.setDaemon(True)
            self.ui.thread.start()
        elif self.ui.reduce_dataset.isChecked():
            print("Reduce dataset...")
            action = 2
            self.ui.thread = threading.Thread(target=call_miner_filter,
                                              args=(in_path1, out_path1, action, len_threshold, name_correction))
            self.ui.thread.setDaemon(True)
            self.ui.thread.start()

        else:
            print("Please select one option.")                                                

    def run_alignment(self):
        """
        connects the running button of Sequence Alignment with the call_mafft function
        :return: None
        """
        in_path2 = self.ui.in_path2.text().strip()
        out_path2 = self.ui.out_path2.text().strip()
        if not os.path.exists(out_path2):
            os.makedirs(out_path2)
        # ali_mod = eval(self.ui.ali_mod.currentText().strip())
        # ali_cmd = self.ui.ali_cmd.toPlainText().strip()
        # ali_add_cho = self.ui.ali_add_cho.currentText().strip()
        # ali_add_path = self.ui.ali_add_path.text().strip()
        ali_alg = self.ui.ali_alg.currentText().strip().split(" ")[0]  # algorithm, auto, add
        ali_thr = self.ui.ali_thr.text().strip()  # thread
        ali_reo = eval(self.ui.ali_reo.currentText().strip())  # reorder
        # ali_add_par = self.ui.ali_add_par.text().strip()
        ali_mod = False
        ali_cmd = ""
        ali_add_cho = ""
        ali_add_path = ""
        ali_add_par = ""

        """in_path, out_path = '', add_choice = '', add_path = '', algorithm = 'auto',
        thread = -1, reorder = True, additional_params = '',
        pure_command_mode = False, pure_command = ''"""
        self.ui.thread = threading.Thread(target=mafft, args=(in_path2, out_path2, ali_add_cho, ali_add_path,
                                                              ali_alg, ali_thr, ali_reo, ali_add_par, ali_mod, ali_cmd))
        print("Running alignment...")
        self.ui.thread.setDaemon(True)
        self.ui.thread.start()

        # print(yiyang.call_mafft.mafft(in_path=in_path2, out_path=out_path2, add_choice=ali_add_cho, add_path=ali_add_path,
        #                               algorithm=ali_alg, thread=ali_thr, reorder=ali_reo, additional_params=ali_add_par,
        #                               pure_command_mode=ali_mod, pure_command=ali_cmd))

    def select_tri_method(self):
        """

        :return:
        """
        tri_method = self.ui.tri_met.currentText()
        if tri_method == "user defined method (set thresholds of non gap, similarity, consistency...)":
            self.ui.tri_gt.setEnabled(True)
            self.ui.tri_st.setEnabled(True)
            self.ui.tri_ct.setEnabled(True)
            self.ui.tri_con.setEnabled(True)
        else:
            self.ui.tri_gt.clear()
            self.ui.tri_st.clear()
            self.ui.tri_ct.clear()
            self.ui.tri_con.clear()
            self.ui.tri_gt.setEnabled(False)
            self.ui.tri_st.setEnabled(False)
            self.ui.tri_ct.setEnabled(False)
            self.ui.tri_con.setEnabled(False)

    def run_trimming(self):
        """
        connects the running button of Alignments Trimming with the call_trimal function
        :return: None
        """
        in_path3 = self.ui.in_path3.text().strip()
        out_path3 = self.ui.out_path3.text().strip()
        if not os.path.exists(out_path3):
            os.makedirs(out_path3)
        tri_met = self.ui.tri_met.currentText().strip().split(" ")[0]
        if tri_met == "user":
            tri_met = ""
        # tri_htm = eval(self.ui.tri_htm.currentText().strip())
        # tri_bpl = eval(self.ui.tri_bpl.currentText().strip())
        tri_htm = True
        tri_bpl = False
        tri_gt = self.ui.tri_gt.text().strip()
        tri_st = self.ui.tri_st.text().strip()
        tri_ct = self.ui.tri_ct.text().strip()
        tri_con = self.ui.tri_con.text().strip()
        # tri_add_par = self.ui.ali_add_par.text().strip()
        # tri_mod = eval(self.ui.tri_mod.currentText().strip())
        # tri_cmd = self.ui.tri_cmd.toPlainText().strip()
        tri_add_par = ""
        tri_mod = False
        tri_cmd = ""
        """
            in_path, out_path='',
           htmlout=True, bp_length=False,
           implement_methods='automated1', gt='', st='', ct='', cons='',
           additional_params='', pure_command_mode=False, pure_command=''
        """
        self.ui.thread = threading.Thread(target=trimal, args=(in_path3, out_path3,
                                                               tri_htm, tri_bpl, tri_met,
                                                               tri_gt, tri_st, tri_ct, tri_con,
                                                               tri_add_par, tri_mod, tri_cmd))
        print("Running trimming...")
        self.ui.thread.setDaemon(True)
        self.ui.thread.start()

        # print(yiyang.call_trimal.trimal(in_path=in_path3, out_path=out_path3,
        #                          htmlout=tri_htm, bp_length=tri_bpl, implement_methods=tri_met,
        #                          gt=tri_gt, st=tri_st, ct=tri_ct, cons=tri_con,
        #                          additional_params=tri_add_par, pure_command_mode=tri_mod, pure_command=tri_cmd))

    def run_concatenation(self):
        """
        connects the running button of Alignments Concatenation with the my_concatenation function
        :return: None
        """
        in_path4 = self.ui.in_path4.text().strip()
        out_path4 = self.ui.out_path4.text().strip()
        if not os.path.exists(out_path4):
            os.makedirs(out_path4)
        self.ui.thread = threading.Thread(target=my_concatenation, args=(in_path4, out_path4))
        print("Running concatenation...")
        self.ui.thread.setDaemon(True)
        self.ui.thread.start()


if __name__ == "__main__":
   # todo: print cannot show?
    # (_OLD_VIRTUAL_PATH)
    # root_path = os.path.abspath(os.path.dirname(__file__))  # running dir

    app = QApplication([])
    main_window = MainWindow()
    main_window.show()
    sys.exit(app.exec_())

    # use uiload
    # app = QApplication([])
    # main_window = MainWindow()
    # main_window.ui.show()
    # sys.exit(app.exec_())




