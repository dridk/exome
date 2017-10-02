from PyQt5.QtWidgets import * 
from PyQt5.QtCore import *
from enum import Enum
import webbrowser
import csv 
import os
import glob
import yaml
from snakehighlighter import *

#------------------------------------------
class PathLineEdit(QWidget):

    class PathType(Enum):
        FILE   = 1
        FOLDER = 2 


    def __init__(self, type = PathType.FILE):
        QWidget.__init__(self)
        self.type = type
        layout = QHBoxLayout()
        self.edit   = QLineEdit()
        self.button = QPushButton("Browse")
        layout.addWidget(self.edit)
        layout.addWidget(self.button) 
        layout.setContentsMargins(0, 0, 0, 0)
        self.setLayout(layout)
        self.button.clicked.connect(self.__browse)

    def __browse(self):
        if self.type == PathLineEdit.PathType.FOLDER:
            filename = QFileDialog.getExistingDirectory(self, "folder path", "path")
            self.edit.setText(filename)


        if self.type == PathLineEdit.PathType.FILE:
            filename = QFileDialog.getOpenFileName(self, "file path", "path")
            self.edit.setText(str(filename[0]))

    @pyqtProperty(str)
    def toolTip(self):
        return self.edit.toolTip()

    @toolTip.setter
    def toolTip(self, value):
        self.edit.setToolTip(value)

    @pyqtProperty(str)
    def placeholder(self):
        return self.edit.placeholderText()

    @placeholder.setter
    def placeholder(self, value):
        self.edit.setPlaceholderText(value)


    def lineEdit(self):
        return self.edit


    @pyqtProperty(str)
    def path(self):
        return self.edit.text()

    @path.setter
    def path(self, value):
        self.edit.setText(value)

                
#------------------------------------------

class BclPage(QWizardPage):
    def __init__(self):
        QWizardPage.__init__(self)

        # Construct Form 
        self.description    = QLabel()
        self.bclPathWidget  = PathLineEdit(PathLineEdit.PathType.FOLDER)
        self.bclSheetWidget = PathLineEdit(PathLineEdit.PathType.FILE)

        self.bclPathWidget.toolTip  = "Folder path to the illumina output"
        self.bclSheetWidget.toolTip = "Data Sheet used for illumina project (*.csv)"

        self.bclPathWidget.placeholder  = "illumina folder ..."
        self.bclSheetWidget.placeholder = "data Sheet ... (*.csv)"
        
        
        formLayout = QFormLayout()
        formLayout.addRow("bcl path", self.bclPathWidget)
        formLayout.addRow("datasheet path", self.bclSheetWidget)
        mainLayout = QVBoxLayout()
        mainLayout.addWidget(self.description)
        mainLayout.addLayout(formLayout)
        self.setLayout(mainLayout)

        self.description.setText("<b>Convert BCL illumina to Fastq </b>")


        # register field 
        self.registerField("bclpath", self.bclPathWidget.lineEdit())
        self.registerField("sheetpath", self.bclSheetWidget.lineEdit())

        # Temp 
        self.bclPathWidget.path = "/DATAS/illumina/170802_NB501647_0006_AHHT5TAFXX"
        self.bclSheetWidget.path = "/DATAS/illumina/170802_NB501647_0006_AHHT5TAFXX/sheet.csv"



#------------------------------------------
class RunPage(QWizardPage):
    def __init__(self, cmd = str() , args = [] ):

        self.cmd  = cmd
        self.args = args 
        self.runComplete = True
        self.toKill = []

        # box layout 
        QWizardPage.__init__(self)
        groupBox = QGroupBox()
        groupBox.setTitle("Output")
        self.console = QPlainTextEdit()
        self.console.setReadOnly(True)
        groupLayout = QVBoxLayout()
        groupLayout.addWidget(self.console)
        groupLayout.setContentsMargins(0, 0, 0, 0)
        groupBox.setLayout(groupLayout)

        self.runButton = QPushButton("Run")
        self.runButton.clicked.connect(self.__run)
        self.cancelButton = QPushButton("Stop")
        self.cancelButton.clicked.connect(self.__stop)
        self.openFolder = QPushButton("Open")
        self.openFolder.setDisabled(True)
        self.openFolder.clicked.connect(self.__open)

        spacer = QWidget()
        spacer.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Preferred)
        runLayout = QHBoxLayout()
        runLayout.addWidget(self.cancelButton)
        runLayout.addWidget(spacer)
        runLayout.addWidget(self.runButton)
        runLayout.addWidget(self.openFolder)


        layout = QVBoxLayout()
        layout.addWidget(groupBox)
        layout.addLayout(runLayout)
        self.setLayout(layout)

    ''' override '''
    def initializePage(self):
        print("Init page")
        self.console.clear()
        self.runComplete = False

    ''' override ''' 
    def isComplete(self):
         return self.runComplete

       
    def __addoutput(self):
        self.console.appendPlainText(bytes(self.proc.readAllStandardOutput()).decode("utf-8"))

    def __runComplete(self):
        self.runComplete = True 
        self.completeChanged.emit()
        self.openFolder.setEnabled(True)
        self.runButton.setEnabled(True)
        #self.cancelButton.setEnabled(False)



    def __run(self):
        self.proc = QProcess(self)
        self.proc.setProcessChannelMode(QProcess.MergedChannels)
        self.proc.readyReadStandardOutput.connect(self.__addoutput)
        self.proc.finished.connect(self.__runComplete)
        self.proc.start(self.cmd, self.args)
        self.console.appendPlainText(self.cmd + " " + " ".join(self.args))
        self.openFolder.setEnabled(False)
        self.runButton.setEnabled(False)
        #self.cancelButton.setEnabled(True)


    def __stop(self):
        self.openFolder.setEnabled(True)
        self.runButton.setEnabled(True)
        print("stop process")
        self.proc.kill()
        self.proc.terminate()

        for p in self.toKill:
            os.system("killall -9 "+p);

        self.console.clear()


    def __open(self):
        webbrowser.open(self.args[0]+"/output")



#------------------------------------------
class BclRunPage(RunPage):
    def __init__(self):
        RunPage.__init__(self)

    ''' override '''
    def initializePage(self):
        self.console.clear()
     
        self.cmd  = os.path.dirname(os.path.realpath(__file__))+"/createfastq.sh"
        self.args = [self.field("bclpath"),self.field("sheetpath")]
        self.toKill = ["bcl2fastq", "gzip"]


#------------------------------------------
class SnakePage(QWizardPage):
    def __init__(self):
        QWizardPage.__init__(self)
        self.fastqpath     = PathLineEdit(PathLineEdit.PathType.FOLDER)
        self.bedpath       = PathLineEdit(PathLineEdit.PathType.FILE)
        self.sampleView    = QTreeWidget()
        self.sampleView.setHeaderLabels(["sample"])
        self.sampleView.setSelectionMode(QAbstractItemView.MultiSelection)
        self.samplesList = []
        self.sampleView.setContextMenuPolicy(Qt.ActionsContextMenu)

        ## add whole check action 
        checkAction   = QAction("check all", self)
        uncheckAction = QAction("un check all", self)

        checkAction.triggered.connect(self.checkAll)
        uncheckAction.triggered.connect(self.uncheckAll)

        self.sampleView.addAction(checkAction)
        self.sampleView.addAction(uncheckAction)

        


        formLayout = QFormLayout()
        formLayout.addRow("fastq folder", self.fastqpath)
        formLayout.addRow("target bed file", self.bedpath)

        mainLayout = QVBoxLayout()
        mainLayout.addLayout(formLayout)
        mainLayout.addWidget(self.sampleView)

        self.setLayout(mainLayout)

        self.registerField("fastqpath*", self.fastqpath.lineEdit())
        self.registerField("bedpath*", self.bedpath.lineEdit())
        self.registerField("samples", self, "samples")

    @pyqtProperty(list)
    def samples(self):
        return self.samplesList 

    @samples.setter
    def samples(self, values):
        self.samplesList = values

    ''' override '''
    def initializePage(self):
        self.samplesList.clear()
        self.samples = []
        self.fastqpath.path = self.field("bclpath")+"/output"
        # Open and read sheet file 
        self.__loadSample()

    ''' override '''
    def validatePage(self):
        self.samples = []
        for i in range(0,self.sampleView.topLevelItemCount()):
            item = self.sampleView.topLevelItem(i)
            if item.checkState(0) == Qt.Checked:
                self.samples.append(str(item.data(0, Qt.UserRole)))

        return len(self.samples) != 0

    def checkAll(self):
        for i in range(0,self.sampleView.topLevelItemCount()): 
            self.sampleView.topLevelItem(i).setCheckState(0,Qt.Checked)

    def uncheckAll(self):
        for i in range(0,self.sampleView.topLevelItemCount()): 
            self.sampleView.topLevelItem(i).setCheckState(0,Qt.Unchecked)

    def __findFastq(self,name):
        return glob.glob(os.path.join(self.fastqpath.path , name+"*.fastq.gz"))


    def __loadSample(self):
        self.samples = []
        self.sampleView.clear()
        sheetpath = self.field("sheetpath")
        if os.path.isfile(sheetpath) is False:
            return

        with open(sheetpath, 'r') as csvfile:
            reader = csv.reader(csvfile)
            save = False
            for row in reader:
                if len(row) > 0:
                    if row[0] == "[Data]":
                        save = True 
                        row = next(reader)
                        row = next(reader)
                if save == True and len(row) > 2:
                    name = row[0]
                    item = QTreeWidgetItem()
                    item.setCheckState(0,Qt.Checked)
                    fastqCount = len(self.__findFastq(name))
                    item.setText(0, "{} ({})".format(name,fastqCount))
                    item.setData(0, Qt.UserRole, name)
                    
                    # if fastqCount >= 9:
                    #     item.setForeground(0,QColor("#00cc00"))
                    # else:
                    #     item.setForeground(0,QColor("#ff6666"))

                    for fq in self.__findFastq(name):
                        sub = QTreeWidgetItem()
                        sub.setText(0,fq)
                        item.addChild(sub)

                    self.sampleView.addTopLevelItem(item)


#------------------------------------------
class SnakeRunPage(RunPage):
    def __init__(self):
        RunPage.__init__(self)
        #self.snakeSyntax = SnakeHighlighter(self.console.document())
        self.configPath = os.path.dirname(os.path.realpath(__file__))+"/../config.yml"
        self.toKill = ["snakemake"]


        ''' override '''
    def initializePage(self):
        self.console.clear()
        self.cmd  = os.path.dirname(os.path.realpath(__file__))+"/../run.sh"
        self.args = ["/tmp/myconfig.yml"]

        self.__createConfig()

    def cleanupPage(self):
        self.console.clear()


    def __createConfig(self): 
        docs = yaml.load_all(open(self.configPath, "r"))
        newCfg = {}
        for doc in docs:
            for k,v in doc.items():
                newCfg[k] = v

        newCfg["SAMPLES"] = self.field("samples")
        newCfg["TARGET_BED_FILE"] = self.field("bedpath")
        newCfg["RAW_FOLDER"] = self.field("fastqpath")
        with open(self.args[0],"w") as newFile:
            yaml.dump(newCfg, newFile, default_flow_style=False)
      
        with open(self.args[0],"r") as newFile:
            for line in newFile:
                self.console.appendPlainText(line)


 