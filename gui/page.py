from PyQt5.QtWidgets import * 
from PyQt5.QtCore import *
from enum import Enum
import webbrowser
import csv 
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
        self.registerField("bclpath", self.bclPathWidget, "path")
        self.registerField("sheetpath", self.bclSheetWidget, "path")


#------------------------------------------
class RunPage(QWizardPage):
    def __init__(self, cmd = str() , args = [] ):

        self.cmd  = cmd
        self.args = args 
        self.runComplete = True

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
        self.openFolder.setDisabled(False)

    def __run(self):
        self.proc = QProcess(self)
        self.proc.setProcessChannelMode(QProcess.MergedChannels)
        self.proc.readyReadStandardOutput.connect(self.__addoutput)
        self.proc.finished.connect(self.__runComplete)
        self.proc.start(self.cmd, self.args)
        self.console.appendPlainText(self.cmd + " " + " ".join(self.args))

    def __stop(self):
        print("stop process")
        self.proc.kill()
        QProcess.execute("killall bcl2fastq")
        QProcess.execute("killall gzip")
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
        self.cmd  = "../scripts/createfastq.sh"
        self.args = [self.field("bclpath"),self.field("sheetpath")]




#------------------------------------------
class SnakePage(QWizardPage):
    def __init__(self):
        QWizardPage.__init__(self)
        self.fastqpath     = PathLineEdit(PathLineEdit.PathType.FOLDER)
        self.bedpath       = PathLineEdit(PathLineEdit.PathType.FILE)
        self.sampleModel   = QStringListModel()
        self.sampleView    = QListView()
        self.sampleView.setModel(self.sampleModel)


        formLayout = QFormLayout()
        formLayout.addRow("fastq folder", self.fastqpath)
        formLayout.addRow("target bed file", self.bedpath)

        mainLayout = QVBoxLayout()
        mainLayout.addLayout(formLayout)
        mainLayout.addWidget(self.sampleView)

        self.setLayout(mainLayout)

    def __loadSample(self):
        pass

    ''' override '''
    def initializePage(self):
        # Open and read sheet file 
        samples = []
        sheetpath = self.field("sheetpath")
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
                    samples.append(str(row[0]))

        self.sampleModel.setStringList(samples)


                    
                        








#------------------------------------------
class SnakeRunPage(RunPage):
    def __init__(self):
        RunPage.__init__(self)

