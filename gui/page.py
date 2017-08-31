from PyQt5.QtWidgets import * 
from PyQt5.QtCore import QObject, pyqtProperty, QProcess
from enum import Enum
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
            self.edit.setText(filename)


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
        self.bclPathWidget  = PathLineEdit(PathLineEdit.PathType.FOLDER)
        self.bclSheetWidget = PathLineEdit(PathLineEdit.PathType.FILE)
        formLayout = QFormLayout()
        formLayout.addRow("bcl path", self.bclPathWidget)
        formLayout.addRow("datasheet path", self.bclSheetWidget)
        mainLayout = QVBoxLayout()
        mainLayout.addLayout(formLayout)
        self.setLayout(mainLayout)

        # register field 
        self.registerField("bclpath", self.bclPathWidget, "path")
        self.registerField("sheetpath", self.bclSheetWidget, "path")


#------------------------------------------
class RunPage(QWizardPage):
    def __init__(self, cmd : str, args:list ):

        self.cmd  = cmd
        self.args = args 

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
        spacer = QWidget()
        spacer.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Preferred)
        runLayout = QHBoxLayout()
        runLayout.addWidget(spacer)
        runLayout.addWidget(self.runButton)

        layout = QVBoxLayout()
        layout.addWidget(groupBox)
        layout.addLayout(runLayout)
        self.setLayout(layout)


    ''' override '''
    def initializePage(self):
        print("Init page")
        self.console.clear()
        self.console.appendPlainText("bcl folder : "+self.field("bclpath"))
        self.console.appendPlainText("data sheet : "+self.field("sheetpath"))
        self.runComplete = False

    ''' override ''' 
    def isComplete(self):
         return self.runComplete

       
    def __addoutput(self):
        self.console.appendPlainText(bytes(self.proc.readAllStandardOutput()).decode("utf-8"))

    def __runComplete(self):
        self.runComplete = True 
        self.completeChanged.emit()

    def __run(self):
        self.proc = QProcess(self)
        self.proc.setProcessChannelMode(QProcess.MergedChannels)
        self.proc.readyReadStandardOutput.connect(self.__addoutput)
        self.proc.finished.connect(self.__runComplete)
        self.proc.start(self.cmd, self.args)

#------------------------------------------
class SnakePage(QWizardPage):
    def __init__(self):
        QWizardPage.__init__(self)

