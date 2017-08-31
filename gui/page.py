from PyQt5.QtWidgets import * 
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

		if self.type == PathLineEdit.PathType.FILE:
			filename = QFileDialog.getOpenFileName(self, "file path", "path")
				
#------------------------------------------

class Bcl2FastqPage(QWizardPage):
	def __init__(self):
		QWizardPage.__init__(self)

		# Construct Form 
		self.bclPathWidget  = PathLineEdit(PathLineEdit.PathType.FOLDER)
		self.bclSheetWidget = PathLineEdit(PathLineEdit.PathType.FILE)
		formLayout = QFormLayout()
		formLayout.addRow("bcl path", self.bclPathWidget)
		formLayout.addRow("datasheet path", self.bclSheetWidget)

		# Construct console output 
		self.console = QPlainTextEdit()
		groupBox  = QGroupBox()
		groupBox.setTitle("output")
		boxLayout = QVBoxLayout()
		boxLayout.setContentsMargins(0, 0, 0, 0)
		boxLayout.addWidget(self.console) 
		groupBox.setLayout(boxLayout)

		# Run button 
		self.setButtonText(QWizard.CustomButton1, "Run")
		self.setCommitPage(True)

		mainLayout = QVBoxLayout()
		mainLayout.addLayout(formLayout)
		mainLayout.addWidget(groupBox)


		self.setLayout(mainLayout)





# ==========================================
class SnakePage(QWizardPage):
	def __init__(self):
		QWizardPage.__init__(self)