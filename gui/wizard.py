from PyQt5.QtWidgets import QWizard, QWizardPage
from page import Bcl2FastqPage, SnakePage

class ExomeWizard(QWizard):
	def __init__(self):
		QWizard.__init__(self)
		self.addPage(Bcl2FastqPage())
		self.addPage(SnakePage())




