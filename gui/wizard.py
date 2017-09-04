from PyQt5.QtWidgets import QWizard, QWizardPage
from page import *

class ExomeWizard(QWizard):
    def __init__(self):
        QWizard.__init__(self)



        # self.addPage(BclPage())
        # self.addPage(BclRunPage())

        self.addPage(SnakePage())
        self.addPage(SnakeRunPage())





