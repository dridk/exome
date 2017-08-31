from PyQt5.QtWidgets import QWizard, QWizardPage
from page import BclPage, RunPage, SnakePage

class ExomeWizard(QWizard):
    def __init__(self):
        QWizard.__init__(self)



        self.addPage(BclPage())
        self.addPage(RunPage("ls",["-l","/home/sacha"]))

        self.addPage(SnakePage())
        self.addPage(RunPage("cal",["1984"]))





