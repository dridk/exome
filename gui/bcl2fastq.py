import sys 
from PyQt5.QtWidgets import * 
from page import * 

class Wizard(QWizard):
    def __init__(self):
        QWizard.__init__(self)
        self.addPage(BclPage())
        self.addPage(BclRunPage())


if __name__ == "__main__":
    app = QApplication(sys.argv)
    w = Wizard()
    w.show()
    sys.exit(app.exec_())