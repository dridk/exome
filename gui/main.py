import sys 
from PyQt5 import QtWidgets
from wizard import ExomeWizard

if __name__ == "__main__":
	app = QtWidgets.QApplication(sys.argv)

	w = ExomeWizard()

	w.show()

	sys.exit(app.exec_())