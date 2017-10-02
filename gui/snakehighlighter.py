from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *

class SnakeHighlighter(QSyntaxHighlighter):
    def __init__(self, doc):
        QSyntaxHighlighter.__init__(self, doc)

        self.rules = []  # (regexp : textformat) tuple

        ruleRegx = QRegularExpression("^[^rule|\\s].+")
        ruleForm = QTextCharFormat()
        ruleForm.setForeground(Qt.darkGreen)

        cmdRegx = QRegularExpression("^[rule|\\s].+")
        cmdForm = QTextCharFormat()
        cmdForm.setForeground(Qt.darkGray)


        self.rules.append((ruleRegx, ruleForm))
        self.rules.append((cmdRegx, cmdForm))

    ''' override '''
    def highlightBlock(self, text):


        for rule in self.rules:
            iterator = rule[0].globalMatch(text)
            while iterator.hasNext():
                match = iterator.next()
                self.setFormat(match.capturedStart(), match.capturedLength(), rule[1])













        