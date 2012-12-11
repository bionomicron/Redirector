# Author: Nicholas Guido
# Date: 4/4/07
# Church lab
# ReportWriter

from core.util.BasicWriter import BasicWriter
from core.util.Report import Report

class ReportWriter(BasicWriter):
    def __init__(self):
        BasicWriter.__init__(self)
        self.headerNames = ["Row Names"]
    
    def write(self, report):
        headerNames = self.headerNames
        headerNames.extend(report.returnColumnNames())
        self.writeLineArray(headerNames)
        rowNames = report.returnRowNames()
        for rowName in rowNames:
            nameArray = [rowName]
            nameArray.extend(report.returnRowArray(rowName))
            self.writeLineArray(nameArray)
            
if __name__=='__main__':
    hash = {"sad":"happy","mad":"calm"}
    stuff = Report()
    stuff.addReportColumnHash("nerd",hash)
    hash2 = {"sad":"crabby","crazy":"cool",}
    hash3 = {"bad":"foolish","great":"hysterical"}
    stuff.addReportColumnHash("loop", hash2)
    stuff.addReportColumnHash("please", hash3)
    print stuff.returnColumnNames()
    print stuff.returnRowNames()
    print stuff.returnRowArray("mad")
    print stuff.returnRowArray("sad")
    
    morestuff = ReportWriter()
    filename = 'testfile'
    morestuff.setFile(filename)
    morestuff.writer(stuff)
    morestuff.closeFile()