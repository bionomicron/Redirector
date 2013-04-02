# Author: Nicholas Guido
# Date: 3/28/07
# Church lab
# Report

from core.util.Cache import SecondOrderCache
from core.util.BasicWriter import BasicWriter

class Report:
    """
    Creates a 2nd matrix of data.
    Largely designed for reading and writing from files.
    """
    
    def __init__(self, blank='N/A'):
        self.data = SecondOrderCache()
        self.blank = blank
        
    def __setitem__(self,key,value):
        self.addColumnHash(key,value)
        
    def __getitem__(self,key):
        result = {}
        for k in self.data.keys():
            if k[1] == key:
                value = self.data[k]
                result[k[0]] = value
        return result
    
    def setEmpty(self,blank):
        """
        @param blank: the value to be used for blank fields
        @type blank: string
        self.blank = blank
        """
        self.blank = blank

    
    def add(self,rowName,columnName,value):
        """
        @param rowName: name of row
        @type rowName: string
        @param columnName: name of column
        @type columnName: string
        @param value: the value to be added
        @type value: string
        """
        self.data.addValue(rowName,columnName,str(value))
        return None
    
    def addElement(self,rowName,columnName,value):
        return self.add(rowName,columnName,value)
    
    def get(self,rowName,columnName):
        """
        @param rowName: name of row
        @type rowName: string
        @param columnName: name of column
        @type columnName: string
        @param value: the value to be added
        @type value: string
        """
        result = self.data.getValue(rowName,columnName)
        return result
    
    def getElement(self,rowName,columnName):
        return self.get(rowName,columnName)

    def addPairList(self,data):
        for ((rowName,columnName),value) in data.items():
            self.addElement(rowName,columnName,value)
        return None
    
    def addColumnHash(self,columnName,columnMap):
        """
        @param columnName: name of column to be added
        @type columnName: string
        @param columnMap: values to be added
        @type columnMap: dict
        """
        if columnMap == None:
            return None
        for rowName in columnMap.keys():
            value = columnMap[rowName]
            if type(value) != type(""):
                value = str(value)
            self.data.addValue(rowName,columnName,columnMap[rowName])
        return None
            
    def extend(self,report):
        """
        @param report: report to be appendent to this one
        @type report: Report
        """
        self.data.extend(report.data)
        
    def getRowNames(self):
        """
        @return: the names of all the rows in the report
        @rtype: string[]
        """
        return self.data.getRowKeys()
        
            
    def returnRowNames(self):
        return self.getRowNames()
        
    def getColumnNames(self):
        """
        @return: the names of all the columns in the report
        @rtype: string[]
        """
        return self.data.getColumnKeys()
        
    def returnColumnNames(self):
        return self.getColumnNames()
        
    def returnRowArray(self,rowName):
        """
        @param rowName: the name of the row of interest
        @type rowName: string
        @return: the values of the row in order from the report
        @rtype: string[]
        """
        columnNames = self.returnColumnNames()
        result = []
        for name in columnNames:
            value = self.data.getValue(rowName,name)
            if value == None:
                result.append(self.blank)
            else:
                result.append(value)
        return result
    
    def getColumn(self,name):
        '''
        @param name: name of column to return
        @return: the values of column of report
        @rtype: {string:string}
        '''
        result = {}
        if name not in self.returnColumnNames():
            return None
        for rowName in self.returnRowNames():
            value = self.data.getValue(rowName,name)
            if value == None:
                continue
            else:
                result[rowName] = value
        return result

class ReportWriter(BasicWriter):
    
    def write(self, report):
        headerNames = ['Row Names']
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
    

    
        