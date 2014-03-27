'''
Author: Graham Rockwell
Date: 9/28/13
Church lab
Report
'''

from util.Cache import SecondOrderCache
from util.BasicWriter import BasicWriter

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
    
    def getRow(self,name):
        '''
        @param name: name of column to return
        @return: the values of column of report
        @rtype: {string:string}
        '''
        result = {}
        if name not in self.getRowNames():
            return None
        for columnName in self.getColumnNames():
            value = self.data.getValue(name,columnName)
            if value == None:
                continue
            else:
                result[columnName] = value
        return result
    
    def getColumn(self,name):
        '''
        @param name: name of column to return
        @return: the values of column of report
        @rtype: {string:string}
        '''
        result = {}
        if name not in self.getColumnNames():
            return None
        for rowName in self.getRowNames():
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
        

from util.FlatFileParser import FlatFileParser
from optparse import OptionParser
from time import time,strftime
    
if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option("-v", "--verbose", action="store_true", dest="verbose", default=False, help="set verbose mode")
    parser.add_option("-f", "--file", dest="input_file", default=None, help="test input file")
    parser.add_option("-o", "--output", dest="output_file", default=None, help="test output file")
    (options,args) = parser.parse_args()  
    
    input_file = options.input_file
    output_file = options.output_file
    
    input_file = 'Variance__Report.txt'
    input_file = 'Variance_FC_00707_Report.txt'
    output_file = 'Original_Report_out_test.txt'
    
    print "Reading report file [%s]" %(input_file)
    s_read_time = time()    
    reader = FlatFileParser()
    report = reader.parseToReport(input_file, keyTag='', header=None, unique=True)
    e_read_time = time() - s_read_time
    print "Reading time [%s]" %(e_read_time)
    
    new_report = Report()
    
    c_names = report.getColumnNames()
    r_names = report.getRowNames()
    s_reload_time = time()
    for r_name in r_names:
        print "reloading [%s]" % (r_name)
        for c_name in c_names:
            value = report.getElement(r_name, c_name)
            if value != None:
                new_report.addElement(r_name, c_name, value)
                
    e_reload_time = time() - s_reload_time
    print "reload time [%s]" % (e_reload_time)
    
    
    s_write_time = time()
    writer = ReportWriter()
    print "Writing report to [%s]" % (output_file)
    writer = ReportWriter()
    writer.setFile(output_file)
    writer.write(new_report) 
    writer.closeFile()
    e_write_time = time() - s_write_time
    print "Writing time [%s]" % (e_write_time)
            
    
    

    
        