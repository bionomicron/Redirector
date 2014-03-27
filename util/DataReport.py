'''
Created on Mar 21, 2014

@author: Graham Rockwell
'''

import collections

class FastReport(object):
    '''
    2 Dimensional flat file data holding object.
    +Trying to improve access and size, arrays for each row indexed to a column map.
    -This may make adding data by columns slow.
    '''

    def __init__(self,column_names=[]):
        '''
        Constructor
        '''
        self._column_map = collections.OrderedDict()
        self._data = {} #map from row to array of values
        self._extendColumns(column_names)
            
    def _extendColumns(self,names=[]):
        for name in names:
            if name not in self._column_map.keys():
                self._column_map[name] = len(self._column_map)
                   
    def getRowNames(self):
        return self._data.keys()
    
    def getColumnNames(self):
        return self._column_map.keys()
        
    def get(self, row_name='', column_name=''):
        return self._data[row_name][self._column_map[column_name]]

    def add(self, row_name, column_name, value):
        if not row_name in self._data:
            self._data[row_name] = [None] * len(self._column_map)
        self._data[row_name][self._column_map[column_name]] = value
        
class SafeReport(FastReport):
    '''
    Tries to have extension and checks for missing row / column keys
    '''
    
    def _updateRow(self,row_name):
        length_diff = len(self._column_map) - len(self._data[row_name])
        if length_diff > 0:
            self._data[row_name].append([None]*length_diff)    
    
    def _updateReport(self):
        for row_name in self._data.keys():
            self.updateRow(row_name)    
    
    def add(self, row_name, column_name, value):
        self._extendColumns([column_name])
        
        if not row_name in self._data:
            self._data[row_name] = [None] * len(self._column_map)
        
        self._updateRow(row_name)

        self._data[row_name][self._column_map[column_name]] = value
 
    def get(self, row_name, column_name):
        if row_name not in self._data.keys() or column_name not in self._column_map.keys():
            return None
        
        self._updateRow(row_name)
        
        return self._data[row_name][self._column_map[column_name]]

    def getColumnNames(self):
        return self._column_map.keys()
    
    def getRowNames(self):
        return self._data.keys()
    
    def addRow(self,row_name,values={}):
        '''
        Safe but slow addition of row to report
        '''
        for (column_name,value) in values:
            self.add(row_name,column_name,value)
    
    def addColumn(self,column_name,values={}):
        '''
        Safe but slow addition of column to report
        '''
        for (row_name,value) in values:
            self.add(row_name,column_name,value)
    
    def extend(self,report):
        '''
        Safe but slow extension of report
        '''
        self._extendColumns(report._column_map.keys())
        for (row_name,values) in report._data.keys():
            for (column_name,col_index) in report._column_map.items():
                if values[col_index] != None:
                    self.add(row_name,column_name,values[col_index])
        return None 
    
import csv
    
class DataReportIO:
    
    def __init__(self,header=None):
        self.header=header
    
    def _check_header(self,header):
        pass
    
    def readReport(self,file_name,key_column,delimiter="\t",quotechar="'"):
        result = SafeReport()
        with open(file_name, 'rb') as csvfile:
            reader = csv.reader(csvfile, delimiter=delimiter, quotechar=quotechar)
            row = reader.next()
            result._extendColumns(row)
            
            if key_column in result._column_map.keys():
                row_key_index = result._column_map.key(key_column)
            else:
                row_key_index = 0
                
            for row in reader:
                row_name = row[row_key_index]
                result._data[row_name] = row
        
        return result
    
    def writeReport(self,report,file_name,delimiter='\t',quotechar="'"):
        
        with open(file_name,'wb') as fh:
            writer = csv.writer(fh, delimiter=delimiter,quotechar=delimiter, quoting=csv.QUOTE_MINIMAL)
            header = report.getColumnNames()
            writer.writerow(header)
            for row_name in report.getRowNames():
                row = report._data[row_name]
                writer.writerow(row)
    
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
    output_file = 'SafeReport_out_test.txt'
    
    reader = DataReportIO()
    read_start_time = time()
    report = reader.readReport(input_file,'row_name',delimiter='\t')
    read_time = time() - read_start_time
    print "read time for %s is [%s]" % (input_file,read_time)
    
    column_names = report.getColumnNames()
    row_names = report.getRowNames()
    print "Reloading the data into a new Report size rows [%s] x columns [%s]" % (len(row_names),len(column_names))
    
    start_fill_time = time()
    new_report = SafeReport()
    row_index = 0
    for r_name in row_names:
        #print"row index [%s]" % (row_index)
        row_index += 1
        
        for c_name in column_names:
            #print "column name [%s]" % (c_name)
            value = report.get(r_name,c_name)
            if value != None:
                new_report.add(r_name,c_name,value)
    end_fill_time = time() - start_fill_time 
    print "Reloading time [%s]" % (end_fill_time) 
    
    print "Writing report [%s]" % (output_file)
    start_write_time = time()
    reader.writeReport(report, output_file, delimiter="\t")
    end_write_time = time() - start_write_time
    print "Writing time [%s]" % (end_write_time)
    
    print "done"