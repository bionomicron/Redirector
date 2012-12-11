#!/usr/bin/env python2.3
'''
@author: Graham Rockwell
@organization: Church Lab Harvard Genetics
@version: 8/22/2009
@note: out of date or with out purpose currently, replaced by biopython tools
Data access object for reading in genetic sequence data
'''

import re        
from cStringIO import StringIO    

class SequenceReader:
    def __init__(self,delimiter = "\t"):
        """
        type delimiter: string
        """
        
        self.delimiter = delimiter
        self.comment = "^>"
        result = {}
    
    def parse(self, fileName):
        """
        type fileName: string
        """
        if fileName == None:
            raise "Null for file name"
        fileHandle = open(fileName, 'r')
        result = {}
        accumulator = StringIO()
        
        #header currently unused
        
        for line in fileHandle:
            if re.search(self.comment,line):
                continue
            
            line = line.rstrip('\n')
            accumulator.write(line)
        
        fileHandle.close()
        result = accumulator.getvalue()
        
        return result