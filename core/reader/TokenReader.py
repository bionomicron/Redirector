'''
@author: Graham Rockwell
@organization: Church lab
@note: currently broken
@version: out of date
'''

import re

class TokenParser:
    
    def __init__(self,delimiter = "\t",listDelimiter = ' ',tokenBreak = "///"):
        self.delimiter = delimiter
        self.listDelimtier = listDelimiter
        self.handle = None
        self.tokenBreak = tokenBreak
        self.emptyToken = ''
        self.currentToken = None
        self.currentValue = None
        self.currentEntry = None
        self.tokenReg = re.compile("^([^" +listDelimiter + "]*)"+ listDelimiter + "*(.*)")
        
    def addToken(self,token):
        self.tokens.append(token)
    
    def setFile(self,fileName):
        self.handle = open(fileName, 'r' )
        
    def close(self):
        self.handle.close()
        
    def parseLine(self,state,line):
        line = line.rstrip()
        if line == self.tokenBreak:
            self.currentEntry[self.currentToken] = self.currentValue
            return "done"
        
        m = self.tokenReg.match(line)
        
        if m == None:
            raise "failed to parse line"
        (tokenId,values) = m.groups()
        
        values = values.strip()
        values = values.strip(';')
        
        if tokenId != '':
            self.currentToken = tokenId
            self.currentValue = values
            state = "new"
        else:
            self.currentValue = self.currentEntry[self.currentToken]
            self.currentValue = self.currentValue + "\n" + values
            state = "open"
            
        self.currentEntry[self.currentToken] = self.currentValue
        
        return state
            
    def getEntry(self):
        state = "start"
        self.currentEntry = {}
        while state != "done":
            line = self.handle.readline()
            #end of file
            if line == '':
                return None
            #get new state
            state = self.parseLine(state,line)
            
        return self.currentEntry
    

from core.util.DatabaseConnection import DatabaseConnection
from core.reader.ConfigParser import GeneralSimpleConfigParser

class LoadToken:
    
    def __init__(self):
        self.connection = None
        self.cursor = None
        
    def configure(self,config):
        connectionFac = DatabaseConnection()
        connectionFac.config(config)
        self.connection = connectionFac.getConnection()
        
    def open(self):
        self.cursor = self.connection.cursor()
        
    def close(self):
        self.cursor.close()
        
    def insertEntry(self,entry):
        code = entry['ENTRY'].strip().split()[0]
        name = entry['NAME'].split("\n")[0]
        sql = "insert ignore into keg (entry,name) values (\"%s\",\"%s\")" % (code,name)

        try:
            self.cursor.execute(sql)
        except Exception,msg:
            print "failed to load" + msg
            raise Exception
    
    
if __name__ == '__main__':
    databaseConfigParser = GeneralSimpleConfigParser()
    configFileName = "database_config.txt"
    configDb = databaseConfigParser.parse(configFileName)
    loader = LoadToken()
    loader.configure(configDb)
    
    
    testfile="testdata/kegg/compound"
    parser = TokenParser()

    parser.setFile(testfile)
    entry = parser.getEntry()
    
    loader.open()
    while entry != None:
        #do something
        try:
            entryCode = entry['ENTRY'].split()[0]
            name = entry['NAME'].split('\n')[0]
            loader.insertEntry(entry)
            print "[%s] %s" % (entryCode,name)
        except Exception,msg:
            print msg
        
        entry = parser.getEntry()
        
    parser.close()
    loader.close()
            
    
    