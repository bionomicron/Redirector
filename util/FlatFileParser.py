#!/usr/bin/env python
'''
@author: Graham Rockwell
@organization: Church Lab Harvard Genetics
@version: 03/04/2013

--Rename to metabolic-model parser
'''

import string, sets, re
from util.Report import Report

class TagedElement(dict):
    
    def __init__(self):
        self.annotation = {}
        
    def __str__(self):
        result = dict.__str__(self)
        result += self.annotation.__str__()
        
class FlatFileParser:
    '''
    @summary: Parent class for parsing flat (delimited) files
    '''
    
    def __init__( self, delimiter='\t', comment= '#', emptyField='na' ):
        self.delimiter = delimiter
        self.requiredHeaderMap = {}
        self.optionalHeaderMap = {}
        self.headerIndex = {}
        
        self.headerLine = 0
        self.startLine = 0
        self.endLine = float("inf")

        self.failLine = True
        self.emptyLine = ''
        self.emptyData = ''
        self.wrapperString = '\"'
        self.comment = comment
        self.emptyField = emptyField
        
        self.fileHandle = None
        self.index = None
        
    def setEmptyLine( self, emptyLine ):
        self.emptyLine = emptyLine
        
    def setDelimiter( self, delimiter ):
        self.delimiter = delimiter

    def getDelimiter( self ):
        return self.delimiter

    def setEmptyField( self, emptyField ):
        self.emptyField = emptyField

    def getEmptyField( self ):
        return self.emptyField
    
    def setWrapper(self,wrapperString):
        self.wrapperString = wrapperString

    def setStartLine( self, startLine ):
        self.startLine = startLine

    def setHeaderLine( self, headerLine ):
        self.headerLine = headerLine
        self.startLine = headerLine

    def setDataNames( self, dataNames ):
        self.dataNames = dataNames

    def resetHeaders( self ):
        self.requiredHeaderMap = {}
        self.optionalHeaderMap = {}
        self.headerIndex = {}
        
    def addRequiredHeaders(self, headerMap):
        self.requiredHeaderMap.update(headerMap)
        
    def addOptionalHeaders(self, headerMap):
        self.optionalHeaderMap.update(headerMap)

    def setRequiredHeaderMap( self, headerMap ):
        self.resetHeaders()
        self.addRequiredHeaders( headerMap )

    def setHeader( self, headers ):
        self.resetHeaders()
        headerMap = {}
        for h in headers:
            headerMap[h] = h
        self.addRequiredHeaders(headerMap)

    def _safeSplit( self, line ):
        '''
        @type line: String
        @return: [String]
        @summary:
        Splits the line into a list
        checks for a wrapper on the elements of the list and removes them
        '''
        line.replace('\n','')
        nline = []
        for v in line.split(self.delimiter):
            v = string.replace(v,self.wrapperString,'')
            v = string.strip(v)
            nline.append(v)
        return nline
    
    def _splitLineCheck( self, line ):
        '''
        @type line: String
        @return: [String]
        @summary:
        Splits the line into a list
        checks for a wrapper on the elements of the list and removes them
        '''
        sline = self._safeSplit( line )
        
        if len( sline ) != len( self.headerIndex ):
            if ( self.failLine ):
                print self.headerIndex
                raise ParseError( "Line should have %d column found %s \n [%s] \n" % ( len( self.headerIndex ), len( sline ) , sline ) )
            else:
                return None
        else:
            return sline
        
    def checkRequiredHeaders( self, line ):
        '''
        @var line: list of Strings to be checked for required headers
        @type line: [String] 
        Checks the header line to see if required headers are present.
        @return boolean
        '''
        rheaders = sets.Set( self.requiredHeaderMap.keys() )
        sLine = sets.Set( line )

        if rheaders.issubset( sLine ):
            return True
        else:
            headerErrorTag = "Expecting headers:[%s]\n found:[%s] \n expected - found (missing): [%s]" % ( '|,|'.join( rheaders ), '|,|'.join( sLine ), '|,|'.join( rheaders - sLine )) 
            raise ParseError(headerErrorTag)
        
        return False

    def _indexHeader( self, line ):
        '''
        @summary:
        Find index of each header in the line.
        Used for matching headers to data.
        '''
        index = 0
        for value in line:
            self.headerIndex[value] = index
            index += 1
        return self.headerIndex

    def checkFieldValue( self, value ):
        '''
        Checks the value of the and strip of extra characters 
        '''
        v = string.strip( value )
        if( string.strip( v ) == self.emptyField or string.strip( v ) == '' ):
            return self.emptyData
        else:
            return v

    def parseRequired( self, line ):
        '''
        #! to be removed
        Parses map of list of data matching to required headers
        '''
        result = {}
        for header in self.requiredHeaderMap.keys():
            dataName = self.requiredHeaderMap[header]
            value = line[self.headerIndex[header]]
            result[dataName] = self.checkFieldValue( value )
        return result

    def parseOptional( self, line ):
        '''
        #! to be removed
        Parses map of list of data matching to required headers
        '''
        result = {}
        for header in self.optionalHeaderMap.keys():
            dataName = self.optionalHeaderMap[header]
            if self.headerIndex.has_key( header ):
                index  = self.headerIndex[header]
                value = line[index]
                result[dataName] = self.checkFieldValue( value )
            else:
                result[dataName] = self.emptyData
        return result
    
    def parseAnnotation(self,line):
        '''
        Parses columns that are present but not required in the output.
        #! currently being reviewed for redundancy.
        '''
        result = {}
        knownHeaders = self.requiredHeaderMap.keys()
        knownHeaders.extend(self.optionalHeaderMap.keys())
        for headerName in self.headerIndex.keys():
            if headerName not in knownHeaders:
                index = self.headerIndex[headerName]
                value = line[index]
                result[headerName] = self.checkFieldValue(value)
        return result
    
    def parseTagedLine( self, line ):
        '''
        @var line: 
        @type line: String 
        @return: [String]
        #! being updated to remove annotation and just return data for all headers found.
        '''
        
        result = TagedElement()
        
        sline = self._splitLineCheck( line )

        if len( sline ) != len( self.headerIndex ):
            print self.headerIndex
            raise ParseError( "Line should have %d column found %s \n [%s] \n" % ( len( self.headerIndex ), len( sline ) , sline ) )
        
        requiredValues = self.parseRequired( sline )
        optionalValues = self.parseOptional( sline )
        annotationValues = self.parseAnnotation( sline )

        result.update( requiredValues )
        result.update( optionalValues )
        result.annotation = annotationValues

        return result

    def isComment( self, line ):
        '''
        Checks a line to see if it is a comment line
        '''
        strip_line = string.strip( line )
        if len( strip_line ) > 0 and strip_line != self.emptyLine:
            cprefix = "^%s.*" % self.comment
            m = re.match(cprefix,strip_line)
            check = m != None
            return check
            #return strip_line[0] == self.comment
        else:
            return True
        
    def getNextLine(self):
        try:
            line = self.fileHandle.next()
        except StopIteration:
            line = None
        return line
        
    def getTagedLine(self):
        result = None
        line = self.getNextLine()
        if line == None:
            result = None
        elif self.isComment( line ):
            result = ""
        elif self.index > self.startLine:
            value = self.parseTagedLine( line )
            result = value
        self.index += 1
        return result
        
    def closeFile(self):
        self.fileHandle.close()
        
    def checkHeader( self, headerLine):

        hLine = self._safeSplit( headerLine )
        hLineSet = sets.Set( hLine )

        if len( hLineSet ) != len( hLine ):
            raise ParseError( "Duplicate column name %s" % ( hLine ) )

        if self.checkRequiredHeaders( hLine ):
            self._indexHeader( hLine )
            return True
        
        return False
    
    def parseHeader(self, headerLine, unique = True):
        hLine = self._safeSplit( headerLine )
        hLineSet = sets.Set( hLine )

        if len( hLineSet ) != len( hLine ) and unique:
            raise ParseError( "Duplicate column name %s" % ( hLine ) )

        if self.requiredHeaderMap != {}:
            if not self.checkRequiredHeaders(hLine):
                raise ParseError("Failed to find required headers [%s].\n  In line [%s}" % (hLine))

        self.setHeader(hLine )
        self._indexHeader( hLine )
        
        return hLine

    def parseFile( self, fileName ):
        '''
        @var fileName:
        @type: String
        @summary: 
        Older file parsing function
        @return: [{header:value}] 
        '''
        
        result = []
        lines = open( fileName, 'r' )
        index = 0
        
        for line in lines:
            if self.isComment( line ):
                pass
            elif index < self.startLine:
                break
            elif index == self.headerLine:
                self.checkHeader( line )
            elif index > self.startLine:
                value = self.parseTagedLine( line )
                result.append( value )
            index += 1
            
        lines.close()
        return result
    
    def parseArray(self,fileName,):
        result = []
        lines = open(fileName,'r')
        index = 1
        for line in lines:
            if self.isComment( line ):
                pass
            elif self.endLine >= index >= self.startLine:
                sLine = self._safeSplit(line)
                result.append(sLine)
            index += 1
        return result
    
    def parseToMap( self, fileName, keyTag, valueTag = None, multi=False ):
        result = {}
        self.startFile(fileName)
        d = self.getTagedLine()
        
        while d != None:
            if d != "":
                keyName = d[keyTag]
                del d[keyTag]
                
                if valueTag:
                    value = d[valueTag]
                else:
                    value = d
                
                if multi and valueTag:
                    if keyName not in result.keys():
                        result[keyName] = []
                    result[keyName].append(value)
                else:
                    result[keyName] = value
            
            d = self.getTagedLine()
            
        self.closeFile()
        return result
        
    def parseAnyToMap(self,header,fileName,keyTag,valueTags=None):
        self.setHeader(header)
        result = self.parseToMap(fileName,keyTag,valueTags)
        return result
    
    def startFile(self,fileName):
        '''
        @summary:
        Find header and checks for required headers.
        Indexes headers 
        Finds first line to start parsing delimited file. 
        '''
        
        self.fileHandle = open(fileName,'r')
        self.index = 0
        line = self.getNextLine()
        while line != None:
            if self.isComment( line ):
                self.index += 1
            elif self.index < self.startLine:
                self.index += 1
            elif self.index >= self.headerLine:
                if not self.checkHeader( line ):
                    raise ParseError("failed to find required headers")
                self.index += 1
                return True
            line = self.getNextLine()
        return False
    
    
    def parseGenericReport(self, fileName, keyTag=None, header = None, unique = True):
        '''
        @var fileName: name of flat (delimited) in text format to be parsed
        @type fileName: String
        @var keyTag: ID of column to be used for report key row
        @type keyTag: String
        @summary: 
        Primary Function
        Parses flat file returns report object.
        '''
        
        result = Report()
        kIndex = None
        lines = open( fileName, 'r' )
        index = 0
        
        for line in lines:
            
            if self.isComment( line ):
                pass
            
            elif self.endLine < index < self.startLine:
                index += 1
                continue
            
            elif index == self.headerLine:
                header = self.parseHeader( line, unique )
                if keyTag in header:
                    kIndex = header.index(keyTag)
                    
            elif self.endLine > index > self.startLine:
                line = line.replace('\n','')
                sLine =self._safeSplit(line)
                if kIndex != None:
                    rName = sLine[kIndex]
                else:
                    rName = str(index)
                for i in range(len(sLine)):
                    if i != kIndex:
                        cName = header[i]
                        value = sLine[i]
                        result.add(rName,cName,value)
                    
            index += 1
            
        lines.close()
        return result
    
    def parseToReport( self, fileName, keyTag, header = None, unique = True):
        '''
        @var fileName: name of flat (delimited) in text format to be parsed
        @type fileName: String
        @var keyTag: ID of column to be used for report key row
        @type keyTag: String
        @summary: 
        Primary Function
        Parses flat file returns report object.
        '''
    
        if header != None:
            self.setHeader(header)
                  
        result = Report()
        self.startFile(fileName)
        d = self.getTagedLine()
        
        while d != None:
            if d != "":
                keyName = d[keyTag]
                del d[keyTag]
                
                for valueTag in d.keys():
                    v = d[valueTag]
                    result.addElement(keyName,valueTag,v)

            d = self.getTagedLine()
            
        self.closeFile()
        return result

#=====================================
#=============ParseError==============
#=====================================
class ParseError( Exception ):
    pass
