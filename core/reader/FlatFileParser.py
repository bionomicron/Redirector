#!/usr/bin/env python
'''
@author: Graham Rockwell
@organization: Church Lab Harvard Genetics
@version: 02/18/2013
'''

import string, sets, re
from core.do.FluxModel import FluxModel, MetabolicNetwork, Reaction, Objective, FluxLimit
from core.do.FluxModelTools import ReactionParser
from core.util.Report import Report

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
        
    def throwCheckError( self ):
        self.failLine = True

    def passCheckError( self ):
        self.failLine = False
        
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

    def addRequiredHeader( self, dataName, header ):
        self.requiredHeaderMap[header] = dataName

    def addOptionalHeader( self, dataName, header ):
        self.optionalHeaderMap[header] = dataName

    def addRequiredHeaders( self, headerMap ):
        self.requiredHeaderMap.update( headerMap )
        
    def addOptionalHeaders( self, headerMap ):
        self.optionalHeaderMap.update( headerMap )

    def setRequiredHeaderMap( self, headerMap ):
        self.resetHeaders()
        self.addRequiredHeaders( headerMap )

    def setHeader( self, headers ):
        self.resetHeaders()
        for h in headers:
            self.addRequiredHeader( h, h )

    def checkRequiredHeaders( self, line ):
        '''
        Checks the header line to see if required headers are present.
        '''
        rheaders = sets.Set( self.requiredHeaderMap.keys() )
        sLine = sets.Set( line )

        if rheaders.issubset( sLine ):
            return True
        else:
            headerErrorTag = "Expecting headers:[%s]\n found:[%s] \n expected - found (missing): [%s]" % ( '|,|'.join( rheaders ), '|,|'.join( sLine ), '|,|'.join( rheaders - sLine )) 
            raise ParseError(headerErrorTag)
        
        return False

    def indexHeader( self, line ):
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
        result = {}
        for header in self.requiredHeaderMap.keys():
            dataName = self.requiredHeaderMap[header]
            value = line[self.headerIndex[header]]
            result[dataName] = self.checkFieldValue( value )
        return result

    def parseOptional( self, line ):
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

    def safeSplit( self, line ):
        '''
        Splits the line into a list
        checks for a wrapper on the elements of the list and removes them
        @return: [String]
        '''
        line.replace( "\n", '' )
        nline = []
        for v in line.split(self.delimiter):
            v = string.replace(v,self.wrapperString,'')
            v = string.strip(v)
            nline.append(v)
        return nline
    
    def checkHeader( self, headerLine ):

        hLine = self.safeSplit( headerLine )
        hLineSet = sets.Set( hLine )
        
        #print hLine

        if len( hLineSet ) != len( hLine ):
            raise ParseError( "Duplicate column name %s" % ( hLine ) )

        if self.checkRequiredHeaders( hLine ):
            self.indexHeader( hLine )        
        
        return True

    
    def parseLine( self, line ):
        '''
        
        '''
        result = TagedElement()
        
        sline = self.splitLineCheck( line )

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
            return strip_line[0] == self.comment
        else:
            #print "Comment: [%s]" % strip_line
            return True
        
    def splitLineCheck( self, line ):
        
        sline = self.safeSplit( line )
        
        if len( sline ) != len( self.headerIndex ):
            if ( self.failLine ):
                print self.headerIndex
                raise ParseError( "Line should have %d column found %s \n [%s] \n" % ( len( self.headerIndex ), len( sline ) , sline ) )
            else:
                return None
        else:
            return sline
        
    def getNextLine(self):
        try:
            line = self.fileHandle.next()
        except StopIteration:
            line = None
        return line
        
    def startFile(self,fileName):
        self.fileHandle = open(fileName,'r')
        self.index = 0
        line = self.getNextLine()
        while line != None:
            if self.isComment( line ):
                self.index += 1
            elif self.index < self.startLine:
                self.index += 1
            elif self.index == self.headerLine:
                self.checkHeader( line )
                self.index += 1
                return True
            elif self.index > self.startLine:
                value = self.parseLine( line )
                self.index += 1
                return False
            line = self.getNextLine()
        return False     

    def getLine(self):
        result = None
        line = self.getNextLine()
        if line == None:
            result = None
        elif self.isComment( line ):
            result = ""
        elif self.index > self.startLine:
            value = self.parseLine( line )
            result = value
        self.index += 1
        return result
        
    def closeFile(self):
        self.fileHandle.close()

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
                value = self.parseLine( line )
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
                sLine = self.safeSplit(line)
                result.append(sLine)
            index += 1
        return result
            
    
    def parseGenericHeader(self, headerLine, unique = True):
        hLine = self.safeSplit( headerLine )
        hLineSet = sets.Set( hLine )

        if len( hLineSet ) != len( hLine ) and unique:
            raise ParseError( "Duplicate column name %s" % ( hLine ) )

        if self.setHeader(hLine ):
            self.indexHeader( hLine )        
        
        return hLine

    def xparseGenericReport(self, fileName, key=None, unique = True):
        '''
        checking for removal
        '''
        result = Report()
        header = None
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
                header = self.parseGenericHeader( line, unique )
                if key in header:
                    kIndex = header.index(key)
                    
            elif self.endLine > index > self.startLine:
                sLine =self.safeSplit(line)
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
    
    def parseToMap( self, fileName, keyTag, valueTag = None, multi=False ):
        result = {}
        self.startFile(fileName)
        d = self.getLine()
        
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
            
            d = self.getLine()
            
        self.closeFile()
        return result
        
    def parseAnyToMap(self,header,fileName,keyTag,valueTags=None):
        self.setHeader(header)
        result = self.parseToMap(fileName,keyTag,valueTags)
        return result
    
    def parseToReport( self, fileName, keyTag, header = None):
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
        d = self.getLine()
        
        if header == None:
            self.parseGenericHeader(headerLine=d, unique=True)
            d = self.getLine()
        
        while d != None:
            if d != "":
                keyName = d[keyTag]
                del d[keyTag]
                
                for valueTag in d.keys():
                    v = d[valueTag]
                    result.addElement(keyName,valueTag,v)

            d = self.getLine()
            
        self.closeFile()
        return result
        

class SmallMoleculeFlatFileParser( FlatFileParser ):
    '''
    @summary: Parser for compound / small molecule flat files
    '''
    
    def __init__(self,delimiter='\t',comment='#'):
        FlatFileParser.__init__( self, delimiter, comment )
        requiredHeaders = {"Abbreviation":"Name", "OfficialName":"OfficialName"}
        optionalHeaders = {"concentration":"concentration","lower":"lower","upper":"upper","compartment":"compartment"}
        
        self.addRequiredHeaders( requiredHeaders )
        self.addOptionalHeaders( optionalHeaders )
    
    def parse( self, file ):
        result = {} #currently reaction are parsed in an unordered form.
        self.startFile(file)
        d = self.getLine()
        
        while d != None:
            if d != "":
                m= {}
                name = d["Name"] 
                if d["compartment"] != '':
                    name = name + "[" + d["compartment"] + "]"
                m["name"] = name
                m["officialName"] = d["OfficialName"]
            for term in self.optionalHeaderMap.keys():
                if term in d.keys():
                    m[term] = d[term]
            d = self.getLine()
            result[name] = m
            
        self.closeFile()
        return result


class ReactionFlatFileParser( FlatFileParser ):
    '''
    @summary: Parser for reaction flat files
    '''

    def __init__( self, delimiter='\t', comment='#' ):
        FlatFileParser.__init__( self, delimiter, comment )
        requiredHeaders = {"Abbreviation":"Name", "OfficialName":"OfficialName", "Equation":"Equation"}
        optionalHeaders = {}

        self.addRequiredHeaders( requiredHeaders )
        self.addOptionalHeaders( optionalHeaders )
    
    def parse( self, rxnfile ):
        rxnParser = ReactionParser()
        result = {} #currently reaction are parsed in an unordered form.
        self.startFile(rxnfile)
        d = self.getLine()
        
        while d != None:
            if d != "":
                reaction = rxnParser.parseEquation(d["Equation"])
                reaction.setId(d["Name"])
                reaction.setName(d["OfficialName"])
                reaction.setEquation(d["Equation"])
                reaction.annotation = d.annotation
                
                result[d["Name"]] = reaction
                
            d = self.getLine()
            
        self.closeFile()
        return result

class MetabolicNetworkFlatFileParser( FlatFileParser ):
    '''
    @summary: parser for metabolic network files
    @note: currently only reaction file matters
    '''

    def parse( self, rxnfile, metabolitefile = ''):
        precursor_predicate = {}
        rxnParser = ReactionFlatFileParser()
        smmParser = SmallMoleculeFlatFileParser()
        #smm = smmParser.parse( metabolitefile )
        rxn = rxnParser.parse( rxnfile )
        result = MetabolicNetwork( 'MetabolicNetwork' )
        
        result.addReactionMap( rxn )
        return result



class  FluxLimitFlatFileParser( FlatFileParser ):
    '''
    @summary: parser for reaction limit flat file
    '''

    def __init__( self, delimiter='\t', comment='#' ):
        FlatFileParser.__init__( self, delimiter, comment )
        self.setHeader( ["Flux", "Lower", "Upper"] )
        self.noLimit = ['None']
        self.posInf = 1e3
        self.negInf = -1e3
    
    def parse( self, fluxlimitfile ):
        fluxLimit = FluxLimit()
        lines = open( fluxlimitfile, 'r' )

        isHeader = True

        for line in lines:
            if not self.isComment( line ):
                if isHeader:
                    self.checkHeader( line )
                    isHeader = False
                else:
                    d = self.parseLine( line )
                    lower = d["Lower"]
                    upper = d["Upper"]
                    
                    if lower in self.noLimit:
                        lower = None
                    elif self.negInf != None:
                        if float(lower) < self.negInf:
                            lower = None
                
                    if upper in self.noLimit:
                        upper = None
                    elif self.posInf != None:
                        if float(upper) > self.posInf:
                            upper = None
                    
                    fluxLimit[d["Flux"]] = ( lower, upper )
        return fluxLimit

class ObjectiveFlatFileParser( FlatFileParser ):
    '''
    @summary: Parser for objective flat file
    '''
    
    def __init__( self, delimiter='\t', comment='#' ):
        FlatFileParser.__init__( self, delimiter, comment )
        self.setHeader( ["Flux", "Coefficient"] )
    
    def parse( self, objectivefile ):
        objective = Objective()
        lines = open( objectivefile, 'r' )

        isHeader = True

        for line in lines:
            if not self.isComment( line ):
                if isHeader:
                    self.checkHeader( line )
                    isHeader = False
                else:
                    d = self.parseLine( line )
                    objective[d["Flux"]] = d["Coefficient"]

        return objective



#===========================================================
#==========FluxAnalysisFlatFileParser=======================
#===========================================================
class FluxAnalysisFlatFileParser( FlatFileParser ):

    def __init__( self, delimiter='\t', comment='#' ):
        FlatFileParser.__init__( self, delimiter, comment )
        self.reactionParser = ReactionParser()
        
        requiredHeaders = {
                            "Abbreviation":"name", 
                            "Equation":"equation", 
                            "Lower Bound":"lowerBound", 
                            "Upper Bound":"upperBound"
                            }
        optionalHeaders = {"objective":"objective"}
        
        self.addRequiredHeaders( requiredHeaders )
        self.addOptionalHeaders( optionalHeaders )
        self.setObjectiveTag( "objective function" )

    def setObjectiveTag( self, objectiveTag ):
        self.objectiveTag = objectiveTag

    def parseReactionLine( self, line ):
        name = line["name"]
        equation = line["equation"]
        reaction = self.reactionParser.parse( name, "", "", equation )
        return reaction

    def parseLimitLine( self, line ):
        lower = line["lowerBound"]
        upper = line["upperBound"]
        return ( lower, upper )

    def parseObjectiveLine( self, line ):
        return -1

    def parse( self, analysisFile ):
        lines = self.parseFile( analysisFile )

        objective = Objective()
        limit = FluxLimit()
        reactions = {}

        for line in lines:
            name = line["name"]
            reaction = Reaction()

            reaction = self.parseReactionLine( line )
            reactions[name] = reaction
            
            
            if line["objective"] == self.objectiveTag:
                objectiveValue = self.parseObjectiveLine( line )
                objective[name] = objectiveValue

            else:
                limitValue = self.parseLimitLine( line )
                if reaction.direction == 'IRREVERSIBLE-LEFT-TO-RIGHT':
                    ( lower, upper ) = limitValue
                    limitValue = ( 0.0, upper )
                limit[name] = limitValue

        return ( reactions, objective, limit )
                
        
        
        

#============================================================
#================FluxModelFlatFileParser=====================
#============================================================
class  FluxModelFlatFileParser( FlatFileParser ):
    
    def __init__(self):
        self.directionCheck = True
        self.posInf = 1e5
        self.negInf = -1e5
        self.networkParser = MetabolicNetworkFlatFileParser()
        self.exchangeParser = FluxAnalysisFlatFileParser()
        self.fluxLimitParser = FluxLimitFlatFileParser()
        self.objectiveParser = ObjectiveFlatFileParser()
        self.smmParser = SmallMoleculeFlatFileParser()
        self.annotationParser = AnnotationFlatFileParser()
        
        self.modelName = ''
        self.basedir = ''
        self.rxnfile = ''
        self.fluxlimitfile = ''
        self.smmfile = ''
        self.objectivefile = ''
        self.geneannotationfile = ''
    
    def checkDirectionality(self,network,fluxLimits):
        for reactionName in network.getOrderedReactionNames():
            reaction = network.getReaction(reactionName)
            if reactionName in fluxLimits.keys():
                limit = fluxLimits[reactionName]
            else:
                #limit = (self.negInf,self.posInf)
                limit = (None,None)

            (lower,upper) = limit
            direction = reaction.getDirection()
            
            if reaction.getDirection() == 'IRREVERSIBLE-LEFT-TO-RIGHT':
                if lower == None:
                    limit = (0.0,upper)
                elif lower < 0.0:
                    limit = (0.0,upper)
                
            if reaction.getDirection() == 'IRREVERSIBLE-RIGHT-TO-LEFT':
                if upper == None:
                    limit = (lower,0.0)
                elif upper > 0.0:
                    limit = (lower,0.0)
                
            fluxLimits[reactionName] = limit
            
        return fluxLimits

    def parse( self, rxnfile, metabolitefile, fluxlimitfile, objectivefile, analysisName, exchangeFile = None):
        #networkParser = MetabolicNetworkFlatFileParser()
        network = self.networkParser.parse( rxnfile, metabolitefile )
        fluxModel = FluxModel( network )
        
        if fluxlimitfile != '' and fluxlimitfile != None:
            fluxLimits= self.fluxLimitParser.parse( fluxlimitfile )
        else:
            fluxLimits = FluxLimit()
            
        if objectivefile != '' and objectivefile != None:
            objective= self.objectiveParser.parse( objectivefile )
        else:
            objective = Objective()
            
        if exchangeFile != '' and exchangeFile != None:
            ( reactions, objective, limit ) = self.exchangeParser.parse(exchangeFile)
            network.addReactionMap(reactions)
            fluxLimits.extend(limit)
            
        if self.directionCheck:
            fluxLimits = self.checkDirectionality(network, fluxLimits)
            
        fluxModel.addAnalysis( analysisName, objective, fluxLimits )
        
        return fluxModel
    
    def generateModel(self,modelName=None):
        if modelName == None:
            modelName = self.modelName
        dir = self.basedir
        
        network = None
        fluxLimits = None
        objective = None
        metaboliteData = None

        if self.rxnfile != '':
            rxnFile =dir + self.rxnfile
            network = self.networkParser.parse(rxnFile)
        else:
            network = MetabolicNetwork('')
        
        if self.smmfile != '':
            smmFile = dir + self.smmfile
            metaboliteData = self.smmParser.parse(smmFile)
        
        if self.fluxlimitfile != '':
            fluxLimitFile = dir + self.fluxlimitfile
            fluxLimits = self.fluxLimitParser.parse(fluxLimitFile)
        else:
            fluxLimits = FluxLimit()
            
        if self.objectivefile != '':
            objectiveFile = dir + self.objectivefile
            objective = self.objectiveParser.parse(objectiveFile)
        else:
            objective = Objective()
                        
        if self.directionCheck:
            fluxLimits = self.checkDirectionality(network, fluxLimits)
        
        #Add data to collector
        if network != None:
            result = FluxModel(network)
            
        if metaboliteData != None:
            result.setMetaboliteData(metaboliteData)
            
        if objective != None and fluxLimits != None:
            result.addAnalysis( modelName, objective, fluxLimits )
            
        return result
    
    def parseGenConfig(self,config,modelName):
        dir = config["basedir"]
        network = None
        fluxLimits = None
        smmData = None
        objective = None
        metaboliteData = None

        if config["rxnfile"] != '':
            rxnFile =dir + config["rxnfile"]
            network = self.networkParser.parse(rxnFile)
        else:
            network = MetabolicNetwork('')
        
        if  config["smmfile"] != '':
            smmFile = dir + config["smmfile"]
            metaboliteData = self.smmParser.parse(smmFile)
        
        if config["fluxlimitfile"] != '':
            fluxLimitFile = dir + config["fluxlimitfile"]
            fluxLimits = self.fluxLimitParser.parse(fluxLimitFile)
        else:
            fluxLimits = FluxLimit()
            
        if config["objectivefile"] != '':
            objectiveFile = dir + config["objectivefile"]
            objective = self.objectiveParser.parse(objectiveFile)
        else:
            objective = Objective()
            
        if config["exchangefile"] != '':
            exchangeFile = dir + config["exchangefile"]
            ( reactions, objective, limit ) = self.exchangeParser.parse(exchangeFile)
            network.addReactionMap(reactions)
            fluxLimits.extend(limit)
            
        if self.directionCheck:
            fluxLimits = self.checkDirectionality(network, fluxLimits)
        
        #Add data to collector
        if network != None:
            result = FluxModel(network)
            
        if metaboliteData != None:
            result.setMetaboliteData(metaboliteData)
            
        if objective != None and fluxLimits != None:
            result.addAnalysis( modelName, objective, fluxLimits )

        return result
    
    def parseConfig(self,config,modelName):
        rxnFile =config.getRxnFile()
        smmFile = config.getSmmFile()
        limitFile = config.getFluxLimitFile()
        objectiveFile = config.getObjectiveFile()
        exchangeFile = config.getExchangeFile()
        
        fluxModel = self.parse(rxnFile,smmFile,limitFile,objectiveFile,modelName,exchangeFile)
        
        return fluxModel


    def parseMultipleAnalysis( self, rxnfile, metabolitefile, fluxlimitfiles, objectivefiles, analysisNames ):
        networkParser = MetabolicNetworkFlatFileParser()
        network = networkParser.parse( rxnfile, metabolitefile )

        fluxModel = FluxModel( network )
        
        fluxLimits = {}
        objective = {}
        analysis_names = []

        fluxLimitParser = FluxLimitFlatFileParser()
        objectiveParser = ObjectiveFlatFileParser()

        for analysisName in analysisNames:
            fluxlimitfile = fluxlimitfiles[analysisName]
            objectivefile = objectivefiles[analysisName]
            
            fluxLimits= fluxLimitParser.parse( fluxlimitfile )
            objective= objectiveParser.parse( objectivefile )

            fluxModel.addAnalysis( analysisName, objective, fluxLimits )

        return fluxModel

    def parseSuperCross( self, rxnfile, metabolitefile, fluxlimitfiles, objectivefiles ):
        analysisNames = []
        newFluxLimitFileList = {}
        newObjectiveFileList = {}
        
        for objectivefile in objectivefiles:
            for limitfile in fluxlimitfiles:
                name = objectivefile + "-" + limitfile
                analysisNames.append( name )
                newFluxLimitFileList[name] = limitfile
                newObjectiveFileList[name] = objectivefile

        fluxModel = self.parseMultipleAnalyis( self, rxnfile, metabolitefile, newFluxLimitFileList, newObjectiveFileList )
        return fluxModel


class AnnotationFlatFileParser( FlatFileParser ):

    def __init__( self, delimiter='\t', comment='#' ):
        FlatFileParser.__init__( self, delimiter, comment )
        requiredHeaders = {"bNumber":"Name", "gene":"gene", "function":"function"}
        optionalHeaders = {}

        self.addRequiredHeaders( requiredHeaders )
        self.addOptionalHeaders( optionalHeaders )
    
    def parse( self, file ):
        result = {}
        self.startFile(file)
        d = self.getLine()
        
        while d != None:
            if d != "":
                result[d["Name"]] = [d["gene"],d["function"]]
            d = self.getLine()
            
        self.closeFile()
        return result


#=====================================
#=============ParseError==============
#=====================================
class ParseError( Exception ):
    pass
