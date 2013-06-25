#!/usr/bin/env python
'''
@author: Graham Rockwell
@organization: Church Lab Harvard Genetics
@version: 03/04/2013

--Rename to metabolic-model parser
'''
from core.do.FluxModel import FluxModel, MetabolicNetwork, Reaction, Objective, FluxLimit
from core.do.FluxModelTools import ReactionParser
from core.reader.FlatFileParser import FlatFileParser

class TagedElement(dict):
    
    def __init__(self):
        self.annotation = {}
        
    def __str__(self):
        result = dict.__str__(self)
        result += self.annotation.__str__()

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
        d = self.getTagedLine()
        
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
            d = self.getTagedLine()
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
        d = self.getTagedLine()
        
        while d != None:
            if d != "":
                reaction = rxnParser.parseEquation(d["Equation"])
                reaction.setId(d["Name"])
                reaction.setName(d["OfficialName"])
                reaction.setEquation(d["Equation"])
                reaction.annotation = d.annotation
                
                result[d["Name"]] = reaction
                
            d = self.getTagedLine()
            
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
                    d = self.parseTagedLine( line )
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
                    d = self.parseTagedLine( line )
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
        d = self.getTagedLine()
        
        while d != None:
            if d != "":
                result[d["Name"]] = [d["gene"],d["function"]]
            d = self.getTagedLine()
            
        self.closeFile()
        return result

#=====================================
#=============ParseError==============
#=====================================
class ParseError( Exception ):
    pass
