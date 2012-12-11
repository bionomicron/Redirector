'''
Created on Dec 21, 2010
@author: Graham Rockwell
@note: Borken
'''

from core.model.LinearModel import LinearModel
from core.reader.FlatFileParser import FlatFileParser
import re


class ParseEquations:
    
    def __init__(self):
        self.identity = None
        self.nullTerm = "Nothing"
        self.defaultLimit = (None,None)
        
        self.overRideIdentityLimit = False
        self.overRideLimit = False
        self.setDefault = True
        
    def parseCoefficient( self, element, direction = "+"):
        coefficient = 1.0
        regex = re.compile( '([0-9\\.E-]+)' )
        m = regex.match( element )
        if m:
            coefficient = float( m.group( 1 ) )
            element = m.string[m.end( 1 ):]
            element = element.strip()
            #print "%s: %f" %(participant, coefficient)
        return ( element, float( coefficient ) )

    def parseEquation( self, equation):
        listOfParticipants = []
        equation.replace( '"', '' )
        if equation.find( ' + ' ) != -1:
            for element in equation.split( ' + ' ):
                p = element.strip()
                value = self.parseCoefficient(p)
                listOfParticipants.append( value )
        else:
            element = equation.strip()
            if element != '' and ( element.find( self.nullTerm ) == -1 ):
                p = element.strip()
                value = self.parseCoefficient(p)
                listOfParticipants.append( value )
        return listOfParticipants
    
    def parseInequality(self,equation):
        regex = re.compile('(.*)([<=>])(.*)')
        m = regex.match(equation)
        if m:
            equation = m.group(1)
            relation = m.group(2)
            limit = m.group(3)
        else:
            equation = equation
            relation = None
            limit = None
        return (equation,relation,limit)
    
    def addData(self,model,rowName,eqn):
        values = self.parseEquation(eqn)
        for (variableName,coeffecent) in values:
            model.addData(rowName,variableName,coeffecent)
        return model
    
    def addRowLimit(self,model,name,relation,limit):
        if relation == "=":
            model.addRowLimit(name,(limit,limit))
        if relation == "<":
            model.addRowLimit(name,(None,limit))
        if relation == ">":
            model.addRowLimit(name,(limit,None))
        if relation == None:
            return None
    
    def addInequality(self,model,name,eqn):
        (equation,relation,limit) = self.parseInequality(eqn)
        self.addData(model,name,equation)
        self.addRowLimit(model,name,relation,limit)
        return model
    
    def setVariableLimits(self,model,limit):
        for name in model.getColumnNames():
            model.addColumnLimit(name,limit)
        return model
    
    def parseObjective(self,eqn):
        result = {}
        values = self.parseEquation(eqn)
        for (variableName,coeffecent) in values:
            result[variableName] = coeffecent
        return result
    
    def applyTerms(self,model,rowName,terms,prefix,target):
        for (tag,coeffecent) in terms:
            if tag == self.identity:
                itag = target
            else:
                itag = prefix % tag
            
            model.addData(rowName,itag,coeffecent)
            
            if (tag == self.identity) & (not self.overRideIdentityLimit):
                continue
            elif self.overRideLimit or (self.setDefault & (itag not in model.getColumnLimits())):
                model.addColumnLimit(itag,self.defaultLimit)
            
        return model
    
    def parseApplication(self,equation):
        regex = re.compile('(.*):(.*)([<=>])(.*)')
        m = regex.match(equation)
        rowTag = None
        relation = None
        limit = None
        
        if m:
            rowTag = m.group(1)
            equation = m.group(2)
            relation = m.group(3)
            limit = m.group(4)
            return (rowTag,equation,relation,limit)
        
        regex = re.compile('(.*)([<=>])(.*)')
        m = regex.match(equation)
        rowTag = None
        if m:
            equation = m.group(1)
            relation = m.group(2)
            limit = m.group(3)

        return (rowTag,equation,relation,limit)
    
    
    def applyRelation(self,model,relationName,value,targets=[''],prefixPattern="%s"):
        '''
        variable pattern
        row name => <relationName><target>
        variable name => <prefixPattern%tag><target>
        '''
        (rowTag,equation,relation,limit) = self.parseApplication(value)
        terms = self.parseEquation(equation)
        for target in targets:
            if rowTag == None:
                rowName = relationName + target
            else:
                rowName = rowTag % target
            prefix = prefixPattern + target
            self.applyTerms(model,rowName,terms,prefix,target)
            self.addRowLimit(model,rowName,relation,limit)
            
        return model

class EquationParser:
    '''
    Parse metabolic reaction equation(s) 
    generate linear model
    '''
    
    def __init__(self):
        self.defaultLimit = (None,None)
        self.identity = "identity"        
        self.setDefault = True
        self.overRideIdentityLimit = False
        self.overRideLimit = False
         
    def parseCoefficient( self, element, direction = "+"):
        coefficient = 1.0
        regex = re.compile( '([0-9\\.E-]+)' )
        m = regex.match( element )
        if m:
            coefficient = float( m.group( 1 ) )
            element = m.string[m.end( 1 ):]
            element = element.strip()
            #print "%s: %f" %(participant, coefficient)
        return ( element, float( coefficient ) )

    def parseEquation( self, equation):
        listOfParticipants = []
        equation.replace( '"', '' )
        if equation.find( ' + ' ) != -1:
            for element in equation.split( ' + ' ):
                p = element.strip()
                value = self.parseCoefficient(p)
                listOfParticipants.append( value )
        else:
            element = equation.strip()
            if element != '' and ( element.find( 'Nothing' ) == -1 ):
                p = element.strip()
                value = self.parseCoefficient(p)
                listOfParticipants.append( value )
        return listOfParticipants
    
    def parseInequality(self,equation):
        regex = re.compile('(.*)([<=>])(.*)')
        m = regex.match(equation)
        if m:
            equation = m.group(1)
            relation = m.group(2)
            limit = m.group(3)
        else:
            equation = equation
            relation = None
            limit = None
        return (equation,relation,limit)
    
    def addData(self,model,rowName,eqn):
        values = self.parseEquation(eqn)
        for (variableName,coeffecent) in values:
            model.addData(rowName,variableName,coeffecent)
        return model
    
    def addRowLimit(self,model,name,relation,limit):
        if relation == "=":
            model.addRowLimit(name,(limit,limit))
        if relation == "<":
            model.addRowLimit(name,(None,limit))
        if relation == ">":
            model.addRowLimit(name,(limit,None))
        if relation == None:
            return None
    
    def addInequality(self,model,name,eqn):
        (equation,relation,limit) = self.parseInequality(eqn)
        self.addData(model,name,equation)
        self.addRowLimit(model,name,relation,limit)
        return model
    
    def setVariableLimits(self,model,limit):
        for name in model.getColumnNames():
            model.addColumnLimit(name,limit)
        return model
    
    def parseObjective(self,eqn):
        result = {}
        values = self.parseEquation(eqn)
        for (variableName,coeffecent) in values:
            result[variableName] = coeffecent
        return result
    
    def applyTerms(self,model,rowName,terms,prefix,target):
        for (tag,coeffecent) in terms:
            if tag == self.identity:
                itag = target
            else:
                itag = prefix % tag
            
            model.addData(rowName,itag,coeffecent)
            
            if (tag == self.identity) & (not self.overRideIdentityLimit):
                continue
            elif self.overRideLimit or (self.setDefault & (itag not in model.getColumnLimits())):
                model.addColumnLimit(itag,self.defaultLimit)
            
        return model
    
    def applyRelation(self,model,relationName,value,targets=[''],prefixPattern="%s"):
        '''
        variable pattern
        row name => <relationName><target>
        variable name => <prefixPattern%tag><target>
        '''
        (equation,relation,limit) = self.parseInequality(value)
        terms = self.parseEquation(equation)
        for target in targets:
            rowName = relationName + target
            prefix = prefixPattern + target
            self.applyTerms(model,rowName,terms,prefix,target)
            self.addRowLimit(model,rowName,relation,limit)
            
        return model
    
    def addMipColumns(self,model,names,tag="%s"):
        for name in names:
            iname = tag % name
            model.setMipColumnName(iname)
        return model


class ReactionEquationParser:

    def __init__( self ):
        '''
        Parses reaction string to reaction object
        '''
        free = (float("-inf"),float("inf"))
        forward = (0.0, float("inf"))
        reverse =  (float("-inf"),0.0)
        self.directionals = {'<==>': free,'<-->': free,'-->': forward, "<=>": free,'<--': reverse}

    def parseCoefficient( self, participant):
        coefficient = 1.0
        name = participant
        m = re.compile( '[ (]?([0-9\\.E-]+)[) ]?[ ]+([^ ]*)' ).match( participant)
        #m = re.compile( '\(([0-9\\.E-]+)\)' ).match( participant )
        if m:
            coefficient = float( m.group( 1 ) )
            if coefficient == 0:
                raise "zero coeffecent when parsing element: %s" % (participant)
            name = m.group(2)
            #participant = m.string[m.end( 1 )+1:]
            name = name.strip()
            #print "%s: %f" %(name, coefficient)
        return (name , float( coefficient ) )

    def parseCompartment( self, participant, compartment ):
        m = re.compile( '(\[[a-z]\])' ).search( participant )
        if m:
            participant = m.string[:m.start( 1 )]
            compartment = m.group( 1 )
            compartment = compartment.strip()
        elif compartment == None:
            compartment = '[c]'

        participant = participant.strip()
        #print "%s to %s" % (participant, physEnt)
        return ( participant, compartment)

    def parseParticipant( self, participant, compartment ):
        ( physEnt, compartment ) = self.parseCompartment( participant, compartment )
        ( physEnt, coefficient ) = self.parseCoefficient( physEnt )
        return ( physEnt, coefficient, compartment )

    def parseParticipants( self, participants, compartment ):
        listOfParticipants = []
        participants.replace( '"', '' )
        if participants.find( ' + ' ) != -1:
            for participant in participants.split( ' + ' ):
                listOfParticipants.append( self.parseParticipant( participant.strip(), compartment ) )
        else:
            participant = participants.strip()
            if participant != '' and ( participant.find( 'Nothing' ) == -1 ):
                listOfParticipants.append( self.parseParticipant( participant.strip(), compartment ) )
        return listOfParticipants

    def parseDirection( self, equation ):
        
        for direction in self.directionals.keys():
            if equation.find( direction ) != -1:
                try:
                    ( left, right ) = equation.split( direction )
                except:
                    print "Failed to parse equation %s" % equation
                return ( left, right, self.directionals[direction] )
        raise Exception("(No direction found) Unable to parse equation: %s " % equation)
    
    def parseEquationPrefix( self, equation ):
        m = re.compile( '^(\[[a-z]\])\s*:?\s*(.*)' ).match( equation )
        prefix = None

        if m:
            prefix = m.group( 1 )
            equation = m.group( 2 )

        return ( equation, prefix )

    def parseConversionType( self, equation):
        proteinClass = None
        m = re.compile( '(\[[a-z]\])\s*:\s*(.*)' ).match( equation )
        if m:
            prefix = m.group( 1 )
            equation = m.group( 2 )
            return ( equation, 'TypeUnspecified' )
        
        return ( equation, 'transport' )
    
    def parseEquation( self, equation, name = ''):
        result = {}
        
        ( newEquation, compartment ) = self.parseEquationPrefix( equation )
        ( newEquation, conversionType ) = self.parseConversionType( newEquation )
        ( leftString, rightString, direction ) = self.parseDirection( newEquation )
        
        left = self.parseParticipants( leftString, compartment )
        right = self.parseParticipants( rightString, compartment )
        
        for (item,value,location) in left:
            tag = item + location
            result[tag] = value * -1.0
            
        for (item,value,location) in right:
            tag = item + location
            result[tag] = value * 1.0
        
        return (result,direction)


class ReactionFileParser( FlatFileParser ):
    '''
    @summary: Parser for directly parsing metabolic model flat file to linear model object 
    '''

    def __init__(self, delimiter="\t", comment="#"):
        FlatFileParser.__init__( self, delimiter, comment )
        requiredHeaders = {"Abbreviation":"Name", "OfficialName":"OfficialName", "Equation":"Equation"}
        optionalHeaders = {}

        self.addRequiredHeaders( requiredHeaders )
        self.addOptionalHeaders( optionalHeaders )
                 
    def parse( self, file ):
        model = LinearModel()
        rxnParser = ReactionEquationParser()
        self.startFile(file)
        d = self.getLine()
        annotation = {}
        
        while d != None:
            if d != "":
                (reaction,direction) = rxnParser.parseEquation(d["Equation"])
                name = d["Name"]
                
                model.addColumn(name, reaction)
                model.addColumnLimit(name, direction)
                annotation[name] = d.annotation
            
            d = self.getLine()
            
        self.closeFile()
        model.annotation = annotation 
        return model 
    
class LimitFileParser( FlatFileParser ):
    
    def __init__( self, delimiter='\t',comment="#"):
        FlatFileParser.__init__( self, delimiter, comment )
        requiredHeaders = {"Flux":"Name", "Lower":"Lower", "Upper":"Upper"}
        optionalHeaders = {}

        self.addRequiredHeaders( requiredHeaders )
        self.addOptionalHeaders( optionalHeaders )


    def parse(self, file):
        self.startFile(file)
        d = self.getLine()
        result = {}
        
        while d != None:
            if d != "":
                name = d["Name"]
                lower = d["Lower"]
                if lower == "None":
                    lower = None
                else:
                    lower = float(lower)
                upper = d["Upper"]
                if upper == "None":
                    upper = None
                else:
                    upper = float(upper)
                
                result[name] = (lower,upper)
                
            d = self.getLine()
            
        self.closeFile()
        
        return result
    
class ObjectiveFileParser( FlatFileParser ):
    
    def __init__( self, delimiter='\t',comment="#"):
        FlatFileParser.__init__( self, delimiter, comment )
        requiredHeaders = {"Flux":"Name", "Coefficient":"Coefficient"}
        optionalHeaders = {}

        self.addRequiredHeaders( requiredHeaders )
        self.addOptionalHeaders( optionalHeaders )


    def parse(self, file):
        self.startFile(file)
        d = self.getLine()
        result = {}
        
        while d != None:
            if d != "":
                name = d["Name"]
                value = d["Coefficient"]
                result[name] = float(value)
                
            d = self.getLine()
            
        self.closeFile()
        
        return result


class MetabolicNetworkParser():
    
    
    def __init__(self):
        self.directionCheck = True
        self.posInf = 1e5
        self.negInf = -1e5
            
        self.modelName = ''
        self.basedir = ''
        self.rxnfile = ''
        self.fluxlimitfile = ''
        self.smmfile = ''
        self.objectivefile = ''
        self.geneannotationfile = ''
        
        self.rxnParser = ReactionFileParser()
        self.limitParser = LimitFileParser()
        self.objParser = ObjectiveFileParser()
    
    def generateModel(self, modelName=None):
        if modelName == None:
            modelName = self.modelName
        dir = self.basedir
        
        model = None
        fluxLimits = None
        objective = None

        if self.rxnfile != '':
            rxnFile =dir + self.rxnfile
            model = self.rxnParser.parse(rxnFile)
        else:
            model = LinearModel()
        
        if self.fluxlimitfile != '':
            fluxLimitFile = dir + self.fluxlimitfile
            fluxLimits = self.limitParser.parse(fluxLimitFile)
        else:
            fluxLimits = {}
            
        if self.objectivefile != '':
            objectiveFile = dir + self.objectivefile
            objective = self.objParser.parse(objectiveFile)
        else:
            objective = {}
            
        if self.directionCheck:
            #fluxLimits = self.checkDirectionality(network, fluxLimits)
            pass
        
        model.addColumnLimits(fluxLimits)
        
        for (key,value) in objective.items():
            model.addObjective(key, value)
            
        for name in model.getRowNames():
            model.addRowLimit(name, (0.0,0.0))
                        
        return model