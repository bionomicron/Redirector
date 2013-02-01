#!/usr/bin/env python
'''
@author: Graham Rockwell
@organization: Church Lab
Updated 20130131
Set of functions for manipulating linear optimization models of metabolic networks
'''

from core.model.LinearModel import LinearModel
from core.model.LPSolver import LPSolver
from core.util.Report import Report
from numpy import array
import re,operator
    
def vectorMultiply(v1,v2):
    result = 0;
    for key in v1.keys():
        if key in v2.keys():
            value1 = v1[key]
            value2 = v2[key]
            result += value1*value2
            print "%s = %s * %s" % (result,value1,value2)
    return result

def vectorMScalar(v1,s):
    result = {}
    for key in v1.keys():
        value = v1[key]
        result[key] = value * s
    return result

def convertMap(data,type):
    result = {}
    for (key,value) in data.items():
        try:
            ivalue = type(value)
            result[key] = ivalue
        except:
            pass
    return result

def reduceVector(data,delta = 0.1,min=1e-4):
    d = array(data)
    d.sort()
    r = []
    previous = None
    for v in d:
        if previous == None:
            r.append(v)
            previous = v
            continue
        if abs(v) < min:
            continue
        if abs((v-previous)/v) > delta:
            r.append(v)
            previous = v
            continue
    return r

class VarabilityAnalysis():
    
    def __init__(self):
        pass
    
    def reactionBounds(self,model,reactionName,objDir = -1.0):
        solver = LPSolver()
        posObjective = {reactionName:objDir}
        negObjective = {reactionName:1.0*objDir}
        
        posPrediction = solver.run(model,posObjective)
        posValue = posPrediction[reactionName]
        
        negPrediction = solver.run(model,negObjective)
        negValue = negPrediction[reactionName]
        
        return (negValue,posValue)

class ValidateModel():
    
    def __init__(self):
        self.varAnalysis = VarabilityAnalysis()
    
    def exportMetabolite(self,model,metaboliteName,max=1e3):
        solver= LPSolver()
        exportTag = "EX_test_export"
        reaction = {metaboliteName:-1.0}
        objective = {exportTag:-1.0}
        
        model.addColumn(exportTag,reaction)
        model.addColumnLimit(exportTag,(0.0,max))        
        prediction = solver.run(model,objective)
        exportValue = prediction[exportTag]
        
        model.removeColumn(exportTag)
        
        return exportValue

    def traceReaction(self,model,reactionName,threshold):
        (negbound,posbound) = self.varAnalysis.reactionBounds(model, reactionName)
        
        if posbound < threshold:
            mMap = model.getColumnValueMap(reactionName)
            for (key,value) in mMap.items():
                eValue = self.exportMetabolite(model, key)
                print "[%s] (%s) => %s" % (key,value,eValue)
        return None
    
    def testPredictions(self,model,predictions):
        rowPrediction = {}
        rowValues = {}
        for key in model.data.keys():
            (rowName,columnName) = key
            value = model.getData(rowName,columnName)
            if rowName not in rowPrediction.keys():
                rowPrediction[rowName] = 0
                rowValues[rowName] = {}
            preditionsValue = predictions[columnName]
            rowPrediction[rowName] += value * preditionsValue
            rowValues[rowName][columnName] = value * preditionsValue
        return rowPrediction    
    
class ViewLinearModel:
    
    def __init__(self):
        self.delimiter = " "
        
    def selectSubMatrix(self, model, rowFilter = None, columnFilter = None, union = False):
        regexRow = None
        regexColumn = None
        
        if rowFilter != None:
            regexRow = re.compile(rowFilter)
        if columnFilter != None:
            regexColumn = re.compile(columnFilter)
        
        result = LinearModel()
        for key in model.data.keys():
            (rowName,columnName) = key
            value = model.getData(rowName,columnName)
            
            if regexRow: rowCheck = regexRow.match(rowName)
            else: rowCheck = False
            if regexColumn: columnCheck = regexColumn.match(columnName)
            else: columnCheck = False
            
            if union:
                if rowCheck and columnCheck:
                    result.addData(rowName,columnName,value)
            else:
                if rowCheck or columnCheck:
                    result.addData(rowName,columnName,value)
                
        return result
    
    
    def viewModelMatrix(self,modelMatrix):
        dataMap = {}
        
        rowNames = modelMatrix.getRowNames()
        columnNames = modelMatrix.getColumnNames()
        for rowName in rowNames:
            for columnName in columnNames:
                value = modelMatrix.getData(rowName,columnName)
                if value == None:
                    valueString = "x.x"
                else:
                    valueString = value
                print "%s%s" % (valueString, self.delimiter),
            if rowName in modelMatrix.getRowLimits().keys():
                (lower,upper) = modelMatrix.getRowLimit(rowName)
                print "[%s][%s]" % (lower,upper),
            print "\n",
        print"done"

class WriteLinearModelToMPS:
    '''
    Currently not working
    '''
    
    def __init__(self):
        self.name = 'ModelMatrix'
        self.objectiveName = 'objective'
        self.smallDelimiter = ' '
        self.delimiter = '     '
        self.fileName = None
        self.fileHandel = None
        self.integerColumns = []

        self.negInf = -1e3        
        self.posInf =  1e3

    def setName(self,name):
        self.name = name

    def setFile(self,file):
        self.fileName = file
        self.fileHandel = open(file,'w')
        
    def setIntegerColumns(self,integerColumns):
        self.integerColumns = integerColumns

        
    def _writeName(self,name):
        line = "NAME" + self.delimiter + self.delimiter + name + '\n'
        self.fileHandel.write(line)
        

    def _writeHeader(self,value):
        line = value + '\n'
        self.fileHandel.write(line)

    def _writeRowTypeLine(self,stringArray):
        line = self.smallDelimiter
        for value in stringArray:
            line += str(value) + self.delimiter
        line += '\n'
        self.fileHandel.write(line)
        

    def _writeColumnTypeLine(self,stringArray):
        line = ''
        for value in stringArray:
            line += self.delimiter + str(value)
        line += '\n'
        self.fileHandel.write(line)
        
        
    def _writeRows(self,modelMatrix):
        rowNames = modelMatrix.getRowNames()
        self._writeHeader('ROWS')
        self._writeRowTypeLine(['N','objective'])
        for name in rowNames:
            
            (lower,upper)  = modelMatrix.getRowLimit(name)
            hasLower = lower != None
            hasUpper = upper != None
            
            if (not hasLower) and (not hasUpper):
                self._writeRowTypeLine(['N',name])
            elif hasLower and hasUpper and (lower != upper):
                print "unable to write range row %s" % name
            elif lower == upper:
                self._writeRowTypeLine(['E',name])
            elif hasLower and not hasUpper:
                self._writeRowTypeLine(['G',name])
            elif hasUpper and not hasLower:
                self._writeRowTypeLine(['L',name])
            else:
                print "nothing found for %s" % (name)
            
        return None
        
                
    def _writeColumns(self,modelMatrix):
        self._writeHeader('COLUMNS')
        objective = modelMatrix.getObjective()
        integerSparseMatrix = []
        
        sparseMatrix = modelMatrix.getSparseMatrix()
        for name in objective.keys():
            value = objective[name]
            sparseMatrix.append(('objective',name,value))
            
        sparseMatrix.sort(key=operator.itemgetter(1))
        
        for (rowName,columnName,value) in sparseMatrix:
            if columnName in self.integerColumns:
                integerSparseMatrix.append((rowName,columnName,value))
            else:
                lineArray = [columnName,rowName,value]
                self._writeColumnTypeLine(lineArray)
        
        if len(integerSparseMatrix) == 0:
            return None
        
        self._writeColumnTypeLine(['integer_start',"\'MARKER\'",'',"\'INTORG\'"])
        
        for (rowName,columnName,value) in integerSparseMatrix:
            lineArray = [columnName,rowName,value]
            self._writeColumnTypeLine(lineArray)
        
        self._writeColumnTypeLine(['integer_end','\'MARKER\'','','\'INTEND\''])
        
        return None

            
    def _writeRowLimits(self,modelMatrix):
        self._writeHeader('RHS')
        rowNames = modelMatrix.getRowNames()
        
        for name in rowNames:
            
            (lower,upper)  = modelMatrix.getRowLimit(name)
            hasLower = lower != None
            hasUpper = upper != None
            
            if hasLower and hasUpper and (lower != upper):
                print "unable to write range row %s" % name
            if hasLower and hasUpper and (lower == upper):
                self._writeColumnTypeLine(['RHS1', name, lower])
            if hasLower and not hasUpper:
                self._writeColumnTypeLine(['RHS1', name, lower])
            if hasUpper and not hasLower:
                self._writeColumnTypeLine(['RHS1', name, upper])
        return None
                
    def _writeBounds(self,modelMatrix):
        self._writeHeader('BOUNDS')
        columnNames = modelMatrix.getColumnNames()
        for name in columnNames:
            (lower,upper) = modelMatrix.getColumnLimit(name)
            if lower != None:
                line = ['LO','BND1', name, lower]
                self._writeRowTypeLine(line)
            if upper != None:
                line = ['UP','BND1', name, upper]
                self._writeRowTypeLine(line)
            if lower == None:
                line = ['LO','BND1', name, self.negInf]
                self._writeRowTypeLine(line)
            if upper == None:
                line = ['UP','BND1', name, self.posInf]
                self._writeRowTypeLine(line)
                
        return None
        
    def close(self):
        self.fileHandel.close()

    def write(self,name,model):
        self._writeName(name)
        self._writeRows(model)
        self._writeColumns(model)
        self._writeRowLimits(model)
        self._writeBounds(model)
        self._writeHeader('ENDATA')    
        
        
class LpEvaluation:
    
    def __init__(self):
        self.verbose = False
        self.originalTag = ""
        self.defaultLimit = (None,None)
    
    def setOriginalTag(self,tag):
        self.originalTag = tag
        
    def parseCoefficient( self, participant):
        coefficient = 1.0
        name = participant
        m = re.compile( '[ (]?([0-9\\.E-]+)[) ]?[ ]+([^ ]*)' ).match( participant)
        #m = re.compile( '\(([0-9\\.E-]+)\)' ).match( participant )
        if m:
            coefficient = float( m.group( 1 ) )
            if coefficient == 0:
                raise "zero coefficient when parsing element: %s" % (participant)
            name = m.group(2)
            name = name.strip()
        
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
        rxn_equality = ['<','>','=',]
        for direction in rxn_equality:
            if equation.find( direction ) != -1:
                try:
                    ( left, right ) = equation.split( direction )
                except Exception:
                    print "failed to parse"
                return ( left, right, direction )
        raise ParseError, "%s could not be parsed" % equation    

    def expression(self,lp,expression,rowName,columnNames,limit):
        result = lp
        for name in columnNames:
            for term in expression.keys():
                coeffecient = expression[term]
                rowIName = rowName + name
                columnIName = term + name
                if columnIName not in result.getColumnLimits().keys():
                    result.addColumnLimit(columnIName,self.defaultLimit)
                result.addData(rowIName,columnIName,coeffecient)
        
        result.addRowLimit(rowName,limit)
        return result
        