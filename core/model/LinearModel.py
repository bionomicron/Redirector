#!/usr/bin/env python

import re
from core.util.Cache import IndexedCache, SecondOrderCache
from core.util.Report import Report
from numpy import mean

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

class LinearModel(SecondOrderCache):
    '''
    Core math matrix object for linear optimization modeling in reflux
    Uses the Second Order Cache object as its foundation
    
    Three main parts:
    1)Sparce matrix of coeffecents
    2)Row and column limits
    3)Objective vector
    
    Minimize: 
        Z = c*x
    Subject to: 
        S*x <=> b
        l < x < u
    
    Extended Linear Control Model
    Two additional components:
    4)Annotations
    5)Targets
    6)ControlMap
    6)NaturalObjective
    7)SyntheticObjective
    
    '''
    
    def __init__(self):
        self.verbose = False
        
        #Data Containers
        self.data = SecondOrderCache()
        self.rowLimits = {}
        self.columnLimits = {}
        self.mipColumns = IndexedCache()
        self.objective = {}
        self.modelName = ''
        
        #Report of model annotation values (dict)
        self.annotation = None
        
        #Control
        self.targets = None
        self.controlMap = None
        self.naturalObjective = None
        self.syntheticObjective = None
        
        #Variables
        self.scale = 1 #! to be removed
        self.defaultLowerLimit = None
        self.defaultUpperLimit = None
        
    def __str__(self):
        return self.data.__str__()
    
    def _addString(self,other):
        print "string"
        
    def _addScalar(self,other):
        print "scalar"

    def _addDict(self,other):
        print "dict"
        
    def _addLinearModel(self,other):
        print "model"
    
    
    def __add__(self,other):
        
        if type(other) == type(""):
            self._addString(other)
        if type(other) == type(0) or type(other) ==type(0.0):
            self._addScalar(other)
        if type(other) == type({}):
            self._addDict(other)
        if type(other) == type(self):
            self._addLinearModel(other)
        
        return self
    
    def __eq__(self,value):
        if type(self) != type(value):
            return False
        
        a = self.data.equals(value.data)
        b = self.rowLimits == value.rowLimits
        c = self.columnLimits == value.columnLimits
        d = self.objective == value.objective
        e = self.mipColumns == value.mipColumns
        
        result = a and b and c and d and e
        return result
           
    def _scaleTuple(self,value):
        '''
        scales a par of values
        used in scaling of whole matrix
        @type value: float
        '''  
        (v1,v2) = value
        if v1 != None:
            v1 = float(v1*self.scale)
        if v2 != None:
            v2 = float(v2*self.scale)
        return (v1, v2)


    def _getDefaultColumnLimits(self):
        '''
        returns the column limits with defaults for values not given
        '''
        columnNames = set(self.getColumnNames())
        limitNames = set(self.columnLimits.keys())
        defaultingNames = columnNames.difference(limitNames)
        result = {}
        for name in defaultingNames:
            if not (self.defaultLowerLimit == None and self.defaultUpperLimit == None):
                result[name] = (self.defaultLowerLimit,self.defaultUpperLimit)
        return result
    
    def _getDefaultRowLimits(self):
        '''
        returns the row limits with defaults for values not given
        '''
        columnNames = set(self.getRowNames())
        limitNames = set(self.rowLimits.keys())
        defaultingNames = columnNames.difference(limitNames)
        result = {}
        for name in defaultingNames:
            if not (self.defaultLowerLimit == None and self.defaultUpperLimit == None):
                result[name] = (self.defaultLowerLimit,self.defaultUpperLimit)
        return result
    
    def _annotateString(self,value,annotationMap,regex,nsep=" "):
        
        result = ''
        tags = re.findall(regex,value)
        for s in tags:
            if s in annotationMap.keys():
                r = annotationMap[s]
            else:
                r = s
            result += nsep + r
        result = result[len(nsep):]
        return result
        
    def _annotateMap(self,data,annotationMap,regex):
        result = {}
        for (key,value) in data.items():
            ikey = self._annotateString(key, annotationMap, regex)
            result[ikey] = value
        return result
    
    def _annotateList(self,data,annotationMap,regex):
        result = []
        for value in data:
            ivalue = self._annotateString(value, annotationMap, regex)
            result.append(ivalue)
        return result
    
    
    def annotateGeneList(self,data,annotationName = "bnumber", regex="[a-zA-Z0-9\(\)]+"):
        result = data
        if self.annotation == None:
            return result
        if annotationName in self.annotation.keys():
                annotationMap = self.annotation[annotationName]
                gMap = annotationMap.getColumn("gene")
                if annotationMap != None:
                    result = self._annotateList(data,gMap,regex)
        return result
    

    def annotateGenes(self,objective,annotationName = "bnumber", regex="[a-zA-Z0-9\(\)]+"):
        result = objective
        if self.annotation == None:
            return result
        if annotationName in self.annotation.keys():
                annotationMap = self.annotation[annotationName]
                gMap = annotationMap.getColumn("gene")
                if annotationMap != None:
                    result = self._annotateMap(objective,gMap,regex)
        return result
    
    def getGeneTargetMap(self):
        if self.controlMap == None:
            return None
        result = {}
        for (r,gs) in self.controlMap.items():
            for g in gs:
                if g not in result.keys():
                    result[g] = set()
                result[g].add(r)
        return result
    
    def setProperty(self,name,value):
        self.properties[name] = value

    def getProperty(self,name):
        return self.properties[name]
        
    def setScalex(self,scale):
        '''
        Sets the scale for the matrix
        Not entirely checked for completeness and usage
        @param scale: the scaling factor for the limits
        @type scale: float
        '''
        self.scale  = scale
        
    def addRowName(self,name):
        '''
        Adds a row name to matrix
        @type name: string
        '''
        self.data.rowCache.addValue(name)
        return None
    
    def addColumnName(self,name):
        '''
        Adds column name to matrix
        @type name: string
        '''
        self.data.columnCache.addValue(name)
        return None
    
    def setMipColumnName(self,name):
        '''
        Sets column with value name as a integer column
        @type name: string
        '''
        self.mipColumns.addValue(name)
        
    def setMipColumnNames(self,names,tag="%s"):
        '''
        @type names: string[]
        '''
        for name in names:
            iname = tag % name
            self.setMipColumnName(iname)
        return None
        
    def getMipColumnNames(self):
        '''
        Returns array of strings of column names which are set to integers
        @rtype: string[]
        '''
        return self.mipColumns.getValues()
    
    def getRowIndex(self,name):
        '''
        Returns index of row of value name
        @type name: string
        @rtype: int
        '''
        return self.data.rowCache.getindex(name)
    
    def getColumnIndex(self,name):
        '''
        Returns index of column of value name
        @type name: string
        @rtype: int
        '''
        return self.data.columnCache.getindex(name)
    
    def getRowNames(self):
        '''
        Returns an list of the row names
        @rtype: string[]
        '''
        return self.data.rowCache.getValues()
    
    def getColumnNames(self):
        '''
        Returns a list of the column names
        @rtype: string[]
        '''
        return self.data.columnCache.getValues()
    
    def addRowLimit(self,rowName,limit):
        '''
        Sets limits for selected row
        @type rowName: string
        @type limit (float,float)
        '''
        self.addRowName(rowName)
        self.rowLimits[rowName] = self._scaleTuple(limit)
    
    def addColumnLimit(self,columnName,limit):
        '''
        Sets limits for selected column
        @columnName: string
        @type limit (float,float)
        '''
        self.addColumnName(columnName)
        self.columnLimits[columnName] = self._scaleTuple(limit)
        
    def getRowLimit(self,name):  
        '''
        returns limit for row
        @rtype: (float,float)
        '''
        if name in self.rowLimits:
            return self.rowLimits[name]
        else:
            return (None,None)
        
    def getColumnLimit(self,name):  
        '''
        returns limit for column
        @rtype: (float,float)
        '''
        if name in self.columnLimits:
            return self.columnLimits[name]
        else:
            return (None,None)
        
    def addRowLimits(self,limits):
        '''
        Sets row limits from map
        @type limits: {string:(float,float)}
        '''
        for key in limits.keys():
            (lower,upper) = limits[key]
            self.addRowLimit(key,(lower,upper))
        
    def addColumnLimits(self,limits):
        '''
        Sets colum limits from map
        @type limits: {string:(float,float)}
        '''
        for key in limits.keys():
            (lower,upper) = limits[key]
            self.addColumnLimit(key,(lower,upper))
        
    def getRowLimits(self):
        '''
        returns map of row limits
        @rtype: {name:(float,float)}
        '''
        result = self.rowLimits.copy()
        defaultLimits = self._getDefaultRowLimits()
        result.update(defaultLimits)
        return result
    
    def getColumnLimits(self):
        '''
        returns map of column limits
        @rtype: {name:(float,float)}
        '''
        result = self.columnLimits.copy()
        defaultLimits = self._getDefaultColumnLimits()
        result.update(defaultLimits)
        return result
    
    def addObjective(self,columnName, value):
        '''
        sets objective coeffecent for a column
        @type columnName: string
        @type value: float
        '''
        self.objective[columnName] = value
        
    def setObjective(self,objectiveMap):
        '''
        Sets objective coeffecents from map
        @type objectiveMap: {sting:float}
        '''
        self.objective = {}
        for key in objectiveMap.keys():
            value = objectiveMap[key]
            self.objective[key]=value
    
    def getObjective(self):  
        '''  
        returns objective map
        @rtype {string,float}
        '''
        return self.objective
    
    def getRowValueMap(self,rowName):
        '''
        returns a dict of row values (column name: value)
        @type rowName: string
        @rtype {string,float}
        '''
        r =  self.data.getRow(rowName)
        #return r
        result = {}
        for key in self.data.keys():
            if key[0] == rowName:
                colName = key[1]
                value = self.data.getValue(key[0],key[1])
                result[colName] = value
        if r != result:
            pass
        return result
    
    def getColumnValueMap(self,colName):
        '''
        returns a dict of row values (column name: value)
        @type rowName: string
        @rtype {string,float}
        '''
        r =  self.data.getColumn(colName)
        #return r
        result = {}
        for key in self.data.keys():
            if key[1] == colName:
                rowName = key[0]
                value = self.data.getValue(key[0],key[1])
                #self.getValue(key[0],key[1])
                result[rowName] = value
        if r != result:
            pass
        return result
    
    def getRowByValue(self,name,function):
        result = {}
        values = self.getRowValueMap(name)
        for (k,v) in values.items():
            if function(v):
                result[k] = v
        return result
    
    def getColumnByValue(self,name,function):
        result = {}
        values = self.getColumnValueMap(name)
        for (k,v) in values.items():
            if function(v):
                result[k] = v
        return result
    
    def getRowValuesFromPred(self,name,values):
        result = {}
        rvalues = self.getRowValueMap(name)
        for (k,v) in rvalues.items():
            vi = None
            if k in values:
                vi = values[k]
            result[k] = (v,vi)
        return result 
        
    def addData(self,rowName,columnName,value):
        '''
        Add data value to model matrix
        @type rowName: string
        @type columnName: string
        @type value: string
        '''
        if value == None or value == 0.0:
            return None
        
        floatValue = float(value)
        self.data.addValue(rowName,columnName,floatValue)
        #self.addValue(rowName,columnName,floatValue)
        
        return (0,0)
    
    
    def addDataCache(self,data):
        '''
        Adds a data cach object to the coeffecent matrix
        @type data: SecondOrderCache {(string,string):float}
        '''
        self.data.extend(data)
        #self.data.rowCache.extend(data.rowCache)
        #self.data.columnCache.extend(data.columnCache)
        
    def getData(self,rowName,columnName):
        '''
        Get the coeffecent of the matrix
        @type rowName: string
        @type columnName: string
        '''
        return self.data.getValue(rowName,columnName)
        #return self.getValue(rowName,columnName)
    
    def removeData(self,rowName,columnName):
        '''
        removes a datapoint by row and column name
        @type rowName: string
        @type columnName: string
        '''
        self.data.removeValue(rowName,columnName)
        
    def removeRow(self,rowName):
        columnNames = self.data.rowMap[rowName]
        self.data.removeRow(rowName)
        if rowName in self.rowLimits:
            del self.rowLimits[rowName]
        for columnName in columnNames:
            rCount = len(self.data.columnMap[columnName])
            if rCount == 0:
                self.removeColumn(columnName)
                #if self.verbose: print "removing column [%s]" % (columnName)
                if columnName in self.columnLimits.keys():
                    del self.columnLimits[columnName]
        return None
    
    def removeColumn(self, columnName):
        rowNames = self.data.columnMap[columnName]
        self.data.removeColumn(columnName)
        if columnName in self.columnLimits: 
            del self.columnLimits[columnName]
        if columnName in self.mipColumns.dataArray: 
            self.mipColumns.removeValue(columnName)
        if self.targets != None:
            if columnName in self.targets:
                self.targets.removeValue(columnName)
        if self.controlMap != None:
            if columnName in self.controlMap.keys():
                del self.controlMap[columnName]
        #print "checking rownames [%s]" % (rowNames)
        for rowName in rowNames:
            cNames = self.data.rowMap[rowName]
            if len(cNames) == 0:
                self.removeRow(rowName)
                #if self.verbose: print "removing row [%s]" % (rowName)
                if rowName in self.rowLimits.keys():
                    del self.rowLimits[rowName]
        return None
    
    def addRow(self,rowName,data):
        '''
        add a row in the form of a dictonary
        @type rowName: string
        @type data: {string,float}
        '''
        for key in data.keys():
            value = data[key]
            self.addData(rowName,key,value)
        
    def addColumn(self,columnName,modelVector):
        '''
        @type columnName: string
        @type data: {string:float}
        '''
        for key in modelVector.keys():
            value = modelVector[key]
            self.addData(key,columnName,value)
            
    def getSparseMatrix(self):
        '''
        @rtype (int,int,float)[]
        '''
        return self.data.getSparseMatrix()
    
    def getSparseMatrixMap(self):
        '''
        return the index matrix
        '''
        return self.data.getIndexMatrix()
    
    def addConstraints(self,model):
        '''
        Add a linear constraints to the current model
        @type model: LinearModel
        '''
        self.data.extend(model.data)
        
        rowLimits = model.getRowLimits()
        columnLimits = model.getColumnLimits()
        
        self.addRowLimits(rowLimits)
        self.addColumnLimits(columnLimits)
        
        self.mipColumns.extend(model.mipColumns)
        
        return None
    
    def addModel(self,model):
        '''
        Add a linear model to the current model
        @type model: LinearModel
        '''
        
        self.data.extend(model.data)
        
        rowLimits = model.getRowLimits()
        columnLimits = model.getColumnLimits()
        
        self.addRowLimits(rowLimits)
        self.addColumnLimits(columnLimits)
        
        self.mipColumns.extend(model.mipColumns)
        
        return None
    
    def extend(self,model):
        '''
        Adda linear model and update the objective function
        @type model: LinearModel
        '''
        self.addModel(model)
        self.objective = model.objective
        
    def multiply(self,modelMatrix,value):
        '''
        Multiply the values of the model by a scalar
        @type value: float
        '''
        result = LinearModel()
        
        result.rowCache.extend(modelMatrix.rowCache)
        result.columnCache.extend(modelMatrix.columnCache)
        
        result.data.extend(modelMatrix.data.multiply(value))
        
        rowLimits = modelMatrix.getRowLimits()
        columnLimits = modelMatrix.getColumnLimits()
        
        result.addRowLimits(rowLimits)
        result.addColumnLimits(columnLimits)
        
    def transpose(self):
        '''
        Get the transpose of the current model
        @rtype: LinearModel
        '''
        result = LinearModel()
        result.scale=self.scale
        result.rowCache=self.columnCache
        result.columnCache=self.rowCache
        result.data = self.data.getTranspose()
        return result
    
    def getIndexMatrix(self):
        '''
        Returns the index matrix of the model
        @rtype: (int,int,float)[]
        '''
        return self.data.getIndexMatrix()
    
    def _reportRow(self,dataValues,limit,predValues=None):
        result = ""
        for (k,v) in dataValues.items():
            p = None
            if predValues != None:
                if k in predValues:
                    p = predValues[k]
            r = "%s(%s)[%s] + " %(v,p,k)
            result += r
        result = result[:-2]
        result += " = (%s,%s)" % (limit[0],limit[1])
        return result
            
    def _vectorValue(self,v1,v2):
        result = 0
        for k in set(v1.keys()).intersection(v2.keys()):
            result += v1[k] * v2[k]
        return result
    
    def floatLimit(self,limit):
        if limit[0] != None:
            r1 = limit[0]
        else:
            r1 = float("-inf")
        if limit[1] != None:
            r2 = limit[1]
        else:
            r2 = float("inf")
        result = (r1,r2)
        return result
    
    def modelReport(self,dir=1,prediction=None):
        report = Report()
        delta = 1e-6
        for rowName in self.getRowNames():
            rLimit = self.floatLimit(self.getRowLimit(rowName))
            rValues =self.getRowValueMap(rowName)
            rString = self._reportRow(rValues,rLimit,prediction)
            report.addElement(rowName,"Type","row")
            report.addElement(rowName,"Equation",rString)
            rLimitS = "(%s,%s)" % (rLimit[0],rLimit[1])
            report.addElement(rowName,"Limit",rLimitS)
            if prediction != None:
                rValue = self._vectorValue(rValues, prediction)
                rValue = round(rValue,6)
                rValid = rLimit[0]-delta < rValue < rLimit[1] + delta
                report.addElement(rowName,"Value",rValue)
                report.addElement(rowName,"Valid",rValid)
        for colName in self.getColumnNames():
            cLimit = self.floatLimit(self.getColumnLimit(colName))
            cValues =self.getColumnValueMap(colName)
            cString = self._reportRow(cValues,cLimit,prediction)
            report.addElement(colName,"Type","column")
            report.addElement(colName,"Equation",cString)
            if prediction != None:
                if colName in prediction.keys():
                    cValue = prediction[colName]
                    cValid = cLimit[0]-delta < cValue < cLimit[1] + delta
                    report.addElement(colName,"Value",cValue)
                    report.addElement(colName,"Valid",cValid)
            
        return report
    
    def _convertToGeneNames(self,reactionName):
        if reactionName not in self.controlMap.keys():
            return None
        geneNames = self.controlMap[reactionName]
        if type(geneNames) == type(""):
            return set([geneNames])
        else:
            return geneNames
        
    def _convertToGeneTag(self,reactionName,geneClusters = None):
        geneNames = self._convertToGeneNames(reactionName)
        if geneNames == None:
            return reactionName
        iGeneNames = set(geneNames)
        for geneName in geneNames:
            if geneClusters != None:
                if geneName in geneClusters.keys():
                    for iGeneName in geneClusters[geneName]:
                        iGeneNames.add(iGeneName)
        geneTag = ''
        for iGeneName in iGeneNames:
                geneTag = geneTag + " " + iGeneName
        geneTag = geneTag[1:]
        
        return geneTag
    
    def getEnzymeControlMap(self):
        result = {}
        for (key,values) in self.controlMap.items():
            for v in values:
                if v not in result.keys():
                    result[v] = set() 
                result[v].add(key)
        return result
    
    def getControlsForNames(self,variableNames):
        result = set()
        for name in variableNames:
            geneNames = self._convertToGeneNames(name)
            if geneNames != None:
                for gName in geneNames:
                    result.add(gName)
        return result
    
    def printGeneObjective(self,iObjective,geneClusters=None):
        iGeneObjective = {}
        iOtherObjective = {}
        igControl = {}
        for rxnName in iObjective.keys():
            rControl = iObjective[rxnName]
            if rxnName not in self.controlMap.keys():
                iOtherObjective[rxnName] = rControl
            else:
                geneTag = self._convertToGeneTag(rxnName, geneClusters)
                
                if geneTag not in iGeneObjective:
                    igControl[geneTag] = []
                igControl[geneTag].append(rControl)
        
        for k in igControl.keys():
            iGeneObjective[k] = mean(igControl[k])
            
        return (iGeneObjective,iOtherObjective)
    
    