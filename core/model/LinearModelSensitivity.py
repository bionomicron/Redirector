'''
Created on Oct 26, 2010
@author: Graham Rockwell
'''
import pickle
import os.path
from core.model.LinearModel import LinearModel
from core.model.LPSolver import LPSolver
from core.util.Report import Report
from core.util.ReportWriter import ReportWriter

from itertools import combinations
from numpy import array, arange, mean


def targetValueCombinations(targets,values):
    result = []
    for target in targets:
        for value in values:
            result.append((target,value))
    return result

def combinationLimitBuilder(model,targetValues,controlValueName):
    if "testControl" in model.getRowNames():
        model.removeRow("testControl")
    for (target,value) in targetValues:
        model.addData("testControl",target,value)
    model.addData("testControl",controlValueName,-1.0)
    model.addRowLimit("testControl",(0.0,0.0))
    model.addColumnLimit("controlValueName",(None,None))
    return model

def addCombinationLimit(model,targetValues,limit):
    '''
    @type targetValues: [(string,float)] list of doubles for reactions and coefficients to be make into a single row in linear optimization model
    '''
    rowTag = ''
    for (target,value) in targetValues:
        rowTag += "(%s:%s)" % (target,value)
    for (target,value) in targetValues:
        model.addData(rowTag,target,value)
    model.addRowLimit(rowTag,limit)
    return model

def enzymeBoundaryControl(model,targets,bounds,objectiveName,productionName,searchSize):
    targets = bounds.keys()
    controlNames = model.getControlsForNames(targets)
    eControlMap = model.getEnzymeControlMap()
    
    controlSubSets = combinations(controlNames,searchSize)
    xcontrolSubSets = set(controlSubSets)
    
    result = Report()
    iter = 0
    for icontrolNames in xcontrolSubSets:
        ibounds = {}
        controlTag = ''
        ienzymeNames = model.annotateGeneList(list(icontrolNames))
        ienzymeNames.sort()
        for ieName in ienzymeNames:
            controlTag += "(%s)" % (ieName)
        for icontrolName in icontrolNames:
            itargets = eControlMap[icontrolName]
            for itarget in itargets:
                if itarget in bounds.keys():
                    ibound = bounds[itarget]
                    ibounds[itarget] = ibound
            #print ("control %s => [%s]") % (icontrolName,ibounds)
            pass
        if len(ibounds) == 0:
            continue
        (fluxMap,oflux,pflux) = findBoundaryProduction(model, ibounds, None, objectiveName, productionName)
        print "[%s] objective %s => production %s" % (controlTag, oflux, pflux)
        iter += 1
        #iEnzymeNames = model.annotateGeneList(icontrolNames) 
        #for icontrolName in iEnzymeNames:
        #    result.add(iter, icontrolName, "active")
        result.add(controlTag,"natural",oflux)
        result.add(controlTag,"production",pflux)
    
    return result
        

def findBoundaryProduction(originalModel, bounds, combiLimits, objectiveName, productionName):
    model = LinearModel()
    model.extend(originalModel)
    if combiLimits != None:
        for (targetValues,limit) in combiLimits:
            model = addCombinationLimit(model, targetValues,limit)
    #objectiveMap = {objectiveName:-1.0}
    objectiveMap = {productionName:1.0}
    model.addColumnLimits(bounds)
    lp = LPSolver()
    fluxMap = lp.run(model,objectiveMap)
    objectiveValue = fluxMap[objectiveName]
    productionValue = fluxMap[productionName]
    return (fluxMap,objectiveValue,productionValue)


def LinearModelVariableBoundarys(originalModel, objectiveName, targets=None, pickleFileName=None, strict=False, minObjectivePercent=None,searchSize=1):
    '''
    Find min and max boundaries for the dependent variables of a linear model. 
    Uses limitMap as constraints on model.
     
    Add cumulative boundary persistence as pickle file function.
    
    @param linearModel: linear model to be analyzed
    @type linearModel: LinearModel
    @param limitMap: map of variable names and lower and upper limits
    @type limitMap: {variableName, (float, float)  
    @param targets: list of variables to perform analysis on
    @type targets: list
    @return: dictionary of lower and upper limits for each variable
    @rtype: {variableName: (lowerlimit, upperlimit)}
    '''
    
    linearModel = LinearModel()
    linearModel.extend(originalModel)
    
    if pickleFileName != None:
        if os.path.isfile(pickleFileName):
            print "loading saved flux boundaries"
            pFile = open(pickleFileName)
            result = pickle.load(pFile)
            pFile.close()
            return result
    
    originalLimits = linearModel.getColumnLimits()
    originalObjectiveLimit = originalLimits[objectiveName]
    objectiveVector = {objectiveName:-1.0}
    
    solver = LPSolver()
    originalValues = solver.run(model=linearModel, objective=objectiveVector)
    originalObjectiveValue = originalValues[objectiveName]
        
    if minObjectivePercent != None:
        minObjectiveValue = originalObjectiveValue * minObjectivePercent
        linearModel.addColumnLimit(objectiveName,(minObjectiveValue,None))
    
    if targets == None:
        targets = linearModel.getColumnNames()
    targets.add(objectiveName)

    print "Searching flux boundaries for %s > %s" % (objectiveName, minObjectiveValue)
    result = {}
    targetValues = targetValueCombinations(targets, [-1.0,1.0])
    controlSubSets = combinations(targetValues,searchSize)
    xcontrolSubSets = set(controlSubSets)
    
    for name in targets:
        nameTag = name
    #   if name not in linearModel.getColumnNames():
    #       print "target %s not found in model"
    #       print "boundary discovery not possible"
    #       continue
        
    #for iTargetValues in xcontrolSubSets:
    #    name = "testControlValue"
    #    nameTag = iTargetValues
    #    linearModel = combinationLimitBuilder(linearModel, iTargetValues, name)
        
        result[name] = (None,None)
        negObjective = {name:1.0}
        posObjective = {name:-1.0}
        
        solver.clearObjective()
        negValues = solver.run(linearModel, objective=negObjective)
        negValue = negValues[name]
        
        solver.clearObjective()
        posValues = solver.run(linearModel, objective=posObjective)
        posValue = posValues[name]
        
        delta = 1e-4
        if posValue != 0 or negValue !=0:
            pass
        
        if False:
            originalValue = originalValues[name]
            originalLimit = originalLimits[name]
        
            if not (negValue - delta <= originalValue <= posValue + delta):
                print "original value not in boundary"
                print "objective value [%s] (%s) %s <= %s <= %s (%s)" % (name, originalLimit[0], negValue, originalValue, posValue, originalLimit[1])
                #print "value neg code %s pos code %s" % (nscode, pscode)
                (negValue,posValue) = (originalLimit[0],originalLimit[1])
            else:
                #print "objective value [%s] (%s) %s <= %s <= %s (%s)" % (name, originalLimit[0], negValue, originalValue, posValue, originalLimit[1])
                pass
                
            if negValue == None or negValue == float("-inf"):
                continue
    
            if posValue == None or posValue == float("-inf"):
                continue
        
        result[nameTag] = (negValue,posValue)
        print "%s < [%s] < %s " % (negValue,nameTag,posValue)
        pass    
    
    if pickleFileName != None:
        pFile = open(pickleFileName,'w')
        pickle.dump(result, pFile)
        pFile.close()
        
    solver.clear()
    del solver
    
    linearModel.addColumnLimit(objectiveName,originalObjectiveLimit)
    #targets.remove(objectiveName)
        
    return result                 
    
class LinearModelSensitivity():
    '''
    Analysis object for linear model sensitivity analysis
    * search the range of independent variables for analysis of sensitivity of the objectives
    * takes in limit map for bounding of search
    * 
    '''
    
    def __init__(self):
        self.verbose = False
        self.minPerturb = 1e-2 # should be greater than 2x delta
        self.delta = 1e-3 # should be smaller than 1/2 minPerturb
        self.deltaRange = 0.1
        self.perturbDivisions = 3
        self.objectiveRange = 0.5
        self.sensivityMaxRange = 100
        self.searchIterMax = 30
        self.sensitivityMin = -200
        self.sensitivityMax = 200
        self.lowerMax = -200
        self.upperMax = 200
        self.maxSearch = 100
        self.fluxBounds = None
        self.solver = LPSolver
        
    def _checkLimit(self,limit):
        (lower,upper) = limit
        if lower == None:
            lower = float("-inf")
        if upper == None:
            upper = float("inf")
        return (lower,upper)
    
    def _magnatueLimits(self,values,limits):
        '''
        @param values: map of names and value for variables in lp problem 
        @type value: map
        @return: returns map of limits which constrain the magnitude of the values to their current level
        @rtype: {string:(float,float)}
        '''
        result = {}
        for (name,value) in values.items():
            (lower,upper) = (0,0)
            if name in limits.keys():
                (ilower,iupper) = limits[name]
                lower = max([lower,ilower])
                upper = min([upper,iupper])
            
            if value == 0:
                result[name] = (0,0)
            elif cmp(value,0) == 1:
                result[name] = (lower,value)
            elif cmp(value,0) == -1:
                result[name] = (value,upper)
                
        return result
    
    def _getAdjust(self,value,limit):
        '''
        Find divisions of range for adjustment values
        @return: list of adjustment to search for the variable
        @rtype: [adjustmentValues]
        '''
        if self.perturbDivisions == None:
            result = arange(-self.deltaRange,self.deltaRange+self.delta,self.delta)*value
            return result
        negInc = 0
        posInc = 0
        if limit[0] != float("-inf"):
            negInc = (value - limit[0])/self.perturbDivisions
        if limit[1] != float("inf"):
            posInc = (limit[1] - value)/self.perturbDivisions
        resultNeg = []
        resultPos = []
        if negInc >  self.delta: #! consider changing to less restrictive
            resultNeg = arange(limit[0] - value, 0, negInc)
        if posInc > self.delta:
            resultPos = arange(0,limit[1] - value, posInc)
        result = []
        result.extend(resultNeg)
        result.extend(resultPos)
        return result
    
    
    def _perturbVar(self, solver, name, value, varLimit, pertLimit):
        '''
        @rtype: {variableValue:sensitivity}
        '''
        result = {}
        
        originalVarLimit = varLimit
        
        varLimit = self._checkLimit(varLimit)
        pertLimit =self._checkLimit(pertLimit)
        adjustLower = max(varLimit[0],pertLimit[0])
        adjustUpper = min(varLimit[1],pertLimit[1])
        adjustLimit = (adjustLower,adjustUpper)
        
        adjust = self._getAdjust(value,adjustLimit)
        originalValues = solver.run()
        primeObjValue = solver.getObjectiveValue()
        
        if pertLimit[1] - pertLimit[0] < self.minPerturb:
            return (None,None)
        
        for a in adjust:
            a = round(a,6)
            l = value + a
            if abs(a) < self.minPerturb:
                continue 
            solver.addColumnLimit(name, l, l)
            newValues = solver.run()
            objValue =  solver.getObjectiveValue()
            
            deltaObjValue = -1.0 * (objValue - primeObjValue)/a
            deltaObjValue = round(deltaObjValue,6)
            #print "[%s = %s] %s = %s - %s / %s" % (name,l,deltaObjValue,objValue,primeObjValue,a)
            result[a] = deltaObjValue
                 
        oVarLimits = {name:originalVarLimit}    
        solver.setColumnLimitsByNames(oVarLimits)

        return result
    
    def _filterSenseValue(self,value):
        if value < self.sensitivityMin:
            return None
        if value > self.sensitivityMax:
            return None
        
    
    def _findSensitvityValue(self,model,objective,name,limit,primeObjValue,value):
        
        originalLimit = model.getColumnLimits()[name]
        model.addColumnLimit(name,(limit,limit))
        
        solver = self.solver()
        solver.verbose = False
        solver.setModel(model)
        solver.clearObjective()
        solver.setObjectiveMap(objective)
        try:
            newValues = solver.run() 
        except:
            if self.verbose: print "Solver crash: failed to solve [%s] sensitivity [%s] = %s" %(objective,name,limit)
            del solver
            model.addColumnLimit(name,originalLimit)
            return None
        
        objValue =  solver.getObjectiveValue()
        solver.clear()
        del solver
        model.addColumnLimit(name,originalLimit)
        
        if objValue == None:
            print"failed to solve [%s] sensitivity [%s] = %s" %(objective,name,limit)
            print"This could represent a serious problem in the model or solver being used"
            solver2 = self.solver()
            model.addColumnLimit(name,(limit,limit))
            solver2.setModel(model)
            scode = solver2.runSimplex()
            prediction = solver2.getPredictionMap()
            tValue = prediction[name]
            objValue = solver2.getObjectiveValue()
            solver2.clear()
            del solver2
            model.addColumnLimit(name,originalLimit)
            print "2nd Solver attempt %s" % objValue
            print "Solver Code %s" % scode
            if scode != 200:
                return None
        else:
            if self.verbose: print"solved [%s] sensitivity [%s] = %s" %(objective,name,limit)
            
        deltaValue = limit - value
        deltaObjValue = -1.0 * (objValue - primeObjValue)
        if deltaValue == 0:
            return None
        deltaObjSense = deltaObjValue / deltaValue 
        deltaObjSense = round(deltaObjSense,6)
        
        return deltaObjSense
            
    def _mergeZone(self,value,range,revQue,searchQue):
        keyValues = array(revQue.keys())
        overlap = (keyValues > value - self.delta) & (keyValues < value + self.delta)
        if sum(overlap) > 0:
            keys = keyValues[overlap]
            (lower,upper) = range
            for key in keys:
                (ilower,iupper) = revQue[key]
                lower = min(lower,ilower)
                upper = max(upper,iupper)
                del revQue[key]
            revQue[mean(keys)] = (lower,upper)
            range = (lower,upper)
        else:
            revQue[value] = range
        
        return (revQue,range)
            
            
        
    def _perturbVarAuto(self, model, objective, name, pertLimit):
        '''
        @rtype: {variableValue:sensitivity}
        '''
        #Re-making solver for every optimization test appears to be critical!
        #Failures here result in invalid optimizations
        
        iterMax = self.searchIterMax
        varLimit = model.getColumnLimits()[name]
        
        solver = self.solver()
        solver.setModel(model)
        solver.clearObjective()
        solver.setObjectiveMap(objective)
        originalValues = solver.run()
        primeObjValue = solver.getObjectiveValue()
        value = originalValues[name]
        solver.clear()
        del solver

        varLimit = self._checkLimit(varLimit)
        pertLimit =self._checkLimit(pertLimit)
        adjustLower = max(varLimit[0],pertLimit[0])
        adjustUpper = min(varLimit[1],pertLimit[1])
        adjustLimit = (adjustLower,adjustUpper)
        
        if pertLimit[0] < varLimit[0] - self.delta or pertLimit[1] > varLimit[1] + self.delta:
            print "Sensitivity search, perturb for %s limits %s larger than natural limit %s" % (name,pertLimit,varLimit) 

        if abs(adjustLimit[1] - adjustLimit[0]) < self.minPerturb:
            return None
        
        if self.maxSearch != None:
            if originalValues[name]- adjustLimit[0] > self.maxSearch:
                print "search space [%s] (%s,%s) too large , limiting search" % (name, adjustLimit[0],adjustLimit[1])
                lower = originalValues[name] - self.maxSearch
                upper = adjustLimit[1]
                adjustLimit = (lower,upper)
            
            if adjustLimit[1] - originalValues[name] > self.maxSearch:
                print "search space [%s] (%s,%s) too large , limiting search" % (name, adjustLimit[0],adjustLimit[1])
                lower = adjustLimit[0]
                upper = originalValues[name] + self.maxSearch
                adjustLimit = (lower,upper)
                    
        searchQue = [(adjustLimit[0],adjustLimit[1])]
        revQue = {}
        resultQue = {}
        resultPointQue = {}
        iter = 0

        while len(searchQue) > 0 and iter < iterMax:
            iter += 1
            (lower,upper) = searchQue.pop(0)
            if self.verbose: print "[%s] sensitivity scan %s <-> %s" % (name,lower,upper)
            
            #This may be a very important step, search to small a space or exceeding the possible solutions may cause very slow solving.
            if abs(upper-lower) < max(self.minPerturb,2*self.delta):    
                if self.verbose: "sensitivity search space to small %s = [%s,%s]" % (name,lower,upper)
                continue 
            if abs(value - lower) < self.minPerturb:
                lower += self.delta
                pass   
            if abs(upper - value) < self.minPerturb:
                upper -= self.delta
                pass
                
            if self.verbose: print "sensitivity: %s: [%s] = %s " % (name,lower,value)
            if lower in resultPointQue.keys():
                lowValue = resultPointQue[lower]
            else:
                lowValue  = self._findSensitvityValue(model,objective,name,lower,primeObjValue,value)
                if lowValue == None:
                    print "no lower sensitivity"
                    lowValue = 0
                lowValue = round(lowValue,4)
                resultPointQue[lower] = lowValue
            if self.verbose: print "lower sensitivity = %s" % (lowValue)
            
            if self.verbose: print "sensitivity: %s: [%s] = %s " % (name,upper,value) #! currently being messed with to check for non solves
            if upper in resultPointQue.keys():
                highValue = resultPointQue[upper]
            else:
                highValue  = self._findSensitvityValue(model,objective,name,upper,primeObjValue,value)
                if highValue == None:
                    print "no upper sensitivity"
                    highValue = 0
                highValue = round(highValue,4)
                resultPointQue[upper] = highValue
            if self.verbose: print "upper sensitivity = %s" % (highValue)
            
            if abs(highValue - lowValue) < self.delta:
                nValue = mean([highValue,lowValue])
                resultQue[(lower,upper)] = nValue
                self._mergeZone(nValue, (lower,upper), revQue, searchQue)
            else:
                (revQue,lrange) = self._mergeZone(lowValue, (lower,lower), revQue, searchQue)
                (revQue,urange) = self._mergeZone(highValue, (upper,upper), revQue, searchQue)
                
                (llower,lupper) = lrange
                (ulower,uupper) = urange
        
                if self.verbose: print "%s -[%s]- %s vs %s -[%s]- %s" %(llower,lowValue,lupper,ulower,highValue,uupper)
                
                if abs(lupper - ulower/2) > self.delta:
                    searchQue.append((lupper,(lupper + ulower)/2 - self.delta))
                if abs(ulower/2+self.delta - ulower) > self.delta: 
                    searchQue.append(((lupper + ulower)/2+self.delta,ulower))
                
        sMap = {}
        for (key,svalue) in revQue.items():
            sMap[svalue] = key
        sKeys = sMap.keys()
        sKeys.sort()
        
        result = sMap
        
        return result
        
    def perturb(self, model, pertubationLimits, objectiveName, targets = None):
        '''
        @return: map of variables to map of flux values and sensitivities
        @rtype: {variableName:{varibleValue:sensitivity}}
        '''
        
        if targets == None: targets = model.getColumnNames()
        autoResult = {}
        varLimits = model.getColumnLimits()
        objective = {objectiveName:-1.0}
        
        for name in targets:
            
            if name not in model.getColumnNames():
                print "Sensitivity testing %s not found in model" % (name)
                continue
            
            pertLimit = pertubationLimits[name]

            #adjustCoeff = self._perturbVar(solver, name, value, varLimit, pertLimit)
            #adjustResult[name] = adjustCoeff
            
            if self.verbose: print "Finding sensitivity coefficients for [%s] [%s]" % (objectiveName,name)
            autoCoeff = self._perturbVarAuto(model, objective, name, pertLimit)
            
            currentLimits = model.getColumnLimits()
            if currentLimits != varLimits:
                print "Sensitivity analysis is altering model limits" 
                model.setColumnLimits(varLimits)
            
            if autoCoeff != None:
                autoResult[name] = autoCoeff
            if self.verbose: print "complete"
            
        return autoResult
    
    def sensitivityAnalysisSampling(self, model, tObjVal, tObjectiveName, targets = None, id=None):
        '''
        @var tPredVal:
        @var tObjVal:
        @var tObjectiveName:
        '''
        sensitivityLibrary = []
        if id == None:
            id = tObjectiveName
        prefix = "_%s_" % id
        
        bPercent = self.objectiveRange
        naturalObjLimits = {tObjectiveName:(tObjVal*(1-bPercent),tObjVal*(1+bPercent))}
        originalLimits = model.getColumnLimits()
        
        #magnitudeLimits = self._magnatueLimits(tPredVal,originalLimits)
        if self.verbose: print "Finding flux boundaries"
        if self.fluxBounds == None:
            fluxBoundaries = LinearModelVariableBoundarys(model, tObjectiveName, targets)
        else: fluxBoundaries = self.fluxBounds
        
        if self.verbose: print "Searching flux space for sensitivity values"
        senseValues = self.perturb(model, fluxBoundaries, tObjectiveName, targets)
        if self.verbose: print "Sensitivity search complete"
        
        #model.addColumnLimits(magnitudeLimits)
        #senseValuesConst = self.perturb(model, tPredVal, fluxBoundaries, tObjectiveName, targets)
        
        model.addColumnLimits(originalLimits)
        
        if self.verbose: print "Processing sensitivity values into control library"
        senseValuesReduced = self.processSensitivity(senseValues)
        #senseValuesConstReduced = self.processSensitivity(senseValuesConst)
        
        sensitivityLibrary.append(("Sensitivity Reduced %s" % id, "Sense_R_" + prefix, senseValuesReduced))
        
        return(fluxBoundaries, sensitivityLibrary, senseValues)
    
    def unifyLibraryValues(self,controlLibrary,factor):
        '''
        unify control library to maximum absolute negative and positive values
        '''
        negControl = {}
        posControl = {}
        
        for (name,prefix,map) in controlLibrary:
            for (key,values) in map.items():
                if values == None or len(values) == 0:
                    continue
                v1 = min(values) * factor
                v2 = max(values) * factor
                
                if v1 < 0 and v1 <= v2:
                    negControl[key] = v1
                if v2 < 0 and v2 < v1:
                    negControl[key] = v2
                if v1 > 0 and v1 >= v2:
                    posControl[key] = v1
                if v2 > 0 and v2 > v1:
                    posControl[key] = v2
        
        return (negControl,posControl)
    
    def multiplyControlLibrary(self,controlLibrary,value):
        result = []
        for (name, prefix, control) in controlLibrary:
            icontrol = {}
            for (key,v) in icontrol.items():
                icontrol[key] = v*value
            result.append((name,prefix,control))
        return result
    
    def writeSensitivityLibraryReport(self, senseLibrary, fileName, report=None):
        if report == None:
            report = Report()
        for (libraryName,lib) in senseLibrary.items():
            report.addColumnHash(libraryName, lib)
        
        writer = ReportWriter()
        writer.setFile(fileName)
        code = writer.write(report)
        writer.closeFile()
        
        return code
        
    
    def processSensitivityCorelation(self,adjustmentMap):
        delta = 0.1
        keys = array(adjustmentMap.keys())
        values = array(adjustmentMap.values())
        #values = values[argsort(keys)]
        result = []
        previous = 0
        for v in values:
            if v == 0:
                if previous != 0:
                    previous = 0
                    result.append(v)
            elif (v - previous)/v > delta:
                previous = v
                result.append(v)
        return result
                 
            
    def processSensitivity(self,sensitivityMap):
        '''
        @var sensitivityMap:  
        @type sensitivityMap:
        @rtype:  
        '''
        result = {}
        
        for (key, map) in sensitivityMap.items():
            if map == None:
                continue
            sValues = self.processSensitivityCorelation(map)
            result[key] = sValues
    
        return result