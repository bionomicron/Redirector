'''
Graham Rockwell
Updated 2012 05 29
Factory for redirector coefficient libraries
'''

from core.model import LinearModel
from core.model.LPSolver import LPSolver
from core.model.LinearModelTools import reduceVector
from core.model.ProcessModel import MatrixTools
from core.model.LinearModelSensitivity import LinearModelSensitivity
from core.util.Report import Report
from core.util.ReportWriter import ReportWriter

from numpy import array
from time import time
from random import Random

import os,pickle


class OptimizationControlFactory:
    
    def __init__(self):
        self.verbose = False
        self.useGeneMap = False
        self.adjustGeneControl = False
        self.compressLibraries = False
        self.usePrices = False
        self.useAdjustedLevels = False
        self.useStepPrices = False
        self.useSensitivity = False
        self.geneMap = None
        self.preload = 0
        self.searchIterMax = 30
        self.maxAbsControl = 100
        self.minAbsControl = 1e-3
    
        self.solverClass = LPSolver
        self.naturalObjectiveName = None
        self.syntheticObjectiveName = None
        
        self.libraryDirectory = "ControlLibraries/"
        self.resultTag = ''
    
        
    def load(self,fileName):
        data = None
        if os.path.exists(fileName):
            if self.verbose: "===loading pre-existing model data=="
            fileHandel = open(fileName,'r')
            data = pickle.load(fileHandel)
            fileHandel.close()
        return data
    
    def filterControlMap(self,control,min=float("-inf"),max=float("inf"),remove = True):
        result = {}
        for (key,value) in control.items():
            if not min < value < max:
                if remove:
                    continue
                elif min > value:
                    value = min
                elif max < value:
                    value = max
            else:
                result[key] = value
        return result        
    
    def multiplySetMap(self,valueMap,value):
        result = {}
        for (key,iValues) in valueMap.items():
            newValues = set()
            for v in iValues:
                newValues.add(v*value)
            result[key] = newValues
        return result 
    
    def fillControlMap(self,controlMap,controlTag,valueMap):
        if controlTag not in controlMap.keys():
            controlMap[controlTag] = {}
        controlSet = controlMap[controlTag]
        for (vkey,value) in valueMap.items():
            if vkey not in controlSet.keys():
                controlSet[vkey] = set()
            controlSet[vkey].add(value)
        return controlMap
    
    def divideMapsDir(self,p,n,m):
        '''
        divides {string:float} map into positive and negative pieces
        '''
        for key in m.keys():
            v = m[key]
            if v > 0.0:
                if key in p.keys():
                    iv = p[key]
                    if v > iv:
                        p[key] = v
                else:
                    p[key] = v
            if v < 0.0:
                if key in n.keys():
                    iv = n[key]
                    if v < iv:
                        n[key] = v
                else:
                    n[key] = v
        return (p,n)
    
    def changeToGeneLibrary(self,library,geneMap=None,dir=0.0):
        '''
        @type geneMap: {recationTag:geneTag}
        
        Flattens library to have the same coefficient for all reactions linked to the same gene control tag
        Changes reaction potentials to be all the same for the same gene group
        Libraries are uni-directional currently
        !Needs to be updated to deal with reaction:set(geneNames) map, be generally better!
        '''
        
        result = {}
        if geneMap == None:
            geneMap = self.geneMap
            
        if geneMap == None:
            return library
    
        if type(geneMap) == type({}):
            reactionGeneMap = {}
            
            for (rxn,genes) in geneMap.items():
                for gene in genes:
                    if gene not in reactionGeneMap.keys():
                        reactionGeneMap[gene] = set()
                    reactionGeneMap[gene].add(rxn)
        
        if type(geneMap) == type([]):
            reactionGeneMap = {}
            
            for (rxn,gene) in geneMap:
                if gene not in reactionGeneMap.keys():
                    reactionGeneMap[gene] = set()
                reactionGeneMap[gene].add(rxn)            
        
        for gene in reactionGeneMap.keys():
            rxns = reactionGeneMap[gene]
            finalValue = 0.0
            rControlValues = []
            
            for rxn in rxns:
                if rxn in library.keys():
                    value = library[rxn]
                    try:
                        value = float(value)
                    except:
                        continue
                    rControlValues.append(value)
            
            if len(rControlValues) > 0:
                minValue = min(rControlValues)
                maxValue = max(rControlValues)
            else:
                minValue = 0.0
                maxValue = 0.0 
                    
            if dir == 0.0:
                if abs(minValue) >= abs(maxValue): finalValue = minValue
                else: finalValue = maxValue
            elif dir == -1.0: finalValue = minValue
            elif dir == 1.0: finalValue = maxValue
                    
            for rxn in rxns:
                if finalValue != 0:
                    result[rxn] = finalValue
                
        return result
    
    def unifyGeneControlLibraries(self,libraries,geneMap=None):
        '''
        unification upon each library in libraries list
        '''
        result = []
        for (name,prefix,controlValues) in libraries:
            iControlValues = self.changeToGeneLibrary(controlValues, geneMap, dir=0.0)
            result.append((name,prefix,iControlValues))
            
        return result
    
    def trimWeightsDirection(self,referenceMap,targetMap,inclusive=True):
        '''
        -Insert description here!
        '''
        result = {}
        for key in targetMap.keys():
            v = targetMap[key]
            if v != 0:
                if key not in referenceMap.keys():
                    result[key] = v
                else:
                    vi = referenceMap[key]
                    try: 
                        vi = float(vi)
                    except:
                        continue
                    if cmp(v,0) != cmp(vi,0) and inclusive:
                        result[key] = v
                    elif abs(v) > abs(vi):
                        result[key] = v 
        return(result)
    
    def trimPrices(self,p,n,m):
        '''
        Trims down prices / coefficients to with in boundaries.
        Not used
        '''
        result = {}
        for key in m.keys():
            v = m[key]
            if v > 0.0 or v < 0.0:
                if key not in p.keys() and key not in n.keys():
                    result[key] = v
                elif key in p.keys():
                    iv = p[key]
                    if v > iv:
                        result[key] = v
                elif key in n.keys():
                    iv = n[key]
                    if v < iv:
                        result[key] = v
        return result
                    
        
    def reducedPrices(self,model,objectiveMap,targetMap):
        '''
        @type model: LinearModel
        @var model: The linear model for which reduced prices will be found for each variable
        @type objectiveMap: {objectiveTag:{var,coeff}}
        @var objectiveMap: key values of objective names and their coeffcient for variables
        @type targetMap: {objectiveTag:ceoff}  
        
        For each objective in the objectiveMap.
         Finds reduced prices for each variable in the model.
        
        Use this method in library update instead of old library update
        '''
        result = {}
        fluxValues = {}
        for oKey in objectiveMap.keys():
            #print "\t Objective [%s]" % (oKey)
            obj = objectiveMap[oKey]
            iModel = LinearModel()
            iModel.addModel(model)
            
            lp = self.solverClass()
            lp.setModel(iModel)
            lp.clearObjective()
            lp.setObjectiveMap(obj)
            lp.runSimplex()
            prediction = lp.getPredictionMap()
            columnStatus = lp.getColumnsStatus()
            fluxValues[oKey] = prediction
            
            #! Finding reduced prices appears to alter model
            #! Hence the model must be copied before reduced prices are found
            statusReport = {}
            for sKey in columnStatus.keys():
                status = columnStatus[sKey]
                prices = lp.getPointReducedCosts({sKey:-1.0},2)
                fluxValue = prediction[sKey]
                statusReport[(sKey,status,fluxValue)] = prices
            
            for tKey in targetMap.keys():
                tValue = targetMap[tKey]
                values = lp.getPointReducedCosts(tValue,2)
                #print "\t target [%s] library size [%s]" % (tKey, len(values))
                result[(oKey,tKey)] = values
                
            lp.clear()
            del lp
        return (fluxValues,result)
    
    def reducedPriceLibrary(self,model,limitsMap,objectiveMap,targetMap):
        '''
        @var limitsMap: set of limits to create the library of reduced prices
        @var objectiveMap: set of objectives to search for reduced prices
         
        @type fluxValuesMap: {condition:{var:value}}
        @var fluxValuesMap: map of flux values for given condition
        @type coeffecentLibrary: (tag,prefix,values:{string,float})
        @var coeffecentLibrary: library of reduced prices for given limits,objectives and targets  
        
        builds libraries of reduced prices for combination of limits, objectives and targets
        '''
        
        p = {}
        n = {}
        fluxValuesMap = {}
        coeffecentLibrary = {}
        for key in limitsMap.keys():
            #print "building library %s" % (key)
            limits = limitsMap[key]
            newModel = LinearModel()
            newModel.addModel(model)
            for name in limits.keys():
                limit = limits[name]
                newModel.addColumnLimit(name,limit)
            (fluxValues, priceMap) = self.reducedPrices(newModel,objectiveMap,targetMap)
            fluxValuesMap[key] = fluxValues
            for (ikey,values) in priceMap.items():
                coeffecentLibrary[(key,ikey[0],ikey[1])]= values
            for pKey in priceMap.keys():
                prices = priceMap[pKey]
                (p,n) = self.divideMapsDir(p,n,prices)
        
        controlMap = {}
        for (key,valueMap) in coeffecentLibrary.items():
            (k1,k2,k3) = key
            pkey = "price_" + k3
            self.fillControlMap(controlMap, pkey, valueMap)
                
        return (fluxValuesMap,coeffecentLibrary,controlMap,p,n)
    
    def stepWiseReducedPrices(self,model,objective,newObjective,targetMap, maxStep = 2000, dim = 2):
        '''
        steps through the simplex solving process for a linear model 
        recording the reduced prices observers as the process proceedes.
        '''
        result = {}
        controlLib = {}
        p = {}
        n = {}
        step = -2
        
        lp = self.solverClass()
        lp.setModel(model)
        lp.clearObjective()
        lp.setObjectiveMap(objective)
        lp.runSimplex()
        pValue = lp.getPredictionMap()
        result[(step,"flux")] = pValue
        result[(step,"rs")] = lp.getRowsStatus()
        result[(step,"cs")] = lp.getColumnsStatus()

        
        for tKey in targetMap.keys():
            tValue = targetMap[tKey]
            values = lp.getPointReducedCosts(tValue,dim)
            #print "\t step [0] target [%s] library size [%s]" % (tKey, len(values))
            result[(step,tKey)] = values
            controlLib[(step,tKey)] = values
        
        lp.clearObjective()
        lp.setObjectiveMap(newObjective)
        lp.setIteration(step)
        
        scode = 0
        step = step + 1
        while step < maxStep and scode != 200:
            if step >= 0:
                lp.setIteration(step)
                scode = lp.runSimplex()
            objValue = lp.getObjectiveValue()
            pValue = lp.getPredictionMap()
            result[(step,"flux")] = pValue
            result[(step,"rs")] = lp.getRowsStatus()
            result[(step,"cs")] = lp.getColumnsStatus()
            for tKey in targetMap.keys():
                tValue = targetMap[tKey]
                values = lp.getPointReducedCosts(tValue,dim)
                #print "\t step [%s] target [%s] library size [%s]" % (step, tKey, len(values))
                result[(step,tKey)] = values
                controlLib[(step,tKey)] = values
                (p,n) = self.divideMapsDir(p,n,values)
            step = step + 1
            
        lp.clear()
        del lp
        
        controlMap = {}
        for (key,valueMap) in controlLib.items():
            (iter,target) = key
            itarget = "step_" + target
            controlMap = self.fillControlMap(controlMap, itarget, valueMap)

        return (result,controlLib,controlMap,p,n)
    
    def findSlope(self,lp,fName,fValue,objValue,iLimit,result,factor,max,useZeros):
        '''
        finds delta Z given a set limit on a objective
        leads to finding dZ/dVi
        '''
        
        delta = 1e-4
        lp.addColumnLimit(fName,iLimit,iLimit)
        sCode = lp.runSimplex()
        fluxValues = lp.getPrediction()
        oValue = lp.getObjectiveValue()
        if abs(fValue) < delta :
            pass
            return result
        if oValue != None and sCode == 200:
            oDiff = (oValue - objValue)
            rDiff = (iLimit - fValue)
            slope = -1.0 * oDiff/rDiff
            fslope = slope * factor
            result[fName] = fslope
        return result

    def findSimpleSlopes(self,model,fluxes,objective,objValue,factor = 1.0,delta = 0.001, useZeros=False,constrict=False):
        '''
        find dz/dv for variables in linear model.
        '''        
        iModel = LinearModel()
        iModel.addModel(model)

        posResult = {}
        negResult = {}
        max = abs(objValue)*2
        step = 0.1
        
        lp = self.solverClass()
        lp.setModel(iModel)
        lp.clearObjective()
        lp.setObjectiveMap(objective)
        lp.runSimplex()
        fluxes = lp.getPredictionMap()
        
        if constrict:
            for name in fluxes.keys():
                value = fluxes[name]
                lower = None
                upper = None
                if value >= 0 :
                    upper = value
                if value <= 0:
                    lower = value
                    
                lp.addColumnLimit(name,lower,upper)
                
        for fName in fluxes.keys():
            fValue = fluxes[fName]
            (lower,upper) = model.getColumnLimit(fName)
            
            if lower == None:
                lower = float("-inf")
            if upper == None:
                upper = float("inf")
            
            if abs(fValue) <= delta:
                
                if abs(lower) != 0.0 and lower != float("-inf"):
                    iLowerLimit = -step
                else:
                    iLowerLimit = -step
                if abs(upper) != 0.0 and upper != float("inf"):
                    iUpperLimit = step
                else:
                    iUpperLimit = step
            else:
                iLowerLimit = fValue * 0.9
                iUpperLimit = fValue * 1.10
            
            negSlope = None
            if iLowerLimit > lower and iLowerLimit != fValue:
                negResult = self.findSlope(lp,fName,fValue,objValue,iLowerLimit,negResult,factor,max,useZeros)
                lp.addColumnLimit(fName,lower,upper)
    
            posSlope = None
            if iUpperLimit  < upper and iUpperLimit != fValue:
                posResult = self.findSlope(lp,fName,fValue,objValue,iUpperLimit,posResult,factor,max,useZeros)
                lp.addColumnLimit(fName,lower,upper)
                
            if lower == float("-inf"):
                lower = None                
            if upper == float("inf"):
                upper = None
            lp.addColumnLimit(fName,lower,upper)
            
        lp.clear()
        del lp
        return (negResult,posResult)
    
    def findBothSimpleSlopes(self,modelMatrix,naturalObjective,syntheticObjective,oPredVal,oBioVal,sPredVal,sSynthVal,cFactor,sFactor):
        '''
        find dz/dv at maxed objective values  
        '''
        
        (sPos,sNeg) = ({},{})
        
        (negBio,posBio) = self.findSimpleSlopes(modelMatrix,oPredVal,naturalObjective,-oBioVal, factor = cFactor)
        (negSynth,posSynth) = self.findSimpleSlopes(modelMatrix,sPredVal,syntheticObjective,-sSynthVal, factor = sFactor)
        
        (sPos,sNeg) = self.divideMapsDir(sPos,sNeg,negBio)
        (sPos,sNeg) = self.divideMapsDir(sPos,sNeg,posBio)
        
        (sPos,sNeg) = self.divideMapsDir(sPos,sNeg,negSynth)
        (sPos,sNeg) = self.divideMapsDir(sPos,sNeg,posSynth)
        
        controlMap = {}
        controlMap = self.fillControlMap(controlMap, "delta_natural", negBio)
        controlMap = self.fillControlMap(controlMap, "delta_natural", posBio)
        controlMap = self.fillControlMap(controlMap, "delta_synth", negSynth)
        controlMap = self.fillControlMap(controlMap, "delta_synth", posSynth)
        
        controlLibrary = []
        controlLibrary.append(("natural_adjust_pos","_bap_",posBio))
        controlLibrary.append(("natural_adjust_neg","_ban_",negBio))
        controlLibrary.append(("synth_adjust_pos","_sap_",posSynth))
        controlLibrary.append(("synth_adjust_neg","_san_",negSynth))
        
        return (sPos,sNeg,controlLibrary,controlMap)
    
    def generateAdjustLibrary(self,model,naturalObjective,syntheticObjective,nfactor,sfactor):
        '''
        generate library from pertibation of fluxes and observed response of objective
        '''
        lo = self.solverClass()
        naturalObjectiveName = naturalObjective.keys()[0]
        syntheticObjectiveName = syntheticObjective.keys()[0]
       
        oPredVal = lo.run(model=model,objective = naturalObjective)
        oBioVal = oPredVal[naturalObjectiveName]
        #oSynthVal = oPredVal[syntheticObjectiveName]
       
        sPredVal = lo.run(objective = syntheticObjective)
        #sBioVal = sPredVal[naturalObjectiveName]
        sSynthVal = sPredVal[syntheticObjectiveName]
        lo.clear()
              
        (pos,neg,controlLibrary,controlMap) = self.findBothSimpleSlopes(model,naturalObjective,syntheticObjective,oPredVal,oBioVal,sPredVal,sSynthVal,nfactor,sfactor)
        
        return (pos,neg,controlLibrary,controlMap)
    
    def sensitivityLibrary(self,model,targets,objectiveName,factor=1.0,prefix =""):
        '''
        find starting state of objective model
        then find sensitivity of objective to values of variables in model 
        '''
        lo = self.solverClass()
        predValues = lo.run(model=model, objective = {objectiveName:-1.0})
        objVal = predValues[objectiveName]
        
        lmSensitivity = LinearModelSensitivity()
        lmSensitivity.verbose = False
        lmSensitivity.objectiveRange = 0.8
        lmSensitivity.searchIterMax = self.searchIterMax #!Makes big difference in search times, probably change back to 30 or after testing phase
        if self.fluxBounds != None:
            lmSensitivity.fluxBounds = self.fluxBounds
        (fluxBoundaries, isLibrary, senseValues) = lmSensitivity.sensitivityAnalysisSampling(model, objVal, objectiveName, targets)
        (nsNeg,nsPos) = lmSensitivity.unifyLibraryValues(isLibrary, factor)
        
        sLibrary = []
        sLibrary.append((prefix + "sense_neg_%s" % (objectiveName),prefix + "_scn_%s_" % (objectiveName),nsNeg))
        sLibrary.append((prefix + "sense_pos_%s" % (objectiveName),prefix + "_scp_%s_" % (objectiveName),nsPos))
     
        controlSet = {}
        for (key,valueMap) in senseValues.items():
            if key not in controlSet.keys():
                controlSet[key] = set()
                for (vkey,value) in valueMap.items():
                    controlSet[key].add(value)
        
        controlSet = self.multiplySetMap(controlSet,factor)
     
        return (fluxBoundaries,sLibrary,controlSet)
   
   
    def unifyLibraries(self,controlLibraries, posControl = {}, negControl = {}):
        '''
        Join libraries and keep maximum positive and minimum negative values
        '''
        
        for (name,prefix,lib) in controlLibraries:
            (posControl,negControl) = self.divideMapsDir(posControl,negControl,lib)
    
        return(posControl,negControl)
    
    def cleanUpControl(self,controlLibraries,targets):
        icl = []
        for (n,c,d) in controlLibraries:
            id = {}
            for t in targets:
                if t in d.keys():
                    id[t] = d[t]
            icl.append((n,c,id))
        controlLibraries = icl
        
        return controlLibraries
    
    def generateCompressedLibraries(self,controlLibraries):
        if self.verbose: print "--Compressing libraries"
        (posControl,negControl) = self.unifyLibraries(controlLibraries, posControl={}, negControl={})
        controlLibraries = []
        controlLibraries.append(("postive_prime", "_pp_", posControl))
        controlLibraries.append(("negative_prime", "_np_", negControl))
        
        return controlLibraries
    
    def _filterControl(self,controlLibraries,targets=None):
        result = []
        for (n,c,d) in controlLibraries:
            id = {}
            for (k,v) in d.items():
                if self.minAbsControl <= abs(v) <= self.maxAbsControl:
                    if targets == None or k in targets:
                        id[k] = v
                    else:
                        pass
            result.append((n,c,id))
        
        return result                
    
    def _libraryUnification(self,controlLibraries):
        result = {} 
        
        tkeys = set()
        for (name,prefix,lib) in controlLibraries:
            keys = lib.keys()
            tkeys = tkeys.union(keys)
        
        for key in tkeys:
            v = 0
            t = 0
            for (name,prefix,lib) in controlLibraries:
                if key in lib.keys():
                    vi = lib[key]
                    v = v + vi
                    t = t + 1
            if v != 0:
                value = v/t
                result[key] = value
            
        return result
        
    def generateControlLibraray(self,model,limitsMap,objectivesMap,targetsMap,naturalObjective,syntheticObjective,nfactor=1.0,sfactor=1.0,targets=None,targetsOnly=False):
        '''
        currently the primary method
        generate a library of sensitivity coefficients from various methods.
        possibly compress it and return it for use with redirector.
        '''
        
        preloadFile = self.libraryDirectory + "ControlLibrary_%s_%s_%s" % (self.naturalObjectiveName,self.syntheticObjectiveName,len(targets))
        
        if self.preload != 0 or False:
            if self.preload != -1:
                preloadFile = self.libraryDirectory + "ControlLibrary_%s_%s_%s" % (self.naturalObjectiveName,self.syntheticObjectiveName,self.preload)
                if self.verbose: print "finding preload file: %s" % preloadFile
            loadedData =  self.load(preloadFile)
            if loadedData != None:
                if self.verbose: print "preloading control libraries %s" % (preloadFile)
                (controlLibraries,controlMap) = loadedData
                return (controlLibraries,controlMap)
        else:
            if self.verbose: print "Not able to preload value = [%s]" % self.preload
        
        controlLibraries = []
        controlMap = {}
        
        if self.usePrices:
            if self.verbose: "Finding reduced price controls"
            (fluxValues,controlLib,priceControlMap,reducePos,reduceNeg) = self.reducedPriceLibrary(model,limitsMap,objectivesMap,targetsMap)
            controlMap.update(priceControlMap)
            
            controlLibraries.append(("positive_price_control","_p_price_r_",reducePos))
            controlLibraries.append(("negative_price_control","_n_price_r_",reduceNeg))
            
        if self.useStepPrices:
            if self.verbose: "Finding simplex step controls"
            (libStepMap,controlLib,stepControlMap,stepPos,stepNeg) = self.stepWiseReducedPrices(model,naturalObjective,syntheticObjective,targetsMap)
            
            controlLibraries.append(("postive_step_control", "_psp_", stepPos))
            controlLibraries.append(("negative_step_control", "_nsp_", stepNeg))
            
            controlMap.update(stepControlMap)
        
        if self.useAdjustedLevels:
            if self.verbose: "Finding small perturbation controls" 
            
            (adjustPos,adjustNeg,adjustControlLibrary,adjustControlMap) = self.generateAdjustLibrary(model,naturalObjective,syntheticObjective,nfactor,sfactor)
            
            controlMap.update(adjustControlMap)
            
            if self.verbose: print "===using adjusted flux level prices=="
            controlLibraries.append(("postive_by_slope", "_ps_", adjustPos))
            controlLibraries.append(("negative_by_slope", "_ns_", adjustNeg))
            
        #cleanup control library for only controls of included targets
        if targetsOnly:
            (controlLibraries) = self.cleanUpControl(controlLibraries, targets)
        
        if self.useSensitivity:
            if self.verbose: print "Finding sensitivity based controls"
            (nfluxBoundaries, nsLibrary,senseNatControlSet) = self.sensitivityLibrary(model, targets, self.naturalObjectiveName, nfactor)
            (sfluxBoundaries, ssLibrary,senseSynthControlSet) = self.sensitivityLibrary(model, targets, self.syntheticObjectiveName, sfactor) 
            
            controlMap["sense_natural"] = senseNatControlSet
            controlMap["sense_synth"] = senseSynthControlSet
            
            controlLibraries.extend(nsLibrary)
            controlLibraries.extend(ssLibrary)
            
            #New library combining sensitivity library for control loading.
            nsControl = [("uni_con_n_sense_","_uncs_",self._libraryUnification(nsLibrary))]
            ssControl = [("uni_con_s_sense_","_uscs_",self._libraryUnification(ssLibrary))]
            
            controlLibraries.extend(nsControl)
            controlLibraries.extend(ssControl)
        
        #check to remove coefficients that are too large or too small
        controlLibraries = self._filterControl(controlLibraries)
 
        #save working control libraries for faster start.
        if preloadFile != '':    
            loadFile = open(preloadFile,'w')
            pickle.dump((controlLibraries,controlMap),loadFile)
            loadFile.close()
               
        return (controlLibraries,controlMap)
    
    
    def generateControl(self,model,targets,targetsOnly=True):
        
        naturalObjective = {self.naturalObjectiveName:-1.0}
        syntheticObjective = {self.syntheticObjectiveName:-1.0}
        
        lo = self.solverClass()
        oPredVal = lo.run(model=model)
        oBioVal = oPredVal[self.naturalObjectiveName]
        oSynthVal = oPredVal[self.syntheticObjectiveName]
        
        if self.verbose: print "Max Natural: Biological [%s], Synthetic [%s]" % (oBioVal,oSynthVal)
        
        sPredVal = lo.run(model=model, objective = syntheticObjective)
        sBioVal = sPredVal[self.naturalObjectiveName]
        sSynthVal = sPredVal[self.syntheticObjectiveName]
        
        lo.clear()
            
        if self.verbose: print "Max Synthetic: Biological [%s], Synthetic [%s]" % (sBioVal,sSynthVal)
            
        objectivesMap = {}
        targetsMap = {}
        limitsMap = {}
        
        #Normalize controls to magnitude of objectives
        
        nFactor = -1.0
        sFactor = 1.0*oBioVal/sSynthVal
        
        if self.verbose:  print "==========Building  control libraries=========="
        
        #Construct reduce prices targets and objectives
        
        if self.usePrices or self.useStepPrices:
            print "--using simplex opportunity costs--"
            objectivesMap["cellular"] = naturalObjective
            objectivesMap["synthetic"] = syntheticObjective
            
            targetsMap["cellular"] = {self.naturalObjectiveName:nFactor}
            targetsMap["synthetic"] = {self.syntheticObjectiveName:sFactor}
            
            limitsMap["None"] = {self.naturalObjectiveName:(None,None)}
            #limitsMap["minBio"] = {objective:(minObjVal,None)}
            #limitsMap["negBio"] = {objective:(-1.0*oBioVal,None)}
            #limitsMap["negSynth"] = {syntheticObjectiveName:(-1.0*sSynthVal,None)}
        
        #!Note change to only return libraries for targets
        (controlLibraries,controlMap) = self.generateControlLibraray(model,limitsMap,objectivesMap,targetsMap,naturalObjective,syntheticObjective,nFactor,sFactor,targets)
                    
        if self.compressLibraries:
            controlLibraries = self.generateCompressedLibraries(controlLibraries)
    
        if self.adjustGeneControl:
            controlLibraries = self.unifyGeneControlLibraries(controlLibraries, self.geneMap)
                 
        return (controlLibraries,controlMap)
    
    def geneControlReport(self,model,controlMap,outputName):
        report = Report()
        geneMap = model.controlMap
        for (target,controls) in geneMap.keys():
            for control in controls:
                report.add(target,control,"control")
        
        writer = ReportWriter()
        writer.setFile(outputName)
        
        writer.write(report)
        writer.closeFile()
    
        return None
    
    def controlReport(self,model,controlLibraries,controlMap,targets,outputName):
        report = Report()
        controlMapList = {}
        if self.geneMap:
            report.addColumnHash("geneticControl",model.controlMap)
        
        for (name, prefix, coeffMap) in controlLibraries:
            report.addColumnHash(name,coeffMap)
            for (target,value) in coeffMap.items():
                if target not in controlMapList.keys():
                    controlMapList[target] = {}
                if value not in controlMapList[target].keys():
                    controlMapList[target][value] = []
                controlMapList[target][value].append(name) 
        
        rControlMap = {}
        tags = controlMap.keys()
        tags.sort()
        for tag in tags:
            valueSetMap = controlMap[tag]
            report.addColumnHash(tag,valueSetMap)
            for (k,valueSet) in valueSetMap.items():
                if k not in rControlMap.keys():
                    rControlMap[k] = {}
                rValueMap = rControlMap[k]
                for v in valueSet:
                    if v not in rValueMap.keys():
                        rValueMap[v] = []
                    rValueMap[v].append(tag)
                rControlMap[k] = rValueMap

        iControlMapList = {}
        for (key,data) in controlMapList.items():
            if data != None:
                list = data.keys()
                l = reduceVector(list)
                iControlMapList[key] = l

        targetNameMap = {}
        for target in targets:
            targetNameMap[target] = "Target"
    
        report.addColumnHash("control target", targetNameMap)        
        report.addColumnHash("complete list of control Library", controlMapList)    
        report.addColumnHash("control Reduced", iControlMapList)
        report.addColumnHash("complete control discovered", rControlMap)
        writer = ReportWriter()
        writer.setFile(outputName)
        
        writer.write(report)
        writer.closeFile()
        
        return None
        
        
    def updateLibrary(self,model,naturalObjectiveName,syntheticObjectiveName,coeffLibrary,factor1,factor2,mfactor=1.0):
        '''
        Used to adjust the values of a set of libraries
        Needs to be worked on so that
        A) model is not altered during update (largely from reduced price discovery)
        B) more that basic reduced prices are discovered if its deemed worth while 
        
        '''
        matrixTools = MatrixTools()
        iModel = LinearModel()
        iModel.addModel(model)
        
        lo = self.solverClass()
        lo.setModel(iModel)
        lo.runSimplex()
        tpv = lo.getPredictionMap()
        
        iPosCoeff = {}
        iNegCoeff = {}

        iPosCoeff = lo.getPointReducedCosts({naturalObjectiveName: factor1},2)
        iNegCoeff = lo.getPointReducedCosts({syntheticObjectiveName: factor2},2)
                
        lo.clear()
            
        iPosCoeff = matrixTools.vectorMScalar(iPosCoeff,mfactor)
        iNegCoeff = matrixTools.vectorMScalar(iNegCoeff,mfactor)
        
        for (name,prefix,value) in coeffLibrary:
            iPosCoeff = self.trimWeightsDirection(value,iPosCoeff,inclusive=True)
            iNegCoeff = self.trimWeightsDirection(value,iNegCoeff,inclusive=True)
                
        if self.geneMap != None:
            iPosCoeff = self.changeToGeneLibrary(iPosCoeff,self.geneMap,dir=1.0)
            iNegCoeff = self.changeToGeneLibrary(iNegCoeff,self.geneMap,dir=-1.0)
            
        return(iPosCoeff,iNegCoeff)
    
    def reduceControl(self,model,objective,target,minTarget,targets=None):
        '''
        filters control library to values absolute value  > min target
        '''
        iobjective = {}
        iobjective.update(objective)
        lp = self.solverClass()
        if targets == None:
            targets = iobjective.keys()
        
        ffluxes = {}
        for t in targets:
            if t not in iobjective.keys():
                continue
            value = iobjective[t]
            del iobjective[t]
            ifluxes = lp.run(model=model, objective = iobjective)
            if abs(ifluxes[target]) < abs(minTarget):
                iobjective[t] = value
            else:
                ffluxes=ifluxes
        
        return (ffluxes,iobjective)
    
    def controlMapReport(self,controlMap):
        '''
        convert control map to report object
        '''
        
        report = Report()
        
        for (key,coeffMap) in controlMap.items():
            colName = str(key)
            report.addColumnHash(colName, coeffMap)
            
        return report
    
    def _reduceVector(self, data, delta = 0.1,min=1e-4):
        '''
        reduce vector to values differing by at least delta percentage
        '''
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
        
    def reduceControlLibrary(self, controlLibrary):
        '''
        Reduce control libraries so that each one contains only values differing by at least some value
        '''
        
        controlMapList = {}
        
        for (name, pefix, coeffMap) in controlLibrary:
            for target in coeffMap.keys():
                if target not in controlMapList.keys():
                    controlMapList[target] = []
                var = coeffMap[target]
                controlMapList[target].append(var)
        
        iControlMapList = {}
        for (key,list) in controlMapList.items():
            if list != None:
                l = self.reduceVector(list)
                iControlMapList[key] = l
                
        return iControlMapList
    
    def flattenControlLibrary(self,controlLibrary,controlPrefixMap):
        result = []
        for (libName,libData) in controlLibrary.items():
            prefix = controlPrefixMap[libName]
            ilibData = {}
            for (t,v) in libData.items():
                if v < 0:
                    ilibData[t] = -1.0
                if v > 0:
                    ilibData[t] = 1.0
            result.append((libName,prefix,ilibData))
        return result
    
    def seedTestLibrary(self,targets,randomize = False, factor = 0.0):
        posLib = {}
        negLib = {}
        rand = Random()
        rand.seed(time())
        
        for t in targets:
            r = 1.0
            if randomize:
                r = rand.random()
            posLib[t] = 1.0 + factor * r
            negLib[t] = -1.0 - factor * r
            
        resultLib = []
        resultLib.append(("flat_pos","_fp_test_",posLib))
        resultLib.append(("flat_neg","_fn_test_",negLib))
        
        return resultLib

    def generateBinaryControl(self,targets,binaryRange):
        resultLib= []
        for n in binaryRange:
            value = pow(2,n)
            controlPos={}
            controlNeg={}
            for t in targets:
                controlPos[t] = value
                controlNeg[t] = -value
                
            nTag = ("%s" % (n)).replace("-","_")
            
            resultLib.append(("binary_pos_b%s" % (nTag),"_p_b%s_" % (nTag), controlPos))
            resultLib.append(("binary_neg_b%s" % (nTag),"_n_b%s_" % (nTag), controlNeg))
        return resultLib