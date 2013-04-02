'''
Updated February 07, 2013
@author: Graham Rockwell
Constructor for Redirector Bilevel Optimization Model
'''

from core.model.LinearModel import LinearModel
from core.model.LPSolver import LPSolver
from core.model.ProcessModel import MatrixTools
from core.model.OptimizationAnalysisTools import OptimizationDebugingTools
from core.model.ConstructRegulationOptimization import ConstructRegulationOptimization
from core.util.Report import Report
from core.util.ReportWriter import ReportWriter
from core.reader.FlatFileParser import FlatFileParser

from time import time,strftime
import pickle,os,re
from numpy import mean


class OptimizationControlRedirector:

    def __init__(self):
        #boolean
        self.verbose = False
        self.useUpdateLibrarys = True
        self.useScip = False
        self.debug = False
        self.useGeneMap = False
        self.iterReport = False
        self.orderGenes = False
        
        #integers
        self.preStart = 0
        self.iterations = 0
        self.preStartSearchSize = 0
        
        #float
        self.naturalScale = 0.2
        self.aneal = 0.0
        self.delta = 0.1
        
        #strings
        self.control = ''
        self.naturalObjectiveName = ''
        self.syntheticObjectiveName = ''
        self.modelName = ''
        self.preStartControl = ''
        
        #Data
        self.targets = []
        self.rGeneMap = None
        self.controlLibraries = None
        self.previousState = {}
        
        #objects / factories
        self.solver = 'GLPK'
        self.modelFactory = None
        self.con = ConstructRegulationOptimization()
        self.controlFactory = None 
        self.debugTools = OptimizationDebugingTools()
        self.solverConfig = ''
        
        #File IO variables:
        self.resultDirectory = ''
        self.analysisDirectory = 'ControlSearches'
        self.resultTag = 'norm'
    
    def _solverFactory(self):
        '''
        This method should be updated to allow for every PuLP solver to be accessed.
        '''
        if self.solver == "PuLP_GLPK":
            solver= LPSolver()
            solver.Mip = False
        elif self.solver == "PuLP_SCIP":
            solver = LPSolver()
            solver.Mip = True
            if self.solverConfig != '': 
                solver.setConfigFile(self.solverConfig)
        else:
            solver= LPSolver()
            print "Using default solver"
        
        return solver
    
    def reconstructControlModel(self,originalModel,controlNames,ltPredV):
        '''
        Rebuilds optimization problem using constructor
        Sets state of optimization to current search location
        '''
        if self.verbose: print "==constructing current model=="
        tmodel = self.con.construct(originalModel,controlNames,self.syntheticObjectiveName)
                
        if len(ltPredV) > 0:
            tmodel = self.con.iterate(tmodel, controlNames, self.con.libdata, ltPredV)
        
        return tmodel

    def runOptimization(self,model,tag=''):
        '''
        Method to select and run solver
        '''
        if self.verbose: print "Using solver [%s]" % (self.solver)
        
        solver = self._solverFactory()
                        
        if self.verbose: print "Preparing optimization"
        solver.setModel(model)
        
        if self.debug:
            solver.writeMPS(self.analysisDirectory + "/RD_log_%s.mps" % tag)
            solver.writeLP(self.analysisDirectory + "/RD_log_%s.lp" % tag)
        
        if self.verbose: print "Start Optimization"             
        sCode = solver.runIntOpt()
            
        return (solver,sCode)
    
    def _getGeneNames(self,reactionName):
        '''
        Returns control tags from reaction names
        '''
        
        if reactionName not in self.rGeneMap.keys():
            return None
        geneNames = self.rGeneMap[reactionName]
        if type(geneNames) == type(""):
            return set([geneNames])
        else:
            return geneNames
        
    def _getGeneTag(self,reactionName,geneClusters = None):
        '''
        Finds control tags using annotation (usually B numbers to Gene IDs)
        '''
        
        geneNames = self._getGeneNames(reactionName)
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
        geneTag = "[" + geneTag[1:] + "]"
        
        return geneTag
            
    
    def printGeneObjective(self,iObjective,geneClusters=None):
        '''
        Returns objective in terms of control tags, usually Gene IDs
        '''
        iGeneObjective = {}
        iOtherObjective = {}
        igControl = {}
        for rxnName in iObjective.keys():
            rControl = iObjective[rxnName]
            if rxnName not in self.rGeneMap.keys():
                iOtherObjective[rxnName] = rControl
            else:
                geneTag = self._getGeneTag(rxnName, geneClusters)
                
                if geneTag not in iGeneObjective:
                    igControl[geneTag] = []
                igControl[geneTag].append(rControl)
        
        for k in igControl.keys():
            iGeneObjective[k] = mean(igControl[k])
            
        return (iGeneObjective,iOtherObjective)
    
    def writeLog(self,startTime,i,iObjective,iGeneObjective,effectValue,searchNeighborhood,ltSynthVal,checkSynthVal):
        logdir = self.analysisDirectory
        runTime = time() - startTime
        naturalControl = iObjective[self.naturalObjectiveName]
        tagName = "rd_%s_%s_%s" % (self.control, self.modelName, self.syntheticObjectiveName)
        logFileName = logdir + "/" + tagName + "_log.txt"
        logFileH = open(logFileName,'a')
        timeStamp = strftime('%Y%m%dT%H%M%S')
        logLine = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (tagName,timeStamp,self.syntheticObjectiveName,"K%s"% searchNeighborhood,"%s"%i,ltSynthVal,checkSynthVal,runTime,naturalControl,effectValue,iObjective)
        if self.useGeneMap:
            mapLine = '\t%s\t%s' % (len(iGeneObjective),iGeneObjective)
            logLine = logLine + mapLine
        logLine = logLine + '\n'
        logFileH.write(logLine)
        logFileH.close()
        
        return None
    
    def getControlNames(self,targets):
        tControls = []
        uControls = set()
        
        for pair in self.con.libdata:
            libName = pair[0]
            libPrefix = pair[1]
            iControls = self.con._getNonZeroRegulationNames(libName,targets)
            tControls.extend(iControls)
            uControls = uControls.union(set(iControls))
            #if self.verbose: print "Library %s [%s] size %s" % (libName,libPrefix,len(iControls))
                     
        controlNames = set(uControls)
        controlNames = controlNames.intersection(self.targets)
        if self.verbose: print "Unique Control Variables size [%s]" % (len(uControls))
        if self.verbose: print "Total Control Variables size [%s]" % (len(tControls))
        return controlNames
        
    def compareTargets(self,model,targets1,targets2):
    
        tmodel1 = self.reconstructControlModel(model,targets1,{})
        tmodel2 = self.reconstructControlModel(model,targets2,{})
        
        keys1 = set(tmodel1.keys())
        keys2 = set(tmodel2.keys())
        
        for k1m in keys1.difference(keys2):
            print "model 1 missing %s" % (k1m)
        for k2m in keys2.difference(keys1):
            print "model 2 missing %s" % (k2m)
    
        return None
    
    def _refactorConstructor(self,con,sep="_"):
        '''
        Needs to be updated to not overwrite input constructor but create new instance.
        '''
        
        rGeneMap = con.rGeneMap
        geneList = con.geneList
        if rGeneMap != None:
            rGeneMap = self.modelFactory.refactorReactionGeneSetMap(rGeneMap, sep)
        if geneList != None:
            geneList = self.modelFactory.refactorGeneList(geneList, sep)
        
        con.rGeneMap = rGeneMap
        con.geneList = geneList
        
        return con
    
    def _refactorControlTargets(self,targets,sep):
        mfactory = self.modelFactory
        
        result = []
        for (lib,value) in targets:
            ivalue = mfactory.orderTags(value, sep)
            itarget = (lib,ivalue)
            result.append(itarget)
        return result
    
    def readAnnotatedTargets(self,targetFile,model,annotationName="bnumber",regex="[a-zA-Z0-9\(\)]+"):
        parser = FlatFileParser()
        report = parser.parseToReport(targetFile, "Row Names", header=["Row Names","Control","Flux","Reaction Control","Flux Changes"])
        controlValues = report.getColumn("Control")
        controlValues = self._deannotateGenes(controlValues, model, annotationName, regex)
        controlScore = report.getColumn("Flux")
        controlScore =  self._deannotateGenes(controlScore, model, annotationName, regex)
        reactionControl = report.getColumn("Reaction Control")
        reactionControl =  self._deannotateGenes(reactionControl, model, annotationName, regex)
        
        return(controlValues,controlScore,reactionControl)
    
    def sliceValueControl(self,value,controlMap):
        '''
        to be refactored
        '''
        result = []
        controlPile = []
        cvalue = value
        
        if value == 0:
            return result
        
        if value > 0:
            dir = 1.0
        
        if value < 0:
            dir = -1.0
        
        
        while cvalue != 0:
            pile = []
            for key in controlMap.keys():
                if cvalue*dir >= key * dir:
                    if key*dir not in controlPile:
                        pile.append(key*dir)
            if len(pile) == 0:
                return result
            ivalue = max(pile)*dir
            clib = controlMap[ivalue]
            result.append(clib)
            controlPile.append(ivalue)
            cvalue = cvalue - ivalue
        
        return result
                    
    def loadControlTable(self,model,controlMap,preStart,directory,controlAdjust=1.0,keyTag="Row Names",coeffTag="Control",productionTag="Flux"):
        '''
        Load flat file of control designs and production values
        '''
        mfactory = self.modelFactory
        (controlTag,searchSize) = self.getPreStartInfo()
        result = {}
        targets = []
        naturalCoefficient = 0
        
        fileName = "Target_Report_T_%s_C_%s_M_%s_N_%s_S_%s_K_%s_I_%s.txt" % (self.resultTag,controlTag,self.modelName,self.naturalObjectiveName,self.syntheticObjectiveName,searchSize,preStart)
        fileName = directory + fileName
        if not os.path.exists(fileName):
            print "file name not found %s" % (fileName)
            return ([],[],0.0)
        
        (controlValues,controlScore,reactionControl) = self.readAnnotatedTargets(fileName,model)
        tMap = model.getGeneTargetMap()
        
        for (geneTag,controlValueString) in controlValues.items():
            geneTag = re.sub("[\]\[]","",geneTag)
            if self.orderGenes:
                orderedTag = mfactory.orderTags(geneTag, sep="_",replace=" ")
            else:
                orderedTag = geneTag
            controlValue = float(controlValueString)
            if geneTag == self.naturalObjectiveName:
                naturalCoefficient = controlValue 
            else:
                icontrolValue = controlValue * controlAdjust
                if controlValue in controlMap.keys():
                    controlLibs = [controlMap[icontrolValue]]
                else:
                    controlLibs = self.sliceValueControl(icontrolValue, controlMap)
                if len(controlLibs) != 0:
                    for controlLib in controlLibs:
                        if controlLib not in result.keys():
                            result[controlLib] = set()
                        if orderedTag not in tMap.keys():
                            continue
                        result[controlLib].add(orderedTag)
                        rTargets = tMap[orderedTag]
                        targets.extend(rTargets)
                else:
                    if self.verbose: print "control value not found [%s]" % (controlValue)
               
        return (result,targets,naturalCoefficient)
    
    def getPreStartInfo(self):
        if self.preStartControl == '':
            controlTag = self.control
        else:
            controlTag = self.preStartControl
        if self.preStartSearchSize == 0:
            searchSize = self.con.controlMax
        else:
            searchSize = self.preStartSearchSize
        return (controlTag,searchSize)
    
    def _setControlList(self,controlTargets):
        activeTargets = {}
        
        for (lib,value) in controlTargets:
            if lib not in activeTargets.keys():
                activeTargets[lib] = set()
            activeTargets[lib].add(value)
        
        iats = []
        for (lib,values) in activeTargets.items():
            iat = (lib,values)
            iats.append(iat)
        activeTargets = iats
    
        return activeTargets
    
    def getPreStart(self,preStart,targets,sep="_"):
        (controlTag,searchSize) = self.getPreStartInfo()   
        preStartName = self.analysisDirectory + "/preStart_%s_%s_%s_%s_%s_%s" % (controlTag,self.modelName, self.naturalObjectiveName, self.syntheticObjectiveName,searchSize,preStart)
        controlTargets = []
        
        if preStart != 0 and os.path.exists(preStartName):
            preStartFile = open(preStartName,"r")
            (constructor,objective,prediction) = pickle.load(preStartFile)
            controlTargets = constructor.getControlNames(prediction.keys(),prediction,prefix = False) 
            preStartFile.close()
        else:
            return None
        
        if self.orderGenes:
            icontrolTargets = self._refactorControlTargets(controlTargets, sep)    
            controlTargets = icontrolTargets
        
        activeTargets = {}
        for (lib,value) in controlTargets:
            if lib not in activeTargets.keys():
                activeTargets[lib] = set()
            activeTargets[lib].add(value)     
                
        return (constructor,objective,prediction,activeTargets)

    def writePreStart(self,objective,prediction,iteration):
        (controlTag,searchSize) = self.getPreStartInfo()
        searchSize = self.con.controlMax
        preStartName = self.analysisDirectory + "/preStart_%s_%s_%s_%s_%s_%s" % (controlTag,self.modelName, self.naturalObjectiveName, self.syntheticObjectiveName, searchSize,iteration)
        preStartFile = open(preStartName,"w")
        pickle.dump((self.con,objective,prediction),preStartFile)
        preStartFile.close()
        return preStartName    
    
    def loadActiveControls(self,model,targets,limitTargets = False,growthParameterAdjustment = -1.0,controlMapAdjust = 1.0):
        preStart = 0
        naturalCoefficient = self.naturalScale
        activeControls = []
        
        if self.preStart != 0:
            loadData = None
            
            if self.control == "flat":
                controlMap = {1.0:"flat_neg",-1.0:"flat_pos"}
            
            if self.control == "binary":
                controlMap = {1.0:"binary_neg_b0",0.5:"binary_neg_b-1",0.25:"binary_neg_b-2",-1.0:"binary_pos_b0",-0.5:"binary_pos_b-1",0.25:"binary_pos_b-2"}
                
            (activeControls,targetsT,naturalCoefficientT) = self.loadControlTable(model,controlMap,self.preStart,directory="../../results/")
            
            if limitTargets:
                targets = targetsT #! target reduction causing problems.
            
            naturalCoefficient = naturalCoefficientT*growthParameterAdjustment
            
            if activeControls != []:
                preStart = self.preStart
            else:
                loadData = self.getPreStart(self.preStart,targets)
            
            if loadData != None:
                if self.verbose: print "Loading preStart condition"
                preStart = self.preStart
                (icon,iObjective,lPrediction,activeControls) = loadData
                naturalCoefficientT = -iObjective[self.naturalObjectiveName]
        else:
            preStart = 0
        
        return (activeControls,targets,naturalCoefficient,preStart)
        
        
    def _intitalizeModelConstructor(self,model,minNaturalValue,primeBounds = None):
        if self.verbose: print "==Formulating Framework=="    
        self.con.delta = 1e-4
        self.con.minObjective = minNaturalValue
        self.con.newObjective = self.syntheticObjectiveName
        self.con.addRegulationLibraries(self.controlLibraries)
        if primeBounds != None:
            self.con.primeControlLimits = primeBounds
        if self.useGeneMap: 
            self.con.setGeneMap(model.controlMap)
            
    def _runOptimization(self,originalModel,model,tmodel,controlNames,runTag,geneClusters=None):
        '''
        Perform optimization Redirector formulation of MILP
        Find targeted genes or reactions
        Print reactions and or genes which are found
        '''
        
        (lt,sCode) = self.runOptimization(tmodel,runTag)
        
        if self.verbose: print "Optimization Completed [%s]" % (sCode)
                
        if sCode != "OPTIMAL" and self.useScip:
            if self.verbose: print "Failure to find optimal solution: breaking search loop"
            return None
        else:
            if self.verbose: print "Optimal solution found"
        
        ltPredV = lt.getMipPredictionMap()
        lt.clear()
                
        (iObjective,iLibObjective) = self.con.generateObjective(originalModel, controlNames, ltPredV)
        
        if self.verbose: print "Current objective[%s]: %s" % (str(len(iObjective)),iObjective)
        if self.verbose: print "Current library source objective[%s]: %s" % (str(len(iLibObjective)),iLibObjective)
        (iGeneObjective,iOtherObjective) = (None,None)
        if self.useGeneMap:
            (iGeneObjective,iOtherObjective) = self.printGeneObjective(iObjective,geneClusters)
            aGeneObjective = self._annotateObjective(iObjective, model, "bnumber", regex="[a-zA-Z0-9\(\)]+")    
            if self.verbose: print "\nCurrent Genetic objective[%s]: %s" % (str(len(aGeneObjective)),aGeneObjective)
            if self.verbose: print "Current Other objective: %s" % (iOtherObjective)
            
        return (ltPredV,iObjective,iGeneObjective)
        
            
    def _getProgressiveParameter(self,oPredVal,sPredVal,iPredVal,iObjective):
        '''
        Calculate progressive growth parameter from current flux state
        '''
        
        zObjective = dict(iObjective)
        naturalCoefficient = zObjective[self.naturalObjectiveName]
        zObjective[self.naturalObjectiveName] = 0.0
        
        oNatVal = oPredVal[self.naturalObjectiveName]
        sNatVal = sPredVal[self.naturalObjectiveName] 
                
        matrixTools = MatrixTools()
        sv = False
        sumValue1 = matrixTools.vectorMultiply(oPredVal,zObjective,verbose = sv)
        sumValue2 = matrixTools.vectorMultiply(iPredVal,zObjective,verbose= sv)
        sumValue3 = matrixTools.vectorMultiply(sPredVal,zObjective,verbose= sv)
                
        rValue = sumValue1 - sumValue2 #redirector function value beyond current vs. natural optimal value
        sEffectValue = sumValue1 - sumValue3 #redirector function value synthetic vs. natural
        bioAdjust = rValue / (oNatVal-sNatVal) # new natural control value
        
        if self.verbose: print "Current Redirector Function [%s] = Z:Vo[%s] - Z:Vc[%s]" % (rValue, sumValue1, sumValue2)
        if self.verbose: print "Control Values: Natural [%s] vs. Current [%s] vs. Synthetic [%s]" % ((oNatVal-sNatVal)*naturalCoefficient, rValue, sEffectValue)
        if self.verbose: print "Adjusted control value [%s] = [%s] / ([%s] - [%s])" % (bioAdjust,rValue,oNatVal,sNatVal)
                                 
        return (rValue,sEffectValue,bioAdjust)
    
    def findCurrentObjective(self,prediction,model,originalModel,controlNames,geneClusters):
        ltPredV = prediction
        ltBioVal = ltPredV[self.naturalObjectiveName]
        ltSynthVal = ltPredV[self.syntheticObjectiveName]
        resultPass = True

        if self.verbose: print "\n======================Optimization Result======================"
        (iObjective,iLibObjective) = self.con.generateObjective(originalModel, controlNames, ltPredV)
        if self.verbose: print "Current objective[%s]: %s" % (str(len(iObjective)),iObjective)
        if self.verbose: print "--> library source objective[%s]: %s" % (str(len(iLibObjective)),iLibObjective)
        (iGeneObjective,iOtherObjective) = (None,None)
        if self.useGeneMap:
            (iGeneObjective,iOtherObjective) = self.printGeneObjective(iObjective,geneClusters)
            aGeneObjective = self._annotateObjective(iObjective, model, "bnumber", regex="[a-zA-Z0-9\(\)]+")    
            if self.verbose: print "Current Genetic objective[%s]: %s" % (str(len(aGeneObjective)),aGeneObjective)
            if self.verbose: print "Current Other objective: %s" % (iOtherObjective)
            
        #--------------------------------------
        #Check prediction
        #--------------------------------------
        if(True):
            lpCheck = self._solverFactory()
            checkPredV = lpCheck.run(model=originalModel,objective = iObjective)
            checkNaturalValue = checkPredV[self.naturalObjectiveName]
            checkSynthVal = checkPredV[self.syntheticObjectiveName]
            lpCheck.clear()
            if self.verbose: print "---Optimization Result check ---"
            if self.verbose: print "Natural [%s] ?= Check [%s] " % (ltBioVal, checkNaturalValue)
            if self.verbose: print "Production [%s] ?= Check [%s]" % (ltSynthVal, checkSynthVal)
            
        return (iObjective,iGeneObjective,ltSynthVal,checkPredV,checkSynthVal) 
        
    def progressiveSearchParameter(self,pred,pPred,naturalOpt,synthOpt,iObjective,naturalCoefficient,pNaturalCoefficient,pEffectValue,maxNatural,maxEffect):
        
        #--------------------------------------------------------
        # Calculate updated natural parameter 
        #--------------------------------------------------------
                
        if self.verbose: print "\n======================Progress Control======================"
        (effectValue,sEffectValue,bioAdjust) = self._getProgressiveParameter(naturalOpt, synthOpt, pred, iObjective)
    
        ltSynthVal = pred[self.syntheticObjectiveName]
        ltBioVal = pred[self.naturalObjectiveName]
        pBioV = pPred[self.naturalObjectiveName]
        pSynthV = pPred[self.syntheticObjectiveName]
        
        delta = self.delta
        
        oBioVal = naturalOpt[self.naturalObjectiveName]
        bioValLower = oBioVal * 0.5 
        synthValLower = synthOpt[self.syntheticObjectiveName] *0.1
        
        #-----------------------------------------
        #Search Loop Conditions
        #-----------------------------------------
        
        if self.verbose: print "===Progressive Growth Search Conditions==="
        if self.verbose: print "Effect %s vs. Current max effect %s" % (effectValue, maxEffect)
        if self.verbose: print "Redirector function control value %s" % (bioAdjust)
        if self.verbose: print "Natural objective value %s vs. Min natural value %s " % (ltBioVal, bioValLower)
        if self.verbose: print "Production value %s vs. Min production %s " % (ltSynthVal, synthValLower)
                
        #-------------------------------------------------------------------
        # Find new progressive control value
        #-------------------------------------------------------------------                
                                
        print "\n=============Finding new progressive control value=============="
        lowSynth = abs(ltSynthVal) < abs(synthValLower)
        lowBio = abs(ltBioVal) < abs(bioValLower)
        staticEffect = abs(effectValue - pEffectValue) < pEffectValue * 2 * delta
        staticNatural = abs(naturalCoefficient - pNaturalCoefficient) < pNaturalCoefficient * 2 * delta
                
        if  staticEffect or lowSynth or lowBio:
            print "-----Finding new progressive control value----"
            print "Adjustment: previous [%s] -> current [%s]" %(pEffectValue, effectValue)
            print "Previous Bio[%s] ->  CurrentBio[%s]  vs.  MinBio[%s] " %(pBioV, ltBioVal,bioValLower)
            print "Previous Synth[%s] -> CurrentSynth[%s]  vs.  MinSynth[%s]" %(pSynthV, ltSynthVal,synthValLower)
            print "Natural control value [%s]" % (naturalCoefficient) 
                    
            if staticEffect:
                print "-----Unchanged Control Contribution-----"
                print "Natural Parameter Change [%s]" % abs(naturalCoefficient - pNaturalCoefficient)
                print "Redirection Function Change [%s]" % abs(effectValue - pEffectValue)  
                        
            if lowBio:
                print "-----Accounting for low natural objective flux-----"               
                pNaturalCoefficient = naturalCoefficient
                naturalCoefficient = bioAdjust * (1.0 + delta)
                if naturalCoefficient > pNaturalCoefficient *1.20: 
                    print "New Natural Control [%s] =  %s / ( %s - %s ):" % (bioAdjust, effectValue, oBioVal, ltBioVal)
                else:
                    naturalCoefficient = pNaturalCoefficient * 1.20
                    print "New Natural Control [%s] = 1.2 * Old Control [%s]" % (naturalCoefficient,pNaturalCoefficient)
                    
            elif lowSynth:
                print "-----Low Production Decreasing natural coefficient-----"
                print "Current Synthetic [%s] < Min Synthetic [%s]" % (abs(ltSynthVal),abs(synthValLower))
                print "Reducing natural control strength"
                        
                if not staticNatural:
                    pNaturalCoefficient = naturalCoefficient
                    naturalCoefficient = pNaturalCoefficient - ( pNaturalCoefficient * delta )
                    print "Small natural coefficient reduction"
                    print "New Natural %s = %s - (%s * %s)" %(naturalCoefficient,pNaturalCoefficient,pNaturalCoefficient,delta)
                else: 
                    pNaturalCoefficient = naturalCoefficient
                    naturalCoefficient = maxNatural + (pNaturalCoefficient - maxNatural)/2
                    print "Mid point natural coefficient reduction"
                    print "New Natural %s = %s + (%s - %s) / 2" %(naturalCoefficient,maxNatural,pNaturalCoefficient,maxNatural)
                                    
            return (effectValue, naturalCoefficient)
    
    def optimizeControl(self,model,minNatural=0.2):
        '''
        Main Function 
        Run: Growth Drive Progressive Target Discovery
        '''
        
        '''
        1. Derive optimal natural and synthetic networks 
           Find maximum Natural Objective
           Find minimum threshold for to use for Natural objective
        '''
       
        # Set variables !working point
        self.rGeneMap = model.controlMap
        self.modelName = model.modelName
        self.targets = model.targets
        geneClusters = model.controlClusters
        
        if self.verbose: print "Finding natural and synthetic parameters of metabolic network."
        
        naturalObjective = {self.naturalObjectiveName:-1.0}
        syntheticObjective = {self.syntheticObjectiveName:-1.0}
        
        lo = self._solverFactory()
        oPredVal = lo.run(model=model,objective = naturalObjective)
        oBioVal = oPredVal[self.naturalObjectiveName]
        oSynthVal = oPredVal[self.syntheticObjectiveName]
            
        minNaturalVal = oBioVal * minNatural
       
        model.addColumnLimit(self.naturalObjectiveName,(minNaturalVal,None))
        if self.verbose: print "Max Natural: Natural [%s] Synthetic: [%s] => min Natural [%s]" % (oBioVal, oSynthVal, minNaturalVal)
        #natObjLimit = model.getColumnLimit(self.naturalObjectiveName)
        
        '''
        2. Find maximum value of Synthetic objective
        '''
        
        lo = self._solverFactory()
        sPredVal = lo.run(model=model,objective = syntheticObjective)
        sBioVal = sPredVal[self.naturalObjectiveName]
        sSynthVal = sPredVal[self.syntheticObjectiveName]
        lo.clear()
        
        if self.verbose: print "Max Synthetic: Natural [%s] Synthetic: [%s]" % (sBioVal, sSynthVal)

        #--------------------------
        # Set values
        #--------------------------
        
        self._intitalizeModelConstructor(model, minNaturalVal, primeBounds=None)
        
        searchNeighborhood = self.con.controlMax
        targets = self.targets
        originalFluxNames = model.getColumnNames()
        
        #-----------------------------------
        # Original Model Objects
        #-----------------------------------
        
        originalModel = LinearModel()
        originalModel.extend(model)
        newModel = LinearModel()
        newModel.extend(model)

        #-----------------------------------------------------------------------------------
        # Progressive target discovery search parameters
        # These controls may be important to tune they should possibly be added to config
        #-----------------------------------------------------------------------------------
        
        pNaturalCoefficient = 0.0
        bioValLower = oBioVal * 0.5 
        synthValLower = sSynthVal *0.1
        delta = self.delta
        effectValue = 0.0
        pBioV = oBioVal
        pSynthV = 0.0
        pEffectValue = 0.0
        
        #-------------------------------------------
        # Progressive Growth Discovery Variables
        #-------------------------------------------
        
        if self.previousState == {}:
            ltPredV = oPredVal
        else:
            ltPredV = self.previousState
        
        maxEffect = 0
        maxNatural = 0
        finalPredictionValue = None
        finalCheckValue = None    
        iObjective = {self.syntheticObjectiveName:-1}
        checkPredV = {}
    
        #----------------------------------------------
        # Check constructor for proper gene ordering
        #----------------------------------------------
        if self.orderGenes:
            self.con = self._refactorConstructor(self.con, sep="_")
        
        #----------------------------------------
        # Find active gene sets before search
        # Use selected starting point
        #----------------------------------------
        
        activeControls = []
        if self.aneal != 0.0:
            if self.verbose: print "Performing Pre-Optimization Target Selection"
            activeControls = self.reductionControl(model, self.targets, sSynthVal, self.controlLibraries)
        
        if self.verbose: print "-----checking pre-start conditions-----"
        (activeControls,targets,naturalCoefficient,preStart) = self.loadActiveControls(model,targets,limitTargets=True)
        if self.verbose: print "Using Prestart [%s] Active controls [%s]" % (preStart,len(activeControls))
                                                      
        if naturalCoefficient != 0:
            pNaturalCoefficient = naturalCoefficient    
        
        
        '''
        ==========================================================
        Begin Primary Search Loop: Progressive Target Discovery
        ==========================================================
        '''
         
        if self.verbose: print "\n==Starting Iterative Search Size[%s] Iterations [%s->%s]==" %(searchNeighborhood,preStart,preStart+self.iterations)
        if self.verbose: print "Number of Targets [%s]" % len(targets)
        if self.verbose: print "Starting natural coefficient [%s]" %  (naturalCoefficient)
        startTime = time()
            
        for i in range(preStart+1,preStart+self.iterations+1):
                iTime = time()
                
                if self.verbose: print "\n=============Iteration [%s]===============" % i    
                if self.verbose: print "Progressive Growth Parameter: [%s]" % naturalCoefficient
                self.con.objectiveCoeffecent = naturalCoefficient
                
                #--------------------------------------
                # Set targets for optimization
                #--------------------------------------
                
                controlNames = self.getControlNames(targets)
                if self.verbose: print "Reaction Target Size [%s]" % (len(controlNames))
                    
                #---------------------------------------
                # reconstruct optimization model
                #---------------------------------------
                
                tmodel = self.reconstructControlModel(originalModel,controlNames,ltPredV)
                if activeControls != [] and activeControls != None:
                    if self.verbose: print "Activating set targets [%s]" % activeControls
                    self.con.activateControl(tmodel,activeControls,on=True)
                    activeControls = None
                    
                #-----------------------------
                # Run optimization
                #-----------------------------
                if self.verbose: print "\n=====================Optimizing Framework====================="
                
                runTag = "%s_%s_%s_k%s_i%s" % (self.control,self.naturalObjectiveName,self.syntheticObjectiveName,searchNeighborhood,i)
                (lt,sCode) = self.runOptimization(tmodel,runTag)
                iRunTime = time() - iTime
                if self.verbose: print "Optimization Completed: Code [%s] Time [%s]" % (sCode,iRunTime)
                
                if sCode == None:
                    if self.verbose: print "Failure to find optimal solution: breaking search loop"
                    break
                else:
                    if self.verbose: print "--> Optimal solution found"
                pPredValue = ltPredV
                ltPredV = lt.getMipPredictionMap()
                
                lt.clear()
        
                #----------------------------------------------------------
                # Process and print current state of optimization search
                #----------------------------------------------------------
                
                (iObjective,iGeneObjective,ltSynthVal,checkPredVal,checkSynthVal) = self.findCurrentObjective(ltPredV,model,originalModel,controlNames,geneClusters)                
                if abs(ltSynthVal- checkSynthVal) > max(abs(ltSynthVal *0.25),1.0):
                    print "Large model prediction inaccuracy"
                    print "==========breaking search loop!=========="
                    continue

                #--------------------------------------------------------
                # Calculate updated natural parameter 
                #--------------------------------------------------------
                                
                (effectValue,naturalCoefficient) = self.progressiveSearchParameter(ltPredV,pPredValue,oPredVal,sPredVal,iObjective,naturalCoefficient,pNaturalCoefficient,pEffectValue,maxNatural,maxEffect)
                pEffectValue = effectValue        
                
                print "-------------Search loop complete-----------"
                print "Setting search position"

        ''''                 
        ================================================
         End Search Loop: Progressive Target Discovery
        ================================================
        '''
        
        if self.verbose: print "\n=======Search Concluded========="
    
        if finalPredictionValue == None:
            finalPredictionValue = ltPredV
        
        (fObjective,iGeneObjective,ltSynthVal,finalPreditionValue,checkSynthVal) = self.findCurrentObjective(finalPredictionValue,model,originalModel,controlNames,geneClusters)
        
        ltReport = self.con.createReport(finalPredictionValue,originalFluxNames,targets,geneClusters)
        #ltReport["confirm"] = finalCheckValue
        
        return (ltReport,oPredVal,sPredVal,finalCheckValue,fObjective)
    
    def testControls(self,model,controlLibraries,objectives,coeffecent=1.0):
        lp = self._solverFactory()
        result = ''
        for (key,values) in controlLibraries:
            ivalues = {}
            for (k,v) in values.items():
                ivalues[k] = coeffecent*v
            prediction = lp.run(model=model, objective = ivalues)
            for obj in objectives:
                v = prediction[obj]
                result = "%s %s[%s]" % (result,obj,v)
            result = result + "\n"
        return result
        
    def controlSearch(self,model,controlLibrary,objective,targetPercent = 0.50,dir=1.0):     
        lp = self._solverFactory()    
        prediction = lp.run(model=model,objective = controlLibrary)
        baseLine = prediction[objective] * targetPercent
        iControl = {}
        iControl.update(controlLibrary)
        for key in controlLibrary.keys():
            del iControl[key]
            prediction = lp.run(model=model,objective = iControl)
            iFlux = prediction[objective]
            if iFlux*dir < baseLine * dir:
                iControl[key] = controlLibrary[key]
        return iControl
    
    def reduceControlMap(self,values,targets,objective):
        ivalues = {}
        ivalues.update(objective)
        for (k,v) in values.items():
            if k in targets:
                vp = max(v)
                vn = min(v)
                if abs(vp) > abs(vn):
                    vf = vp
                else:
                    vf = vn
                ivalues[k] = vf*-1.0
                
        return ivalues
        
    def adjustControlMap(self,control,dir=-1.0):
        result = {}
        for (k,v) in control.items():
            result[k] = dir*v
        return result
        
    def reductionControl(self, model, targets, syntheticMaxValue, controlLibraries):
            
        if self.verbose: print "===Testing objective control=="
        activeControlObjectives = []
        
        if self.verbose: print "Anealing starting from %s of max" % self.aneal
        
        minVal = syntheticMaxValue * self.aneal
        if self.verbose: print "Original Synthetic Flux %s * min percent %s := min anealing flux %s" % (syntheticMaxValue,self.aneal,minVal)
        
        for (key,prefix,values) in controlLibraries:
            ivalues = self.adjustControlMap(values)
            ivalues.update({self.naturalObjectiveName:-0.02})
            
            if self.syntheticObjectiveName in ivalues.keys():
                del ivalues[self.syntheticObjectiveName]
            
            lo = self._solverFactory()
            iPredVal = lo.run(model=model, objective = ivalues)
            iBioVal = iPredVal[self.naturalObjectiveName]
            iSynthVal = iPredVal[self.syntheticObjectiveName]    
            lo.clear
            
            if self.verbose: print "Objective Library [%s] size [%s]: Biological [%s] Synthetic: [%s]" % (key, len(ivalues), iBioVal, iSynthVal)
            #if self.verbose: print "Control: [%s]" % (ivalues)
            
            if iSynthVal > minVal:
                itargets = iPredVal.keys()
                (minFlux, minControl) = self.controlFactory.reduceControl(model, ivalues, target=self.syntheticObjectiveName, minTarget = minVal, targets= itargets)
                if self.verbose: print "Control reduction to [%s] flux %s [%s]" % (len(minControl),self.syntheticObjectiveName, minFlux[self.syntheticObjectiveName])
                if len(minControl) > 0:
                    if self.verbose: print "Control set %s" % (minControl)
                for (k,v) in minControl.items():
                    minControl[k] = v * -1.0
                activeControlObjectives.append((key,minControl))
                
            else:
                if self.verbose: print "Production level too low for pre-processing"
        
        result = []        
        for (libName,control) in activeControlObjectives:
            if self.rGeneMap != None:
                controlList = self.geneControlList(control, model.controlMap)
            else:
                controlList = control.keys()
            if self.verbose: print "Gene Control [%s] => (%s)[%s]" %(libName,len(controlList),controlList)
            result.append((libName,controlList))
        
        return result
    
    def geneControlList(self,objective,rGeneMap):
        result = set()
        for (key,value) in objective.items():
            if key in rGeneMap.keys():
                values = rGeneMap[key]
                result = result.union(values)
        return result
    
    def writeToLog(self,fileName,logLine,startTime,endTime):
        logFileH = open(fileName,'a')
        logFileH.write(logLine)
        logFileH.close()
        return 0;
    #! should be moved to some kind of string utility package
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
    
    def _annotateGenes(self,objective,model,annotationName,regex):
        result = objective
        if annotationName in model.annotation.keys():
                annotationMap = model.annotation[annotationName]
                gMap = annotationMap.getColumn("gene")
                if annotationMap != None:
                    result = self._annotateMap(objective,gMap,regex)
        return result
    
    def _deannotateGenes(self,objective,model,annotationName,regex):
        result = objective
        if annotationName in model.annotation.keys():
                annotationMap = model.annotation[annotationName]
                gMap = annotationMap.getColumn("gene")
                bMap = {}
                for (key,value) in gMap.items():
                    bMap[value] = key
                if annotationMap != None:
                    result = self._annotateMap(objective,bMap,regex)
        return result
        
    def _annotateObjective(self,objective,model,annotationName,regex):
        aValue = {}
        if self.useGeneMap:
            (iGeneObjective,iOtherObjective) = self.printGeneObjective(objective)
            iGeneObjective = self._annotateGenes(iGeneObjective, model, annotationName, regex) 
            aValue.update(iGeneObjective)
            aValue.update(iOtherObjective)
        return aValue
    
    def _condenseObjective(self,objective):
        result = {}
        for rxn in objective.keys():
            adjust = objective[rxn]
            geneTag = self._getGeneTag(rxn, geneClusters=None)
            if geneTag not in result.keys():
                result[geneTag] = {}
            result[geneTag][rxn] = adjust
        return result
    
    def _controlScore(self,objective,originalPrediction,newPrediction):
        fluxDiff = {}
        k1 = originalPrediction.keys()
        k2 = newPrediction.keys()
        keys = set(k1).union(k2)
        keys = keys.intersection(objective.keys())
        
        for k in keys:
            v1 = 0.0
            v2 = 0.0
            
            if k in k1:
                v1 = originalPrediction[k]
            if k in k2:
                v2 = newPrediction[k]
            
            diff = v2-v1
            fluxDiff[k] = diff
        
        result = {}    
        for key in objective.keys():
            if key in fluxDiff.keys():
                adjust = fluxDiff[key]
                geneTag = self._getGeneTag(key, geneClusters=None)
                if geneTag not in result.keys():
                    result[geneTag] = {}
                result[geneTag][key] = adjust
                    
        return result
    
    def _condenseMap(self,controlMap,fun):
        result = {}
        for (key,value) in controlMap.items():
            v = value.values()
            result[key] = fun(v)
        return result

    def _directionCount(self,map):
        result = {}
        for (key,values) in map.items():
            iresult = {"neg":0,"zero":0,"pos":0}
            for (ikey,value) in values.items():
                fvalue = float(value)
                if fvalue == 0.0:
                    iresult["zero"] += 1
                if fvalue < 0.0:
                    iresult["neg"] += 1
                if fvalue > 0.0:
                    iresult["pos"] += 1
            result[key] = iresult
        return result
                
    def writeAnnotatedTargets(self,objective,model,annotationName,regex,iteration,oPrediction={},nPrediction={},outputFileName=None):
        
        if outputFileName == None:
            outputFileName = "Target_Report_T_%s_C_%s_M_%s_N_%s_S_%s_K_%s_I_%s.txt" % (self.resultTag,self.control,self.modelName,self.naturalObjectiveName,self.syntheticObjectiveName,self.con.controlMax,iteration)
            
        aValue = self._annotateObjective(objective, model, annotationName, regex)
        controlValues = self._condenseObjective(objective)
        controlValues = self._annotateGenes(controlValues, model, annotationName, regex)
        controlScores = self._controlScore(objective,oPrediction, nPrediction)
        controlScores = self._annotateGenes(controlScores, model, annotationName, regex)
        controlScore = self._condenseMap(controlScores,sum)
        scoreCount = self._directionCount(controlScores)
        
        report = Report()
        report.addColumnHash("Control", aValue)
        report.addColumnHash("Flux", controlScore)
        report.addColumnHash("Reaction Control", controlValues)
        report.addColumnHash("Flux Changes", controlScores)
        report.addColumnHash("Flux Count", scoreCount)

        oName = self.resultDirectory + outputFileName
        try:
            writer = ReportWriter()
            writer.setFile(oName)
            writer.write(report) 
            writer.closeFile()
            if self.verbose: print "Report written [%s]" % (outputFileName)
        except:
            if self.verbose: print "Failed to write report [%s] " % (outputFileName)
        
        return None
    
    def writeReport(self,ltReport,fluxModel, oPredVal,sPredVal,neighborhood,iterations,fObjective,resultDir):
        if self.verbose: print "Preparing report ..."        
        report = Report()
        reactionNames = fluxModel.network.getOrderedReactionNames()
        reactionMap = fluxModel.network.getReactionMap()
        
        for reactionName in reactionNames:
            reaction = reactionMap[reactionName]
            equation = reaction.getEquation()
            pathway = reaction.getAnnotation("Subsystem")
            name = reaction.getName()#Currently sensitivity values can be too large for control factors.
            report.addElement(reactionName,"name",name)
            report.addElement(reactionName,"equation", equation)
            report.addElement(reactionName,"Subsystem", pathway)
            
        #--------------------
        #Analysis in report
        #--------------------
        report.addColumnHash("Original", oPredVal)
        report.addColumnHash("Synthetic", sPredVal)
        report.extend(ltReport)    
        report.addColumnHash("Final Objective", fObjective)
    
        #--------------------
        #Write Report
        #--------------------
        if self.verbose: print "Writing report ..."
        outputFileName = resultDir + "RD_" + self.control + "_" + self.modelName + "_" + self.syntheticObjectiveName + "_" + str(neighborhood) + "_" + str(iterations) +  "_" + strftime('%Y%m%dT%H%M%S') + ".txt"  
        
        writer = ReportWriter()
        writer.setFile(outputFileName)
        writer.write(report) 
        writer.closeFile()
        if self.verbose: print "Report Written [%s]" % (outputFileName)
        
        return report
    
    def testOpt(self,tmodel,originalModel,targets,oPrediction=None,testName = "test",fullReport = False):
        
        try:
            (lt,sCode) = self.runOptimization(tmodel)
        except:
            print "----------------optimization test failure-----------------------"
            return None
                    
        ltPredV = lt.getMipPredictionMap()
        ltBioVal = ltPredV[self.naturalObjectiveName]
        ltSynthVal = ltPredV[self.syntheticObjectiveName]
        ltRowStatus = lt.getRowsStatus()
        ltColumnStatus = lt.getColumnsStatus() 
        targets = originalModel.getColumnNames()
        
        if sCode != "OPTIMAL" and sCode != 200:
            if self.verbose: print "Failure to find optimal solution Code: %s" % (sCode)
            dumpFileName = "Test_Optimization_Failure_%s_%s" % (testName,sCode)
            lt.writeMPS(dumpFileName + ".mps")
        else:
            dumpFileName = "Test_Optimization_success_%s_%s" % (testName,sCode)
        
            
        if fullReport:
            
            modelReport = tmodel.modelReport(prediction=oPrediction)
            modelReport.addColumnHash("status",ltRowStatus)
            modelReport.addColumnHash("status",ltColumnStatus)
            writer = ReportWriter()
            writer.setFile(dumpFileName + ".csv")
            writer.write(modelReport) 
            writer.closeFile()
            
        lt.clear()
                
        (iObjective,iLibObjective) = self.con.generateObjective(originalModel, targets, ltPredV)
        if self.verbose: print "Current library source objective[%s]: %s" % (str(len(iLibObjective)),iLibObjective)
        if self.verbose: print "Current objective[%s]: %s" % (str(len(iObjective)),iObjective)
        if self.verbose: print "Current Natural [%s] Synthetic [%s]" % (ltBioVal,ltSynthVal)
        (iGeneObjective,iOtherObjective) = (None,None)
        if self.useGeneMap:
            (iGeneObjective,iOtherObjective) = self.printGeneObjective(iObjective)
            if self.verbose: print "Current Genetic objective[%s]: %s" % (str(len(iGeneObjective)),iGeneObjective)
            if self.verbose: print "Current Other objective: %s" % (iOtherObjective)
        
        return (ltPredV,iObjective)
    
    
    def testDesignOptimization(self,originalModel,fullTest = False):
        #---------------------------------------
        # adjust model
        #---------------------------------------
        controlMin = self.con.controlMin
        controlMax = self.con.controlMax
        targets = originalModel.targets
        
        self.con.controlMin = 0
        self.con.controlMax = 0
        
        controlNames = self.getControlNames(self.targets)
        tmodel = self.reconstructControlModel(originalModel,controlNames,{})        
        #-----------------------------
        # Run optimization
        #-----------------------------
        
        #No controls active:
        if self.verbose: print "Testing optimization framework with no controls in on state"
        objective = self.testOpt(tmodel, originalModel, targets)
        if objective == None:
            return None
        
        #short version of test
        if not fullTest:
            return True
        
        #Try each control
        for (name,prefix) in self.con.controlPrefixMap.items():
            
            if self.useGeneMap:
                controlMap = self.con.getLibraryGeneControlMap(name)
            else:
                controlMap = self.con.regulationLibrary[name]
            
            if controlMap == {}:
                continue
            self.useScip = True
            
            self.con.controlMin = 0
            self.con.controlMax = 1.0
            
            tmodel = self.reconstructControlModel(originalModel,targets,{})
            (predV,obj) = self.testOpt(tmodel, originalModel, targets)
            
            controlNames = controlMap.keys()            
            for controlName in controlNames:
                controlValues = controlMap[controlName]
                
                if self.verbose: print "\nConstructing optimization check [%s] %s => %s" % (name,controlName,controlValues)
                tmodel = self.con.setControl(tmodel,prefix + "y__",[controlName],on = True, binary = True)        
                self.testOpt(tmodel, originalModel,targets,testName = prefix+controlName+"_boundary")
                tmodel = self.con.setControl(tmodel,prefix + "y__",[controlName],on = False, binary = True)
        
            self.con.controlMin = 0
            self.con.controlMax = 0    
            tmodel = self.reconstructControlModel(originalModel,targets,{})
            
            for controlName in controlMap.keys():    
                controlValues = controlMap[controlName]
            
                if self.verbose: print "\nConstructing optimization check [%s] %s => %s" % (name,controlName,controlValues)                    
                tmodel = self.con.setControl(tmodel,prefix + "y__",[controlName],on = True)        
                self.testOpt(tmodel, originalModel, targets,testName = prefix+controlName+"_active")
                tmodel = self.con.setControl(tmodel,prefix + "y__",[controlName],on = False)
                
        self.con.controlMin = controlMin
        self.con.controlMax = controlMax

        return None
