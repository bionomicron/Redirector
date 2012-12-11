#Graham Rockwell
#7 7 2009
#Module for analysis of pathway design data

from core.reader.FlatFileParser import FluxModelFlatFileParser,FlatFileParser, AnnotationFlatFileParser
from core.model.ModelFactory import ModelFactory
from core.model.LPSolver import LPSolver
from core.model.ConstructRegulationOptimization import ConstructRegulationOptimization
from core.model.OptimizationControlRedirector import OptimizationControlRedirector
from core.model.ProcessModel import ProcessGeneticObjective, ProcessModel, MatrixTools, ModelProperties

from core.util.Config import ReflectionConfig
from core.util.Report import Report
from core.util.ReportWriter import ReportWriter
from math import sqrt
from itertools import combinations
from optparse import OptionParser
import pickle, os

#from pydot import Dot, Node, Edge, Cluster

def getPreStart(fileName):
    
    if os.path.exists(fileName):
        fileHandle = open(fileName,"r")
        (constructor,objective,prediction) = pickle.load(fileHandle)
        fileHandle.close()
    else:
        return None
    return (constructor,objective,prediction)

def parseUserTargets(controlFile,model):
    parser = FlatFileParser()
    gtMap = model.getGeneTargetMap()
    report = parser.parseToReport(controlFile, "GeneName", header=["GeneName","library","prefix","bnumber","value"])
    rowNames = report.returnRowNames()
    controlNames = []
    controlIds = set()
    controlLibrariesMap = {}
    
    for rowName in rowNames:

        libName = report.getElement(rowName, "library")
        prefix = report.getElement(rowName, "prefix")
        geneID = report.getElement(rowName, "bnumber")
        value = report.getElement(rowName, "value")
        
        rxnControlMap = {}
        geneID = geneID.replace(" ","_")
        if geneID not in gtMap.keys():
            continue
        for rxn in gtMap[geneID]:
            rxnControlMap[rxn] = value
        
        controlID = (libName,prefix,geneID)
        
        controlNames.append(geneID)
        controlIds.add(controlID)
    
        if (libName,prefix) not in controlLibrariesMap.keys():
                controlLibrariesMap[(libName,prefix)] = {}
        controlLibrariesMap[(libName,prefix)].update(rxnControlMap)
    
    controlLibraries =[]
    for (key,value) in controlLibrariesMap.items():
        (libName,prefix) = key
        controlLibraries.append((libName,prefix,value))
    
    return (controlNames,controlIds,controlLibraries)
    
def parseUserControl(controlFile,model,config,minObjVal):
    
    (controlNames,controlIds,controlLibraries) = parseUserTargets(controlFile,model) 
    
    con = ConstructRegulationOptimization()
    config.reflect("Redirector Model",con)
    con.delta = 1e-4
    con.minObjective = minObjVal
    con.controlMax = options.searchNeighborhood
    con.newObjective = syntheticObjectiveName
    
    con.addRegulationLibraries(controlLibraries)
    if modelMatrix.controlMap != None: 
        con.setGeneMap(modelMatrix.controlMap)
    
    return (con,controlIds)

def parseControlSave(originalModel,controlTags,prefix,k,n,dir='',orderGenes=False):
    '''
    Parse files for control targets.
    '''
    
    con = None
    controlList = list()
    result = set()
    names = []
    
    for controlTag in controlTags:
        for zk in k:
            for zn in n:
                iname =  prefix % (controlTag,zk,zn)
                names.append(iname)
        
    print "Parsing [%s] persisted control optimizations" % (len(names))
    finalConstructor = ConstructRegulationOptimization()
    
    for name in names:
        iname = dir + name
        persist = getPreStart(iname)
        if persist == None:
            continue
        else:
            print "Parsed: [%s]" %(iname)
        (con,iObjective,predictions) = persist
        targets =  iObjective.keys()
        controlNames = con.getControlNames(targets,predictions)
        controlList.append(controlNames)
        result = result.union(controlNames)
        finalConstructor.extendControl(con)
        
    return (finalConstructor,controlList,result)

def simpleTestOpt(con,originalModel,controlTargets,naturalObjectiveName,syntheticObjectiveName,verbose = False,naturalCoefficient=-0.05):
    '''
    Optimiztaion Tester for depdenency anlaysis
    Construciton object contains union of all libraries which are used by target sets.
    
    '''
    gRMap = con.getGeneTargetMap()
    regulationLibrary = con.regulationLibrary
    objective = {}
    
    for (lib,prefix,target) in controlTargets:
        if gRMap != None:
            if target in gRMap.keys():
                targets = gRMap[target]
            else:
                print "target %s not found in gene map" % (target)
                targets = [target]
        else:
            targets = [target]
        for target in targets:
            if lib in regulationLibrary.keys():
                if target in con.regulationLibrary[lib].keys():
                    controlValue = con.regulationLibrary[lib][target]
                    objective[target] = -1.0 * controlValue
    
    iObjective = dict(objective)
    objective[naturalObjectiveName] = naturalCoefficient
            
    lp = LPSolver()
    prediction = lp.run(originalModel, objective= objective)

    rObjective = iObjective
    if verbose: print "Current objective[%s]: %s" % (str(len(iObjective)),iObjective)
    (iGeneObjective,iOtherObjective) = (None,None)
    if originalModel.controlMap != None:
        (iGeneObjective,iOtherObjective) = originalModel.printGeneObjective(iObjective)
        if verbose: print "Current Genetic objective[%s]: %s" % (str(len(iGeneObjective)),iGeneObjective)
        if verbose: print "Current Other objective: %s" % (iOtherObjective)
        rObjective = iGeneObjective
        
    lp.clear()
    del lp
    
    return (prediction,rObjective)
    
            
def controlThreshold(values,exclude):
    result = {}
    for v in values:
        for (e,level) in exclude:
            e = set(e)
            if len(e.difference(v)) == 0:
                if v in result.keys():
                    ilevel = result[v]
                    if level > ilevel:
                        result[v] = level
                else:
                    result[v] = level  
    return result

def loadPersist(tag):
    if os.path.exists(tag):
        file = open(tag,"r")
        value = pickle.load(file)
        file.close()
        return value
    else:
        return None

def persist(tag,value):
    persistFile = open(tag,"w")
    pickle.dump(value,persistFile)
    persistFile.close()
    #print "persisting dependency analysis"
    return tag
        
def findControlDependencies(con,originalModel,controls,exclusion,options,searchSize=0,targetPrecent = 0.40,verbose = False,usePersist=False):
    name = "".join(originalModel.modelName)
    naturalObjectiveName  = options.bioObj
    syntheticObjectiveName = options.synthObj
    controlTag = options.control
    
    persistTag = "Dependency_persist_%s_%s_%s_%s_%s_%s" % (controlTag,name,naturalObjectiveName,syntheticObjectiveName,str(searchSize),str(targetPercent))
    delta = 1e-4
    iter = 0
    iterSave = 100
    
    completed = set()
    controlSet = set()
    controlMap = []
    report = Report()
    
    persistValue = None
    if usePersist:
        persistValue = loadPersist(persistTag)
        if persistValue != None:
            print "loading persisted dependency search"
            (completed,controlMap,controlSet,report) = persistValue 
    
    controlSubSets = combinations(controls,searchSize)
    xcontrolSubSets = set(controlSubSets)
    xcontrolSubSets = xcontrolSubSets.difference(completed)
    
    print "Control dependency depth [%s] search sets [%s]" % (searchSize,len(xcontrolSubSets))
            
    lo = LPSolver()
    sObjective = {syntheticObjectiveName:-1.0}
    iPred = lo.run(model= originalModel, objective = sObjective)
    synthVal = iPred[syntheticObjectiveName]
    targetSVal = synthVal * targetPrecent
        
    #print "Synthetic Target [%s] = %s * %s" % (targetSVal,synthVal,targetPercent)
        
    controlMin = con.controlMin
    controlMax = con.controlMax
        
    con.controlMin = 0
    con.controlMax = 0
        
    targets = []
    gRMap = con.getGeneTargetMap()
        
    for (lib,prefix,target) in controls:
        if gRMap != None:
            if target in gRMap.keys():
                rTargets = gRMap[target]
                targets.extend(rTargets)
        else:
            targets.append(target)        
        
    #-----------------------------
    # Run optimization
    #-----------------------------
    iter = 0
    for controlSub in xcontrolSubSets:
        iter += 1
        completed.add(controlSub)
        
        #Analysis of control
        (prediction,objective) = simpleTestOpt(con, originalModel, controlSub, naturalObjectiveName,syntheticObjectiveName,verbose=False,naturalCoefficient=-0.5)
                
        #Save results        
        syntheticFlux = prediction[syntheticObjectiveName]
        #if verbose: print "Objective [%s]" % (objective)
        #if verbose: print "Flux [%s]" % (prediction)
            
        if syntheticFlux > targetSVal:    
            controlThresholds = controlThreshold([controlSub], exclusion)
            if controlSub in controlThresholds.keys():
                if syntheticFlux <= controlThresholds[controlSub] + delta:
                    #if verbose: print "no improvement"
                    continue
                
            if verbose: print "Synthetic Production [%s]" % (syntheticFlux)    
            
            controlResultMap = originalModel.annotateGenes(objective,annotationName = "bnumber", regex="[a-zA-Z0-9\(\)]+")
            controlTag = str(controlResultMap.keys())                
            
            controlSV = (controlSub,syntheticFlux)
            controlSet.add(controlSV)
            controlSV = (controlResultMap,syntheticFlux)
            controlMap.append(controlSV)
            
            #if verbose: print "Control Set [%s]" % (controlSub)
            if verbose: print "Control Map [%s]" % (controlResultMap)
                
            report.add(controlTag, "control values", controlResultMap)
            report.add(controlTag, "production", syntheticFlux)
                
        if iter % iterSave == 0:
            print "persisting iteration %s" % (iter)
            persistValue = (completed, controlSet,controlMap, report)
            persist(persistTag, persistValue)
                     
    con.controlMin = controlMin
    con.controlMax = controlMax
        
    #Persist dependency search
    if usePersist:
        persistValue = (completed, controlSet,controlMap, report)
        persist(persistTag, persistValue)

    return (controlMap,controlSet,report)
    
def scanControlDependencies(con,originalModel,controls,options,searchRange,targetPrecent = 0.80,verbose = False):
    report = Report()
    controlMaps = []
    controlSets = set()
    for z in searchRange:
        (controlMap,controlSet,iReport) = findControlDependencies(con,originalModel,controls,controlSets,options,z,targetPrecent,verbose)
        controlSets = controlSets.union(controlSet)
        controlMaps.extend(controlMap) 
        report.extend(iReport)
    
    return (controlMaps,report)

def productionLevels(originalModel,options,minBio=0.20,verbose=False):
    naturalObjectiveName  = options.bioObj
    syntheticObjectiveName = options.synthObj
    
    cellularObjective = {naturalObjectiveName:-1.0}
    syntheticObjective = {syntheticObjectiveName:-1.0}
    
    lo = LPSolver()
    lo.setModel(originalModel)
    oPredVal = lo.run(objective = cellularObjective)
    oBioVal = oPredVal[naturalObjectiveName]
    oSynthVal = oPredVal[syntheticObjectiveName]

    if verbose: print "Max cellular: Biological [%s] Synthetic: [%s]" % (oBioVal, oSynthVal)
    
    #---------------------------------------
    # Values for Max Synthetic Objective
    #---------------------------------------
    
    sPredVal = lo.run(objective = syntheticObjective)
    sBioVal = sPredVal[naturalObjectiveName]
    sSynthVal = sPredVal[syntheticObjectiveName]
    lo.clear()
    
    if verbose: print "Max synthetic: Biological [%s] Synthetic: [%s]" % (sBioVal, sSynthVal)

    #-----------------------------------
    # Values for minimum Natural objective
    #-----------------------------------
    
    minObjVal = oBioVal*0.20
    originalModel.addColumnLimit(naturalObjectiveName,(minObjVal,None))
    lo = LPSolver()
    s2PredVal = lo.run(model=originalModel,objective = syntheticObjective)
    s2BioVal = s2PredVal[naturalObjectiveName]
    s2SynthVal = s2PredVal[syntheticObjectiveName]
    
    lo.clear()
    
    return (oPredVal,sPredVal,s2PredVal)

def controlDependencies(originalModel,resultTag,options,targetPercent,searchDepth=3,verbose=False,usePersist=True):
    naturalObjectiveName  = options.bioObj
    syntheticObjectiveName = options.synthObj
    k = options.searchNeighborhood
    n = options.iterations
    
    '''
    * Check for persisted control dependencies
    '''
    
    persistTag= "C_%s_M_%s_N_%s_S_%s_k%s_n%s_s%s" % (options.control,options.modelName,naturalObjectiveName,syntheticObjectiveName,k,n,str(searchDepth))
    persistFileName = "ControlSearches/Persist_Cdependency_%s.bac" % persistTag    
    if os.path.exists(persistFileName) and usePersist:
        persistFile = open(persistFileName,"r")
        (dependencies,dReport) = pickle.load(persistFile)
        persistFile.close()
        return (dependencies,dReport)
    
    '''
    * Find synthetic flux value and calculate minimum flux threshold 
    '''
    
    (oPredVal,sPredVal,s2PredVal) = productionLevels(originalModel, options, minBio=0.20)
    s2BioVal = s2PredVal[naturalObjectiveName]
    s2SynthVal = s2PredVal[syntheticObjectiveName]
    minObjVal = s2SynthVal*targetPercent
    
    '''
    Retrieve targets from persisted designs or from flat file
    #!Fix this section to allow for user input and other input in better command line input.
    '''
    
    if verbose: print "Min Biological Objective Max Synthetic: Biological [%s], Synthetic [%s]" % (s2BioVal,s2SynthVal)
    controlFile = options.useControl
    if controlFile == "T":
        controlFile = "ExpControlTargets_%s.csv" % (syntheticObjectiveName)
    elif controlFile != '':
        (con,targets) = parseUserControl(controlFile,originalModel,config,minObjVal)
        targetList = [targets]
    else:
        for controlTag in options.control.split(","):    
            print "loading control %s" % (controlTag)
        
        searchRange = range(1,searchDepth+1)
        prefix = "preStart_%s_%s_%s_%s_%s_%s" % ("%s",options.modelName,naturalObjectiveName,syntheticObjectiveName,"%s","%s")
        kRange = range(1,k+1)
        nRange = range(0,n+1)
            
        print "Search depth [%s], control optimizations for k[1->%s] and n[0->%s]" % (str(searchDepth),k,n)
        controlTags = options.control.split(",")
        (con,targetsList,targets) = parseControlSave(originalModel,controlTags,prefix,kRange,nRange, dir="ControlSearches/")
        
    maxTargetSize = 0 
    for targetset in targetsList:
        tsize = len(targetset)
        if tsize > maxTargetSize:
            maxTargetSize = tsize
    totalTargets = len(targets)
        
    print "Largest Design Set [%s], Total targets [%s]" % (maxTargetSize,totalTargets)
    print "Performing dependency analysis on [%s] t
        result = result.union(controlNames)
        finalConstructor.extendControl(con)
        
    return (finalConstructor,controlList,result)

def simpleTestOpt(con,originalModel,controlTargets,naturalObjectiveName,syntheticObjectiveName,verbose = False,naturalCoefficient=-0.05):
    '''
    Optimiztaion Tester for depdenency anlaysis
    Construciton object contains union of all libraries which are used by target sets.
    
    '''
    gRMap = con.getGeneTargetMap()
    regulationLibrary = con.regulationLibrary
    objective = {}
    
    for (lib,prefix,target) in controlTargets:
        if gRMap != None:
            if target in gRMap.keys():
                targets = gRMap[target]
            else:
                print "target %s not found in gene map" % (target)
                targets = [target]
        else:
            targets = [target]
        for target in targets:
            if lib in regulationLibrary.keys():
                if target in con.regulationLibrary[lib].keys():
                    controlValue = con.regulationLibrary[lib][target]
                    objective[target] = -1.0 * controlValue
    
    iObjective = dict(objective)
    objective[naturalObjectiveName] = naturalCoefficient
            
    lp = LPSolver()
    prediction = lp.run(originalModel, objective= objective)

    rObjective = iObjective
    if verbose: print "Current objective[%s]: %s" % (str(len(iObjective)),iObjective)
    (iGeneObjective,iOtherObjective) = (None,None)
    if originalModel.controlMap != None:
        (iGeneObjective,iOtherObjective) = originalModel.printGeneObjective(iObjective)
        if verbose: print "Current Genetic objective[%s]: %s" % (str(len(iGeneObjective)),iGeneObjective)
        if verbose: print "Current Other objective: %s" % (iOtherObjective)
        rObjective = iGeneObjective
        
    lp.clear()
    del lp
    
    return (prediction,rObjective)
    
            
def controlThreshold(values,exclude):
    result = {}
    for v in values:
        for (e,level) in exclude:
            e = set(e)
            if len(e.difference(v)) == 0:
                if v in result.keys():
                    ilevel = result[v]
                    if level > ilevel:
                        result[v] = level
                else:
                    result[v] = level  
    return result

def loadPersist(tag):
    if os.path.exists(tag):
        file = open(tag,"r")
        value = pickle.load(file)
        file.close()
        return value
    else:
        return None

def persist(tag,value):
    persistFile = open(tag,"w")
    pickle.dump(value,persistFile)
    persistFile.close()
    #print "persisting dependency analysis"
    return tag
        
def findControlDependencies(con,originalModel,controls,exclusion,options,searchSize=0,targetPrecent = 0.40,verbose = False,usePersist=False):
    name = "".join(originalModel.modelName)
    naturalObjectiveName  = options.bioObj
    syntheticObjectiveName = options.synthObj
    controlTag = options.control
    
    persistTag = "Dependency_persist_%s_%s_%s_%s_%s_%s" % (controlTag,name,naturalObjectiveName,syntheticObjectiveName,str(searchSize),str(targetPercent))
    delta = 1e-4
    iter = 0
    iterSave = 100
    
    completed = set()
    controlSet = set()
    controlMap = []
    report = Report()
    
    persistValue = None
    if usePersist:
        persistValue = loadPersist(persistTag)
        if persistValue != None:
            print "loading persisted dependency search"
            (completed,controlMap,controlSet,report) = persistValue 
    
    controlSubSets = combinations(controls,searchSize)
    xcontrolSubSets = set(controlSubSets)
    xcontrolSubSets = xcontrolSubSets.difference(completed)
    
    print "Control dependency depth [%s] search sets [%s]" % (searchSize,len(xcontrolSubSets))
            argets" % (len(targets))
   
    (dependencyMap,dReport) = scanControlDependencies(con, modelMatrix, targets, options, searchRange,targetPercent,verbose)
    persistFile = open(persistFileName,"w")
    pickle.dump((dependencyMap,dReport),persistFile)
    persistFile.close()
        
    return (dependencyMap,dReport)

def cleanNodeNames(nodes):
    keys = nodes.keys()
    keys.sort()
    used = set()
    result = {}
    for key in keys:
        targetNames = nodes[key]
        targetNames = targetNames.difference(used)
        used = used.union(targetNames)
        result[key] = targetNames
    return result

def loadControlTable(fileName,keyTag="Row Names",controlTag="control values",productionTag = "production"):
    '''
    Load flat file of control designs and production values
    '''
    
    if not os.path.exists(fileName):
        print "file name not found %s" % (fileName)
        return ([],None)
        
    header = [keyTag,controlTag,productionTag]
    parser = FlatFileParser()
    report = parser.parseToReport(fileName, keyTag=keyTag, header=header)
    result = []
    for rowName in report.returnRowNames():
        controlMapString = report.get(rowName, controlTag)
        controlMap = eval(controlMapString)
        production = report.get(rowName, productionTag)
        controlMapValue = (controlMap,production)
        result.append(controlMapValue)
    return (result,report)

def loadTargetTable(fileName,keyTag="Row Names",controlTag="control values"):
    '''
    Load flat file of control targets
    '''
    if not os.path.exists(fileName):
        #print "file name not found %s" % (fileName)
        return ({},None)
    else:
        pass
        #print "parsing %s" % (fileName)
    
    header = [keyTag,controlTag]
    parser = FlatFileParser()
    report = parser.parseToReport(fileName, keyTag=keyTag, header=header)
    controlMap = {}
    for rowName in report.returnRowNames():
        keyName= rowName
        controlValue = report.get(rowName, controlTag)
        controlMap[keyName] = controlValue
    
    return (controlMap,report)

def checkControlMap(controlMap,filter):
    for (key,value) in controlMap.items():
        if not filter(value):
            return False
    else: return True

def loadControlTables(fileNames,keyTag,controlTag,productionTag,filter):
    controlMaps = []
    for fileName in fileNames:
        (icontrolMap,report) = loadTargetTable(fileName, keyTag, controlTag)
        if len(icontrolMap) != 0 and checkControlMap(icontrolMap,filter):
            #print "control design: %s" % (icontrolMap)
            controlMaps.append(icontrolMap)     
    return controlMaps

def loadControlTableRange(fileNameTag,controls,sObjective,k,n,keyTag,controlTag,productionTag,filter):
    kRange = range(1,k+1)
    nRange = range(0,n+1)
    fileNames = []
    
    for ki in kRange:
        for ni in nRange:
            for ci in controls:
                iname = fileNameTag % (ci,sObjective,ki,ni)
                fileNames.append(iname)
        
    controlMaps = loadControlTables(fileNames,keyTag,controlTag,productionTag,filter)

    return controlMaps
            
def processControlMaps(controlMaps,minProduction=0.0,bestTargets=[]):
    uniqueTargets = set()
    targetValues = {}
    maxTargetSize = 0
    bestScore = 0
    bestMaps = []
    
    for controlMap in controlMaps:
        iscore = 0
        productionValue = 1 #! (incase method is updated to check production levels of control maps.
        
        if productionValue > minProduction:
            targets = controlMap.keys()
            uniqueTargets = uniqueTargets.union(targets)
            size = len(targets)
            
            if size > maxTargetSize:
                maxTargetSize = size
            
            for (key,value) in controlMap.items():
                if key not in targetValues:
                    targetValues[key] = set()
                targetValues[key].add(value)
                
            
                for k1 in bestTargets:
                    if k1 in key:
                        iscore += 1
                
            if iscore == bestScore and bestScore > 0:
                bestMaps.append(controlMap)
                
            if iscore > bestScore:
                bestScore = iscore
                bestMaps = [controlMap]
                
    bestMaps = (bestMaps,bestScore)
    return (uniqueTargets,targetValues,maxTargetSize,bestMaps)

def loadVariousControlTable(fileNameTag,controls,sObjectiveTags,k,n,filter=lambda x: True):
    '''
    Parse control table files
    '''
    keyTag = "Row Names"
    controlTag = "Control"
    productionTag = ''
    #bestTargets = 'accA,fadD,fadE'.split(",")
    #bestTargets = 'aceE,fabH,fabD,accA,fabA,fabZ,fabB,fabF,fadA,fadI,fadB,fadD'.split(",")
    #mTargets = "pgi,pfkA,eno,fbaA,tpiA,gapA,pgk,gpmA".split(",")
    #bestTargets.extend(mTargets)
    bestTargets = 'accA,lpdA,aceE,aceF,gapA,pgk,fumB,fumC,mdh,acnA,scpC,sucC,sucD,sdhA,sdhB'.split(",")
    #bestTargets = 'lpdA,aceE,aceF,'.split(",")
    
    summaryReport = Report()
    targetReport = Report()
    
    print "best targets %s" % (bestTargets)
    for sObjective in sObjectiveTags:
        controlMaps = loadControlTableRange(fileNameTag, controls, sObjective, k, n, keyTag, controlTag, productionTag,filter)
        (uniqueTargets,targetValues,maxTargetSize,(bestMaps,score)) = processControlMaps(controlMaps, minProduction=0.0,bestTargets=bestTargets)
        print "Production %s => %s" % (sObjective,score)
        
        for map in bestMaps:
            print "%s" % (map)
        uTargetSize = len(uniqueTargets)
        summaryReport.add(sObjective, "Unique Target Number", uTargetSize)
        summaryReport.add(sObjective, "Largest Design", maxTargetSize)
        
        for (key,value) in targetValues.items():
            targetReport.add(key, sObjective, value)
        
    return (summaryReport,targetReport)
        
def buildDependencyControlGraph(dependencies,sizes=[1,2],productionTarget=1.25):
    from pydot import Dot, Node, Edge, Cluster
    
    dgraph = Dot(graph_type='graph',fontname='Verdana',splines="line",maxiter="500")
    #splines="line",
    nodes = {}
    edges = {}
    
    for (controlMap,pvalue) in dependencies:
        targetNames = []
        size = len(controlMap)
        if size not in sizes:
            continue
        for (targetName,cvalue) in controlMap.items():
            #targetTag = "%s(%s)" % (targetName,cvalue)
            targetTag = targetName
            targetNames.append(targetTag)
            if size not in nodes.keys():
                nodes[size] = set()
            nodes[size].add((targetName,cvalue))            
        edgeCombos = combinations(targetNames,2)
        for value in edgeCombos:
            key = list(value)
            key.sort()
            (e1,e2) = key
            if (e1,e2) not in edges.keys() and (e2,e1) not in edges.keys():
                edge = Edge(e1,e2)
                edges[(e1,e2)] = edge
            else:
                #print "dup key %s %s" % (e1,e2)
                pass
            
    nodes = cleanNodeNames(nodes)
    
    for (key,nodeValues) in nodes.items():
        if key == 1:
            ishape = "rectangle"
        elif key == 2:
            ishape = "oval"
        elif key == 3:
            ishape = "hexagon"
        else:
            ishape = "circle"
        
        for (name,value) in nodeValues:
            icolor = "black"
            if value > 0:
                icolor = "red"
            if value < 0:
                icolor = "green"
            targetTag = "%s(%s)" % (name,value)
            targetTag = name       
            dgraph.add_node(Node(targetTag,shape=ishape,color=icolor))
            #print "Node: [%s] : size[%s] value[%s] shape[%s]" % (name,key,value,ishape)
        
    for (key,edge) in edges.items():
        dgraph.add_edge(edge)
    
    return dgraph 
        

def findBestTargetSets(targetFileObjectives,k,n):        
    targetdir = "../../results/"
    #reportTag = "Target_Report_M_%s,_N_%s_S_%s_K_%s_I_%s.txt"
    modelTag = "iAF1260,iAF1260_FattyAcids"
    #naturalTag = "Biomass"
    controls = ["flat","cp","T_cp_C_sense","T_cp_C_flat","T_cp_binary"]
    #controls = ["flat","cp","T_cp_C_flat"]
    #controlFilter = lambda x : x in [1.0,-1.0]
    controlFilter = lambda x: True
            
    #fileNameTag = targetdir + "Target_Report_T_flat_C_flat_iAF1260,iAF1260_FattyAcids_N_Biomass_S_%s_K_%s_I_%s.txt"
    fileNameTag = targetdir + "Target_Report_%s_M_iAF1260,iAF1260_FattyAcids_N_Biomass_S_%s_K_%s_I_%s.txt"

    sObjectiveTags = targetFileObjectives.split(",")
    (targetSummaryReport,targetProductReport) = loadVariousControlTable(fileNameTag, controls, sObjectiveTags, k, n, filter = controlFilter)
        
    targetReportFileName = "ControlSearches/Control_Target_Summary_Report_%s_%s_%s_%s_k%s_n%s.csv" % (resultTag,options.modelName,naturalObjectiveName,syntheticObjectiveName,k,n)
    print "Writing target report %s" % targetReportFileName
    writer = ReportWriter()
    writer.setFile(targetReportFileName)
    writer.write(targetSummaryReport) 
    writer.closeFile()
        
    targetProductReportFileName = "ControlSearches/Control_Target_Product_Report_%s_%s_%s_%s_k%s_n%s.csv" % (controls,modelTag,naturalObjectiveName,syntheticObjectiveName,k,n)
    print "Writing target product report %s" % targetProductReportFileName
    writer = ReportWriter()
    writer.setFile(targetProductReportFileName)
    writer.write(targetProductReport) 
    writer.closeFile()

    return None

if __name__  == "__main__":
    parser = OptionParser()
    
    parser.add_option("-v", "--verbose",
                      action="store_true", 
                      dest="verbose", 
                      default=False,
                      help="print status messages to stdout")
     
    parser.add_option("-c", "--configfile", 
                      dest="configFileName", 
                      default="model_config.txt",
                      help="configuration file", 
                      metavar="FILE")
    
    parser.add_option("-m","--modelname", 
                      dest="modelName", 
                      default="",
                      help="name of model from configuration file", 
                      metavar="String")
    
    parser.add_option("--section",
                      dest="subSections",
                      default = '',
                      help = "Comma separated list of sections of the model files to use")
    
    parser.add_option("-n","--configNames",
                      type = "string",
                      dest = "configNames",
                      default = "Default",
                      help = "comma separated list of configuration settings to include",
                      )  
    
    parser.add_option("-b","--bioObjective", 
                      type = "string",
                      dest="bioObj", 
                      default = 'Biomass',
                      help="name of biological objective reaction", 
                      metavar="String")
    
    parser.add_option("-s","--synthObjective",
                      type = "string",
                      dest="synthObj", 
                      default='',
                      help="name of synthetic objective reaction", 
                      metavar="String")
    
    parser.add_option("--sn", "--searchNeighborhood",
                      dest = "searchNeighborhood",
                      default = 1,
                      type = int,
                      help = "size of search neighborhood",
                      metavar = "int")
    
    parser.add_option("--iter", "--searchIterations",
                      dest = "iterations",
                      default = 1,
                      type = int,
                      help = "maximum number of iterations, if blank no max",
                      metavar = "int")
    
    parser.add_option("-r", "--reduced",
                      action="store_true", 
                      dest="reduce", 
                      default=False,
                      help="use reduced model for analysis")
    
    parser.add_option("--depth",
                      dest="searchDepth",
                      default=3,
                      help="search depth",
                      metavar="int")
                  
    parser.add_option("-a", "--analysisfile", 
                      dest="analysisFile", 
                      default="",
                      help="Tab delimited flat file of analysis", 
                      metavar="FILE")
    
    parser.add_option("--useControl", 
                      dest="useControl", 
                      default='',
                      help="Tab delimited flat file of use control library and targets", 
                      metavar="FILE")

    parser.add_option("--control", 
                      dest="control", 
                      default='',
                      help="Comma separated list of control tags to use", 
                      metavar="FILE")

    parser.add_option("--loadFlatDep", 
                      dest="loadFlatDep", 
                      default='',
                      help="tab delimited file of dependency groups", 
                      metavar="FILE")

    parser.add_option("--loadTargetFiles", 
                      dest="loadTargetFiles", 
                      action="store_true", 
                      default=False,
                      help="load target report files for dependency", 
                      metavar="FILE")

    parser.add_option("--makeGraph", 
                      dest="makeGraph", 
                      action="store_true",
                      default=False,
                      help="turn on undirected graph making", 
                      metavar="FILE")

    parser.add_option("--targetFileObjectives", 
                      dest="targetFileObjectives",
                      default = "",
                      help="List of objectives to get target count information for", 
                      metavar="FILE")
    
    parser.add_option("-o","--outputfile", 
                      dest="outputFileName",
                      default = "rd_output.csv",
                      help="name of report FILE", 
                      metavar="FILE")
    
    # Parse options
    (options,args) = parser.parse_args()    
    
    # Parse main configuration file
    config = ReflectionConfig()    
    config.readfp(open("Redirector.config"))

    #---------------------------
    # configure preset analysis
    #---------------------------
    
    #options.configNames = "Default"
    #options.configNames = "iAF1260 Export"
    #options.configNames = "Test iAF1260"
    #options.configNames = "Test Simple Model"
    #options.configNames = "iAF1260 Export"
    
    configNames = options.configNames.split(",")
    configNames.append("Redirector Model")
    for name in configNames:
        config.merge(name,"Redirector",append=True)
    
    #----------------------------------------
    # reflect options from configuration
    #----------------------------------------

    config.reflect("Redirector",options)
    config.load("Redirector", options.__dict__,override=False)
    
    #-----------------------------
    # Build Storage Directories
    #-----------------------------
    
    if not os.path.exists("ControlLibraries"):
            os.makedirs("ControlLibraries")
    if not os.path.exists("ControlSearches"):
            os.makedirs("ControlSearches")
    
    #----------------------------
    # Parse Inputs
    #---------------------------

    verbose        = options.verbose 
    modelName      = options.modelName
    naturalObjectiveName  = options.bioObj
    syntheticObjectiveName = options.synthObj
    searchDepth = int(options.searchDepth)
    targetFileObjectives = options.targetFileObjectives
    k = options.searchNeighborhood
    n = options.iterations
    targetPercent = 0.40
    usePersist=False
    
    #options.userControl = "T"
    #options.loadTargetFiles = True
    #options.controlFile = "T"
    
    '''
    #----------------------------------------------------
    # Initialize and set values for tools and tfactories
    #----------------------------------------------------
    '''
    verbose = True
    naturalObjective = {naturalObjectiveName:-1.0}
    syntheticObjective = {syntheticObjectiveName:-1.0}
    resultTag = config.get("Redirector", "resultTag")
    protectedTargets = set(naturalObjective.keys()).union(syntheticObjective.keys())
    
    if verbose: print "Process Analysis Version 1.0"
    if verbose: print "Model names: [%s]" % (modelName)
    if verbose: print "Synthetic objective: [%s]" % (syntheticObjectiveName)
    if verbose: print "Parsing data files for [%s]" % (modelName)
    
    '''
    * 0. Process target files into overview table
    '''
    targetFileObjectives = options.targetFileObjectives
    if targetFileObjectives != '':
        findBestTargetSets(targetFileObjectives, k, n)      
    
    '''
    I. Parse data files and configuration settings
    '''
    
    modelNames = modelName.split(",")
            
    modelFactory = ModelFactory()
    config.reflect("Redirector",modelFactory)
    modelFactory.protectedTargets = protectedTargets
    
    if verbose: print "----------------Loading Metabolic Models---------------"
    (fluxModel,modelMatrix,reducer,geneReduction) = modelFactory.loadModel(modelNames)    
    
    targets = modelMatrix.targets
    
    persistTag= "%s_%s_%s_%s_k%s_n%s_s%s" % (resultTag,options.modelName,naturalObjectiveName,syntheticObjectiveName,k,n,str(searchDepth))
    persistFileName = "ControlSearches/Persist_Cdependency_%s.bac" % persistTag
    
    '''
    II. Dependency analysis
    
    if option search depth != 0:
        Perform dependency analysis
    '''
    '''
    Scan persisted designs
    or
    Load controls from previously created flat file
    '''
    
    if searchDepth == 0:
        pass
    elif options.loadFlatDep != '':
        keyTag = "Row Names"
        controlTag = 'control values'
        productionTag = 'production'
        (dependencyMaps,dReport) = loadControlTable(options.loadFlatDep, keyTag, controlTag, productionTag)
    else:
        (dependencyMaps,dReport) = controlDependencies(modelMatrix,resultTag,options,targetPercent,verbose=True,usePersist=usePersist)
        
    if options.outputFileName != '':
        reportFileName = "ControlSearches/Control_dependency_%s_%s_%s_%s_k%s_n%s_s%s.csv" % (resultTag,options.modelName,naturalObjectiveName,syntheticObjectiveName,k,n,str(searchDepth))
    else:
        reportFileName = options.outputFileName
    
    if searchDepth != 0: 
        print "Writing dependency report"
        writer = ReportWriter()
        writer.setFile(reportFileName)
        writer.write(dReport) 
        writer.closeFile()
        
    '''
    Build undirected graph of dependency control sets
    '''
    
    if options.makeGraph:
        dgraph = buildDependencyControlGraph(dependencyMaps)
        graphFileName =  "ControlSearches/Dep_Graph_%s" % persistTag
        if verbose: "writing undirected dependency graph [%s]" % graphFileName
        dgraph.write_ps(graphFileName + ".eps")
        dgraph.write(graphFileName + ".dot")
        
    