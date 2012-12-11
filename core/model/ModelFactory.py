'''
Created on Dec 17, 2010
@author: Graham Rockwell
Factory for generating 
 Linear optimization models
 Control optimization models
'''

from core.do.FluxModel import FluxModel
from core.model.LinearModel import LinearModel
from core.reader.FlatFileParser import FlatFileParser, FluxModelFlatFileParser 
from core.reader.LinearModelParser import MetabolicNetworkParser
from core.util.Config import ReflectionConfig
from core.util.Report import Report
import re, os, pickle

class CreateLinearModel:
    '''
    Object for transforming FluxModel objects into LinearProblem models
    '''
    
    def __init__(self):
        pass
        #self.scale = 1
    
    def _parseReactionBounds( self, fluxModel, analysisName):
        #result = FluxLimit()
        result = {} #! testing removal of Flux Limit    
        
        limits = fluxModel.getFluxLimits(analysisName)
        for name in limits.keys():
            upper = fluxModel.getUpperByName( analysisName, name )
            lower = fluxModel.getLowerByName( analysisName, name )
            
            if upper != None:
                upper = float(upper)
            if lower != None:
                lower = float(lower)
            
            if upper == fluxModel.posInf: 
                xupper = None
            else: xupper = upper
                
            if lower == fluxModel.negInf: 
                xlower = None
            else:xlower = lower
 
            result[name] = (xlower,xupper)
                
        return result

    def _parseNetwork(self,network,limits):
        result = LinearModel()
        result.data.rowCache.extend(network.speciesCache)
        result.data.columnCache.extend(network.reactionCache)
        result.addDataCache(network.dataCache)

        rowNames = network.getOrderedSpeciesNames()
        columnNames = network.getOrderedReactionNames()
        
        for rowName in rowNames:
            result.addRowLimit(rowName,(0.0,0.0))    
                    
        for columnName in columnNames:
            if columnName in limits.keys():
                (lower,upper) = limits[columnName]
                result.addColumnLimit(columnName,(lower,upper))
            
        return result
            
    def fromFluxModel(self,fluxModel,analysisName):
        network = fluxModel.getMetabolicNetwork()
        limits = self._parseReactionBounds(fluxModel,analysisName)
        objective = fluxModel.getObjective(analysisName)
        annotation = fluxModel.annotation
        
        result = self._parseNetwork(network, limits)
        result.setObjective(objective)
        result.annotation = annotation 
        
        return result

class ModelFactory:
    
    def __init__(self):
        #Data holders
        self.configFileName =''
        self.modelConfigFile = ''
        self.protectedTargets = []
        self.subSections = ''
        self.targets = ''
        self.targetSlice = ''
        self.modelAnnotationFile = None #! not used yet
        self.maxControlTagLen = 100
        
        #boolean controls
        self.verbose = False
        self.debug = False
        self.preLoad = False
        self.useGeneMap = False
        self.useReduced = False
        self.geneFill = False
        self.orderGenes = False
        
    def load(self,fileName):
        data = None
        if os.path.exists(fileName):
            if self.verbose: "===loading pre-existing model data=="
            fileHandel = open(fileName,'r')
            data = pickle.load(fileHandel)
            fileHandel.close()
        return data
    
    def parseAnnotationFiles(self,annotationFiles):
        parser = FlatFileParser()
        result = {}
        for (name,keyTag,fileName) in annotationFiles:
            annotation = parser.parseGenericReport(fileName, keyTag, unique = True)
            result[name] = annotation
        return result
    
    def joinAnnotation(self,a1,a2):
        result = {}
        result.update(a1)
        for key in a2.keys():
            if key not in result.keys():
                result[key] = Report()
            result[key].extend(a2[key])
        return result
                    
    def parseModel(self,modelNames,finalName="core"):
        '''
        Generate linear optimization model object
        using names from configuration set in ModelFactory
        @var modelNames: list of the names of the metabolic models
        @type modelNames: [String]
          
        '''
        modelConfig = ReflectionConfig()
        modelConfig.readfp(open(self.modelConfigFile))
        
        fluxModel = None
        #fluxModel = FluxModel()
        pName = ''
        annotation = {}
        for mName in modelNames:
            if mName != '':
                lmParser = MetabolicNetworkParser()
                modelConfig.reflect(mName,lmParser)
                xFluxModel = lmParser.generateModel(mName)
                
                modelParser = FluxModelFlatFileParser()
                modelParser.setWrapper("\'")
                modelConfig.reflect(mName, modelParser)
                jFluxModel = modelParser.generateModel(mName)
                
                if fluxModel == None:
                    fluxModel = jFluxModel
                    pName = fluxModel.objectives.keys()[0]
                else:
                    fluxModel.extend(pName,mName,jFluxModel)
                    pass
                if modelConfig.has_option(mName, "annotationFiles"):
                    annotationFiles = modelConfig.get(mName,"annotationfiles")
                    annotationFiles = eval(annotationFiles)
                    iAnnotation = self.parseAnnotationFiles(annotationFiles)
                    annotation = self.joinAnnotation(annotation, iAnnotation)
                                
        #---------------------------
        # Parse Linear Model
        #---------------------------
        
        modelName = modelNames[0]
        modelMaker = CreateLinearModel()
        model = modelMaker.fromFluxModel(fluxModel, modelName)
        model.annotation = annotation
        model.modelName = ",".join(modelNames)
        
        #---------------------------
        # Select subsections
        #---------------------------
        subsections = self.subSections
        if self.subSections != '':
            subSections = self.subSections.split(',')
            subValues= fluxModel.getMetabolicNetwork().getAnnotation("Section")
            for (name,subValue) in subValues.items():
                if subValue != None and subValue not in subSections:
                    if self.verbose: print "removing %s" % (name)
                    model.removeColumn(name)
                    
        validColumns = model.getColumnNames()
                            
        #-------------------------    
        #Set Filter Target List
        #------------------------
        targets = []
        nonTargets = []
        
        targetTags = fluxModel.getMetabolicNetwork().getAnnotation("ValidTarget")
        for (name,value) in targetTags.items():
            if value == "TRUE":
                targets.append(name)
            if value =="FALSE":
                nonTargets.append(name)
        
        targets = set(validColumns).intersection(targets) 
                
        model.targets = targets
        
        return (fluxModel,model)
    
    def sortTags(self,tags,sep="_"):
        q = []
        result = ""
        ior = False
        iand = False
        
        if len(tags) == 1:
            return tags[0]
            
        while len(tags) != 0:
            t = tags.pop(0)
            if t == "(":
                v = self.sortTags(tags)
                v = "(" + (v)
                q.append(v)
                q.sort()
            elif t == "or":
                if iand:
                    v = (sep+"and"+sep).join(q)
                    result += v
                    q = []
                iand = False
                ior = True        
            elif t == "and":
                if ior:
                    v = (sep+"or"+sep).join(q)
                    result += v 
                    q = []
                ior = False
                iand = True
            elif t == ")":
                if ior:
                    v = (sep+"or"+sep).join(q)
                if iand:
                    v = (sep+"and"+sep).join(q)
                result += sep+v+sep+")"
                return result
            else:
                q.append(t)
                q.sort()
        
        result = q[0]    
        return result
    
    def orderTags(self,tagString,sep="_",replace=" "):
        itagString = tagString.replace(replace,sep)
        itagString = re.sub("%s+"%sep,",",itagString)
        tags = itagString.split(",")
        result = self.sortTags(tags, sep=sep)
        if len(tags) != 1:
            #print "old %s -> new %s" % (tagString,result)
            pass
        return result

    def genValidString(self,v,sep="_",replace=" ",orderTags=False):
        iv = re.sub('[^\w()_]',sep,v)
        if orderTags:
            iv = self.orderTags(iv, sep, replace)
        return iv
    
    def _checkRefactoring(self,ivalues,value,ivalue):
        if ivalue in ivalues.keys():
            if ivalues[ivalue] != value:
                print "String refactoring error old %s != new %s" % (ivalues[ivalue],value)
                ivalue = ivalue + "_1"
        ivalues[ivalue] = value
        return ivalues 
    
    def refactorGeneList(self,data,sep="_"):
        result = []
        for t in data:
            it = self.genValidString(t, sep="_",orderTags=self.orderGenes)
            result.append(it)
        return result
    
    def refactorGeneMap(self,data,sep="_",replace = " "):
        result = {}
        ikeys = {}
        ivalues = {}
        for (key,value) in data.items():
            value = value.replace(replace,sep)
            ikey = self.genValidString(key,sep,False)
            ivalue = self.genValidString(value,sep=sep,replace=replace,orderTags=self.orderGenes)
            if not self.orderGenes:
                ikeys = self._checkRefactoring(ikeys, key, ikey)
                ivalues = self._checkRefactoring(ivalues, value, ivalue)
            
            result[ikey] = ivalue
            
        return result
    
    def refactorReactionGeneSetMap(self,data,sep="_"):
        result = {}
        for (key,values) in data.items():
            ivalues = set()
            for v in values:
                iv = self.orderTags(v, sep)
                ivalues.add(iv)
            result[key] = ivalues
        return result
            
    
    def parseControlMap(self,controlMap,targets,controlFilter = lambda x: True):
        '''
        @type controlMap is variable -< control set
        @type targets: list of targets
        @return: geneMap {var:control}, rGeneMap {var:set(controls)} rGenecurrently only produces gene group -< reaction list
        #! update to remove insure removal of blank gene names, currently sometimes causes problems
        #! consider parsing each gene group into multiple controls
        
        '''
        geneMap = {}
        filterCount = 0
        for (key,value) in controlMap.items():
            if not controlFilter(value):
                filterCount += 1
                #if self.verbose: print "removing control tag (too long %s) [%s]" % (len(value),value)
                continue

            if key in targets:
                if value != '': #! consider removing
                    geneMap[key] = str(value)
        if self.verbose: print "control tags removed [%s]" % (filterCount)
        
        #sort and clean up control tags            
        geneMap = self.refactorGeneMap(geneMap)
        
        rGeneMap = {}
        cGeneMap = {}
        for (key,value) in geneMap.items():
            rGeneMap[key] = set([value])
            if value not in cGeneMap.keys():
                cGeneMap[value] = set()
            cGeneMap[value].add(key) 
        return (geneMap,rGeneMap)
            
            
    def reduce(self,modelMatrix,geneMap,rGeneMap,loadFileName):
        if self.verbose: print "==Performing Reduction=="        
        oMSize = len(modelMatrix.getRowNames())
        oRSize = len(modelMatrix.getColumnNames())
        
        rTargets = modelMatrix.getColumnNames()
        preferedTargets = set(rTargets).difference(modelMatrix.targets)
        
        irGeneMap = rGeneMap
        if self.preLoad and loadFileName != '':
            (reducer,modelMatrix,rGeneMap) = self.load(loadFileName)
        else:
            reducer = ReduceLp()
            reducer.verbose = self.verbose
            reducer.protectedTargets = self.protectedTargets
            reducer.preferedTargets = preferedTargets
            reducer.geneMap = geneMap
            reducer.rGeneMap = rGeneMap    
            
            modelMatrix = reducer.reduceColumns(modelMatrix)  
            modelMatrix = reducer.reduce(modelMatrix,rTargets)
            modelMatrix = reducer.reduceRows(modelMatrix)
        
        geneCluster = None
        if self.useGeneMap:
            (genePairs,rGeneMap,geneCluster) = reducer.reduceGeneMap(geneMap,irGeneMap)
        
        #checkReduction for gene Map
        if rGeneMap != None:
            colNames = modelMatrix.getColumnNames()
            for (rName,gSet) in rGeneMap.items():
                if rName in reducer.removedColumns:
                    del rGeneMap[rName]
                    
            for (rName,gSet) in rGeneMap.items():
                if rName not in colNames and rName :
                    print "--Failed to find control reaction in model!: [%s] : [%s]" % (rName,gSet)
                    del rGeneMap[rName]
                    if rName in reducer.removedColumns:
                        print "reaction member of removed columns"
        
        iMSize = len(modelMatrix.getRowNames())
        iRSize = len(modelMatrix.getColumnNames())
        
        if self.verbose:
            print "Original metabolites (%s) reactions (%s)" % (oMSize,oRSize)
            print "Reduced metabolites (%s) reactions (%s)" % (iMSize,iRSize)
            
        itargets = set(modelMatrix.getColumnNames()).intersection(modelMatrix.targets)
        modelMatrix.targets = itargets
        
        if loadFileName != '':    
            loadFile = open(loadFileName,'w')
            pickle.dump((reducer,modelMatrix,rGeneMap),loadFile)
            loadFile.close()
                
        return (reducer,modelMatrix,rGeneMap,geneCluster)
    
    def geneCluster(self,rGeneMap,geneCluster):
        '''
        find clusters of linked controls
        --Needs to be updated to include all the workings of gene cluster generation.
        '''
        for (rxnName, geneSet) in rGeneMap.items():
            for geneName in geneSet:
                if geneName not in geneCluster.keys():
                    geneCluster[geneName] = set([geneName])
                geneCluster[geneName].add(geneName)
                
        return (geneCluster)    
    
    def sliceTargets(self,modelMatrix):
        targets = []    
    
        if self.targets == '' or self.targetSlice != '':
            targets = modelMatrix.targets
            for target in self.protectedTargets:
                if target in targets:
                    targets.remove(target)
        
        tSlice = self.targetSlice
        if tSlice != '':
            if self.verbose: print "reducing target pool to %s" % (tSlice)
            tSlice = tSlice.split(",")
            a = int(tSlice[0])
            b = int(tSlice[1])
            targets = list(targets)
            targets = targets[a:b]
            if self.verbose: print "Targets: [%s]" % (targets)
        
        if self.targets != '':
            otargets = self.targets.split(",")
            targets = set(targets)
            itargets = targets.union(otargets)
            targets = list(itargets)
            print "Using user specified targets [%s]" % (len(targets))
            if len(self.protectedTargets) != 0:
                protectedTargets = list(self.protectedTargets)
                protectedTargets.extend(targets)
        
        itargets = set(targets)
        if self.useGeneMap:
            #!use only targets for which there are genes
            cTargets = modelMatrix.controlMap.keys()
            itargets = itargets.intersection(cTargets)
        else:
            itargets = itargets.intersection(modelMatrix.getColumnNames())
        targets = itargets
    
        return targets
        
    def loadModel(self,modelNames):
        '''
        Primary function of factory
        '''
        
        #Parse FluxModel object
        (fluxModel,modelMatrix) = self.parseModel(modelNames)
        
        if self.verbose: print "Reactions: %s" % (len(modelMatrix.getColumnNames()))
        if self.verbose: print "Targets %s" % (len(modelMatrix.targets))
            
        #--------------------
        #Get Gene Mappings5
        #--------------------
        geneMap = None
        rGeneMap = None
        if self.useGeneMap:
            controlMap = fluxModel.getMetabolicNetwork().getAnnotation("GeneAssociation")
            (geneMap,rGeneMap) = self.parseControlMap(controlMap,modelMatrix.targets,controlFilter = lambda x : len(x) <= self.maxControlTagLen)
            if self.verbose: print "Control Map %s" % (len(rGeneMap))
        else:
            rGeneMap = {}
            for target in modelMatrix.targets:
                rGeneMap[target] = target
        #----------------------------
        # Reduce Model #! remove reduction functions: complex and don't seem to give much improvement.
        #---------------------------- 
        reducedFileName = "reducedFile_" + str(modelNames)
        reducer = None
        geneCluster = {}
        
        if self.preLoad and os.path.exists(reducedFileName) and False:
            (reducer,modelMatrix,rGeneMap) = self.load(reducedFileName)
        elif self.useReduced:
            if self.verbose: "---Print performing model reductions--"
            originalModel = LinearModel()
            originalModel.extend(modelMatrix)
            (reducer,modelMatrix,rGeneMap,geneCluster) = self.reduce(modelMatrix,geneMap,rGeneMap,reducedFileName)
            #!self.checkReduction(originalModel, modelMatrix) #! this method needs work to produce more informative and readable results
            
        #---------------------------------
        # finish populating geneCluster
        #---------------------------------
        if rGeneMap != None:
            (geneCluster) = self.geneCluster(rGeneMap, geneCluster)
            
        modelMatrix.controlMap = rGeneMap
        targets = self.sliceTargets(modelMatrix)
        modelMatrix.targets = targets
        
        return (fluxModel,modelMatrix,reducer,geneCluster)

    def getReactionGeneMap(self,data):
        result = {}
        
        for key in data.keys():
            value = data[key]
            value = value.replace('(','')
            value = value.replace(')','')
            value = value.replace(' ','')
            value = value.replace('and',' ')
            value = value.replace('or',' ')
            terms =  value.split(' ')
            result[key] = terms
            
        return result
    
    def mapMultipleAnnotations(self,data,catigory=None,sep=' '):
        result = {}

        for key in data.keys():
            value = data[key]
            if catigory != None:
                value = value[catigory]
            terms = re.split(sep,value)
            result[key] = terms
        
        return result
    
    def reverseAnnotation(self,data, list=False):
        result = {}
        for (key,value) in data.items():
            if list:
                for v in value:
                    if v not in result.keys():
                        result[v] = set()
                    result[v].add(key)
            else:
                result[value]=key
        return result
    
    def mergeMAnnotations(self,a1,a2):
        result = {}
        
        for k1 in a1.keys():
            v1 = a1[k1]
            for vi in v1:
                if vi in a2.keys():
                    v2 = a2[vi]
                    if k1 not in result.keys():
                        result[k1] = set()
                    result[k1] = result[k1].union(v2)
                
        return result
