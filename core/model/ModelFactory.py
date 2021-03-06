'''
Created on Dec 17, 2010
Updated on Feb 18, 2013
@author: Graham Rockwell
Factory for generating 
 Linear optimization models
 Control optimization models
'''

from core.model.LinearModel import LinearModel
from util.FlatFileParser import FlatFileParser
from core.reader.FluxModelParser import FluxModelFlatFileParser 
from core.reader.LinearModelParser import MetabolicNetworkParser
from util.Config import ReflectionConfig
from util.Report import Report
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
            annotation = parser.parseGenericReport(fileName, keyTag, unique=True)
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
                    
    def parseModel(self,modelNames):
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
                #xFluxModel = lmParser.generateModel(mName)
                
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
    
    def numericalRefactorModel(self,model,floor=1.0,minValue=1e-3,maxFactor=100,verbose=True):
        '''
        Method for numerically refactoring a linear model
        Currently only has floor option.
        '''
        columnNames = model.getColumnNames()
        
        for cName in columnNames:
            if cName in model.getObjective().keys():
                continue
            
            rMap = model.getColumnValueMap(cName)
            rMin = abs(min(rMap.values(), key = lambda x:abs(x)))
            factor = floor/rMin
            
            if factor > maxFactor: factor = maxFactor
             
            if rMin < floor:
                #if verbose: print "original reaction: [%s]" % rMap    
                if verbose: print "refactoring %s [%s]" % (cName,factor)
                
                for rName in rMap.keys():
                    value = rMap[rName]
                    if abs(value) < minValue:
                        model.removeData(rName,cName)
                    else:
                        iValue = value * factor
                        model.addData(rName,cName,iValue)
                rMap = model.getColumnValueMap(cName)
                #if verbose: print "final reaction [%s]" % rMap
                    
                #solver = LPSolver()
                #oPredVal = solver.run(model=model)
                #print "check production [%s]" % oPredVal["Biomass"]
                #print "done"
        
        return model
            
    
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
        
        #Parse flux model files
        (fluxModel,modelMatrix) = self.parseModel(modelNames)
        
        if self.verbose: print "Reactions: %s" % (len(modelMatrix.getColumnNames()))
        if self.verbose: print "Targets %s" % (len(modelMatrix.targets))
        
        #---------------------------------------------------------------------------
        # Adjust or refactor model to reduce presence of very small coefficents
        #---------------------------------------------------------------------------
        coefficentFloor = 1.0
        if self.verbose: print "Refactoring coefficients of linear model to [%s]" % (coefficentFloor)
        modelMatrix = self.numericalRefactorModel(modelMatrix, floor=coefficentFloor, verbose=True)
            
        #--------------------
        #Get Gene Mappings
        #--------------------
        geneMap = None
        rGeneMap = None
        controlMap = None
        if self.useGeneMap:
            controlMap = fluxModel.getMetabolicNetwork().getAnnotation("GeneAssociation")
            (geneMap,rGeneMap) = self.parseControlMap(controlMap,modelMatrix.targets,controlFilter = lambda x : len(x) <= self.maxControlTagLen)
            if self.verbose: print "Control Map %s" % (len(rGeneMap))
        else:
            rGeneMap = {}
            for target in modelMatrix.targets:
                rGeneMap[target] = target

        #---------------------------------
        # finish populating geneCluster
        #---------------------------------
        geneCluster = {}
        if controlMap != None:
            (geneCluster) = self.geneCluster(controlMap, geneCluster)
            
        modelMatrix.controlMap = rGeneMap #important to use this after filtering and fixing gene set.
        modelMatrix.controlCluster = geneCluster
        targets = self.sliceTargets(modelMatrix)
        modelMatrix.targets = targets
        
        return (fluxModel,modelMatrix)

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