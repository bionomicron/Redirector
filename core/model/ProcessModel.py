#Graham Rockwell
#7.6.2009

from core.model.LinearModel import LinearModel
from core.model.ModelFactory import CreateLinearModel
from core.reader.FlatFileParser import FluxModelFlatFileParser
from core.util.Report import Report

import re, os, pickle

class MatrixTools:
    
    def __init__(self):
        pass
    
    def vectorMultiply(self,v1,v2,verbose = False):
        result = 0;
        s = ''
        for key in v1.keys():
            if key in v2.keys():
                value1 = v1[key]
                value2 = v2[key]
                result += value1*value2
                s += " + %s[%s]" % (key,value1*value2)
                #print "%s = %s * %s" % (result,value1,value2)
        if verbose:
            print "%s = %s" % (result,s[3:])
        return result
    
    def vectorMScalar(self,v1,s):
        result = {}
        for key in v1.keys():
            value = v1[key]
            result[key] = value * s
        return result
    
    def matrixVectorM(self,m,v):
        result = {}
        delta = 1e-3
        for rName in m.getRowNames():
            v1 = m.getRow(rName)
            if v1 == None:
                pass
            if v1 != None:
                iv = vectorMultiply(v,v1)
                result[rName] = iv
        return result

    
class ProcessGeneticObjective:
    
    def init(self):
        pass
    
    def translateGeneNumber(self,term,annotation):
            
        group = term.replace("_"," ")
        geneGroup = term.replace("_","").replace("and","_").replace("or","_")
        geneGroup = geneGroup.replace("(","").replace(")","").split("_")
    
        for geneNumber in geneGroup:
            if geneNumber in annotation.keys():
                iName = annotation[geneNumber]["gene"]
                group = group.replace(geneNumber,iName)
                
        return group
    
    def processObjective(self,geneticObjective,annotation,nb,iter,report=None,clusterData=None):
        verbose = True
        if report == None:
            report = Report()

        names = {}
        function = {}
        adjustment = {}
        tag = {}
        
        clusterTags = {}
        for key in geneticObjective.keys():
            value = geneticObjective[key]
            
            checkName = ''
            if key == checkName:
                pass
            
            group = key.replace("_"," ")
            geneGroup = key.replace("_","").replace("and","_").replace("or","_")
            geneGroup = geneGroup.replace("(","").replace(")","").split("_")
   
            for geneNumber in geneGroup:
                adjustment[geneNumber] =value 
                if geneNumber in annotation.keys():
                    iName = annotation[geneNumber]["gene"]
                    names[geneNumber] = iName
                    function[geneNumber] = annotation[geneNumber]["function"]
                    group = group.replace(geneNumber,iName)
        
            for geneNumber in geneGroup:
                tag[geneNumber] = group
            
            if clusterData != None:
                if key not in clusterData.keys():
                    continue
                
                cluster = clusterData[key]
                clusterTag = str(cluster).replace("set([","").replace("])","").replace(", ",",").replace(" ","").replace("'","")
                
                clusterTag = clusterTag.replace("_"," ")

                clusterGroup = clusterTag.replace("_","").replace("and","_").replace("or","_")
                clusterGroup = clusterGroup.replace("(","").replace(")","").replace(" ","").replace(",","_")
                clusterGroup = clusterGroup.split("_")
                
                for geneNumber in clusterGroup:
                    if geneNumber in annotation.keys():
                        iName = annotation[geneNumber]["gene"]
                        clusterTag = clusterTag.replace(geneNumber,iName)
                
                verbose = False
                for geneNumber in geneGroup:
                    clusterTags[geneNumber] = clusterTag
                    
                    if verbose: 
                        print "--- cluster---"
                        print geneNumber
                        print names[geneNumber]
                        print tag[geneNumber]
                        print str(cluster)
                        print clusterGroup
                        print clusterTag


        report["gene name"] = names
        report["gene ID"] = tag
        if clusterData != None:
            report["gene cluster"] = clusterTags        
        report["function"] = function
        adjName = "%s_%s" % (nb,iter)
        report[adjName] = adjustment
        return report
    
    def stringToMap(self,value):
        result = eval(value)
        return result
    
    def processAnalysis(self,analysisData,annotationData,nbKey,iterKey,objKey,clusterData=None):
        if analysisData == None:
            return None
        
        report = Report()
        targets = {}
        maxIter = 0
        maxNb = 0
        for key in analysisData.keys():
            data = analysisData[key]
            nb = int(data[nbKey])
            iter = int(data[iterKey])
            targetString = data[objKey].rstrip()
            #print targetString
            target = eval(targetString)
            targets[(nb,iter)] = target
        
        if clusterData != None:
            clusterTargets = {}
            for key in targets.keys():
                target = targets[key]
                newTarget = {}
                newTarget.update(target)
                for tk in target.keys():
                    if tk in clusterData.keys():
                        for clusterTarget in clusterData[tk]:
                            iValue = target[tk]
                            newTarget[clusterTarget] = iValue
                clusterTargets[key] = newTarget
            targets =clusterTargets
            
        if annotationData != None:    
            k = targets.keys()
            k.sort()

            for (nb,iter) in k:
                target = targets[(nb,iter)]
                report = self.processObjective(target,annotationData,nb,iter,report,clusterData)

        return (targets,report)
            
        

class ProcessModel:
    
    def __init__(self):
        self.configFileName =''    
 
        self.verbose = False
        self.preLoad = False
        self.useGeneMap = False
        self.useReduced = False
        self.model = None
        
        self.geneFill = False
        
    def load(self,fileName):
        data = None
        if os.path.exists(fileName):
            if self.verbose: "===loading pre-existing model data=="
            fileHandel = open(fileName,'r')
            data = pickle.load(fileHandel)
            fileHandel.close()
        return data
    
    def getFluxModel(self,configFileName,modelNames):
        gconfigParser = GeneralMultiConfigParser()
        gconfigs = gconfigParser.parse(configFileName)
        modelParser = FluxModelFlatFileParser()
        modelParser.setWrapper("\'")
        
        fluxModel = None
        pName = ''
        iFluxModel = LinearModel()
        for mName in modelNames:
            if mName != '':
                gconfig = gconfigs[mName]
                iFluxModel = modelParser.parseGenConfig(gconfig, mName)
                if fluxModel == None:
                    fluxModel = iFluxModel
                    pName = fluxModel.objectives.keys()[0]
                else:
                    fluxModel.extend(pName,mName,iFluxModel)
        
        #---------------------------
        # Parse Linear Model
        #---------------------------
        
        modelName = modelNames[0]
        modelMaker = CreateLinearModel()
        modelMaker.setScale(1)
        modelMatrix = modelMaker.fromFluxModel(fluxModel, modelName)
        
        #-------------------------    
        #Set Filter Target List
        #------------------------
        nonTargets = []
        targetList = fluxModel.getMetabolicNetwork().getAnnotation("ValidTarget")
        for name in targetList.keys():
            value = targetList[name]
            if value =="FALSE":
                nonTargets.append(name)
        
        return (fluxModel,modelMatrix,nonTargets)
    
    def refactorStringMap(self,data):
        result = {}
        for (key,value) in data.items():
            ikey = key.replace(" ","_")
            ivalue = value.replace(" ","_")
            result[ikey] = ivalue
        return result
    
    def parseGeneMap(self,fluxModel,nonTargets):
        '''
        gene map is rxn key to gene name value
        currently only produces gene group -< reaction list
        #! update to remove insure removal of blank gene names
        '''
        reactionGeneMap = fluxModel.getMetabolicNetwork().getAnnotation("GeneAssociation")
        geneMap = {}
        for (key,value) in reactionGeneMap.items():
            #if key not in nonTargets and value != '':
            if key not in nonTargets:
                geneMap[key] = str(value)
            elif key not in nonTargets:
                if self.geneFill:
                    geneMap[key] = str("gene_"+key)
                    if self.verbose: print "filling in gene for reaction %s" % (key)
            else:
                if self.verbose: print "removing reaction gene pair [%s:%s]" % (key,value)
                    
        geneMap = self.refactorStringMap(geneMap)
        rGeneMap = {}
        for (key,value) in geneMap.items():
            #vs = set([value])
            rGeneMap[key] = set([value])
        
        return (geneMap,rGeneMap)
            
    def parseGenePairs(self,reducer,geneMap,rGeneMap):
        #-----------------------------
        # Adjust for model reduction
        # -gene pairs lists relation between reactions and genes  many genes to many reduced reactions
        # -point gene conglomerate to geneNames of memebers
        #-----------------------------
        genePairs = []
        geneCluster = {}
        
        if self.useReduced:
            (genePairs,rGeneMap,geneCluster) = reducer.reduceGeneMap(geneMap,rGeneMap)
        
        #-------------------------------------
        # Transfer remainder of gene mapping 
        # to gene pair list
        #-------------------------------------
        for (rxnName,geneNames) in rGeneMap.items():
            for geneName in geneNames:
                genePairs.append((rxnName,geneName))
                if geneName not in geneCluster.keys():
                    geneCluster[geneName] = set([geneName])
                geneCluster[geneName].add(geneName)
                
        return (genePairs,geneCluster)
    
    
    def reduce(self,modelMatrix,geneMap,rGeneMap,nonTargets,objective,syntheticObjectiveName,loadFileName):
        if self.verbose: print "==Preforming Reduction=="        
        
        reducer = ReduceLp()
        reducer.preferedTargets = nonTargets
        reducer.geneMap = geneMap
        reducer.rGeneMap = rGeneMap
        rTargets = modelMatrix.getColumnNames()
        rTargets.remove(objective)
        rTargets.remove(syntheticObjectiveName)
        modelMatrix = reducer.reduceColumns(modelMatrix)  
        modelMatrix = reducer.reduce(modelMatrix,rTargets)
        modelMatrix = reducer.reduceRows(modelMatrix)
            
        loadFile = open(loadFileName,'w')
        pickle.dump((reducer,modelMatrix),loadFile)
        loadFile.close()
        return (reducer,modelMatrix)
        
        
    def loadModel(self,modelNames,objective,syntheticObjectiveName):
        #------------------
        # Parse Flux Model
        # Set model matrix
        #------------------
        (fluxModel,modelMatrix,nonTargets) = self.getFluxModel(self.configFileName,modelNames)
        
        #update to allow for more complex objective
        modelName = str(modelNames)
        modelMatrix.setObjective({objective:-1})
        
        nonTargets.append(objective)
        nonTargets.append(syntheticObjectiveName)
                
        if self.verbose: print "Non Targets %s" % (len(nonTargets))
            
        #--------------------
        #Get Gene Mappings
        #--------------------
        geneMap = None
        rGeneMap = None
        if self.useGeneMap:
            (geneMap,rGeneMap) = self.parseGeneMap(fluxModel,nonTargets)


        #----------------------------
        # Reduce Model
        #---------------------------- 
        reducedFileName = "reducedFile_" + modelName + "_" + objective + "_" + syntheticObjectiveName
        reducer = None
        if self.preLoad and os.path.exists(reducedFileName):
            (reducer,modelMatrix) = self.load(reducedFileName)        
        elif self.useReduced:
            (reducer,modelMatrix) = self.reduce(modelMatrix,geneMap,rGeneMap,nonTargets,objective,syntheticObjectiveName,reducedFileName)
            
        #---------------------------------
        # Populate gene reation pair list
        #---------------------------------
        genePairs = None
        geneCluster = None
        if geneMap != None:
            (genePairs,geneCluster) = self.parseGenePairs(reducer,geneMap,rGeneMap)
            
        return (fluxModel,modelMatrix,reducer,nonTargets,genePairs,geneMap,rGeneMap,geneCluster)

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
            
    
class OptimizationPriceAnalysis:
    
    def init(self):
        pass
    
    def findAllPrices(self,model,targets,dim=2):
        result = {}
        lo = LinearOptimization()
        lo.setModel(model)
        lo.runSimplex()
        for i in targets:
            prices = lo.getPointReducedCosts({i: -1.0},dim)
            for j in prices.keys():
                term = (i,j)
                result[term] = prices[j]
        lo.clear()
        return result
    
    def findTargetPrices(self,model,objectives,targets,dim=2,coeff = -1.0):
        result = Report()
        lo = LinearOptimization()
        lo.setModel(model)
        
        for obj in objectives:
            lo.clearObjective()
            lo.setObjectiveMap({obj:-1})
            lo.runSimplex()        

            for i in targets:
                name = obj + "_" + i
                result[name] =  lo.getPointReducedCosts({i: coeff},dim)
                
        lo.clear()
        
        return result
    
    def findReducedPrices(self,model,objective,target,dim=2):
        lo = LinearOptimization()
        lo.setModel(newModel)
        lo.clearObjective()
        lo.setObjective(objective)
        lo.runSimplex()
        prices = lo.getPointReducedCosts(targets,dim)
        return prices

class ModelProperties:
    
    def __init__(self):
        pass
    
    def findMinMax(self, model, objectiveName, minValue, targets):
        positive = {}
        negative = {}
        imodel = LinearModel()
        imodel.extend(model)
        if minValue < 0:
            imodel.addColumnLimit(objectiveName,(None,minValue))
        else:
            imodel.addColumnLimit(objectiveName,(minValue,None))
        
        lp = LinearOptimization()
        lp.setModel(model)
        for t in targets:
            lp.clearObjective()
            lp.setObjectiveMap({t: -1})
            lp.runSimplex()
            fluxes = lp.getPredictionMap()
            lowValue = fluxes[t]
            negative[t] = lowValue
            
            lp.clearObjective()
            lp.setObjectiveMap({t: 1})
            lp.runSimplex()
            fluxes = lp.getPredictionMap()
            highValue = fluxes[t]
            positive[t] = highValue
            
        return (positive,negative)
    
        