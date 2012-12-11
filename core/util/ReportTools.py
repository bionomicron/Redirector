from core.reader.FlatFileParser import FlatFileParser
from core.util.Report import Report
from core.util.ReportWriter import ReportWriter
from core.model.ModelFactory import ModelFactory
from time import re, strftime



class ReportTools:
    
    def __init__(self,model = None):
        self.model = model
        self.annotations = None
        
        #Boolean controls
        self.orderGenes = False
    
    def generateLimitReport(self,limits):
        result = Report()
        upper = {}
        lower = {}
        for key in limits.keys():
            (lower[key],upper[key]) = limits[key]
        result.addReportColumnHash("LowerBound",lower)
        result.addReportColumnHash("UpperBound",upper)
        
        return result
    
    def generateLongNameReport(self,fluxModel):
        report = Report()
        for name in fluxModel.getMetabolicNetwork().getOrderedReactionNames():
            reaction = fluxModel.getMetabolicNetwork().getReaction(name)
            longName = reaction.getName()
            equation = reaction.getEquation()
            report.addElement(name,"LongName",longName)
            report.addElement(name,"Equation",equation)
        return report
    
    def generateReport(self,fluxModel, modelName, predictions, limits = False, longName = False,):
        result = Report()
        if longName:
            longNamesReport = self.generateLongNameReport(fluxModel)
            result.extend(longNamesReport)
        result.addReportColumnHash("FluxPrediction",predictions)
        if limits:
            limits = fluxModel.getLimits(modelName)
            limitReport = self.generateLimitReport(limits)
            result.extend(limitReport)

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
    
    def readAnnotatedTargets(self,targetFile,model,annotationName,regex,bNames = False,orderGenes = False):
        parser = FlatFileParser()
        modelFactory = ModelFactory()
        modelFactory.orderGenes = orderGenes
        report = parser.parseToReport(targetFile, "Row Names", header=["Control","Flux","Reaction Control","Flux Changes"])
        controlValues = report.getColumn("Control")
        if bNames:
            controlValues = self._deannotateGenes(controlValues, model, annotationName, regex)
        if self.orderGenes:
            controlValues = modelFactory.refactorStringMap(controlValues, sep="_")
        controlScore = report.getColumn("Flux")
        if bNames:
            controlScore =  self._deannotateGenes(controlScore, model, annotationName, regex)
        
        return(controlValues,controlScore)    
                
    def writeAnnotatedTargets(self,objective,model,annotationName,regex,iteration,oPrediction={},nPrediction={},outputFileName=None):
        
        if outputFileName == None:
            outputFileName = "Target_Report_T_%s_C_%s_M_%s_N_%s_S_%s_K_%s_I_%s.txt" % (self.resultTag,self.control,self.modelName,self.naturalObjectiveName,self.syntheticObjectiveName,self.con.controlMax,iteration)
            
        aValue = self._annotateObjective(objective, model, annotationName, regex)
        controlValues = self._condenseObjective(objective)
        controlValues = self._annotateGenes(controlValues, model, annotationName, regex)
        controlScores = self._controlScore(objective,oPrediction, nPrediction)
        controlScores = self._annotateGenes(controlScores, model, annotationName, regex)
        controlScore = self._condenseMap(controlScores,sum)
        
        report = Report()
        report.addColumnHash("Control", aValue)
        report.addColumnHash("Flux", controlScore)
        report.addColumnHash("Reaction Control", controlValues)
        report.addColumnHash("Flux Changes", controlScores)

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
    
    def writeReport(self,ltReport,fluxModel, oPredVal,sPredVal,s2PredVal,neighborhood,iterations,fObjective,resultDir):
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
        report.addColumnHash("Synthetic mBio", s2PredVal)
        report.extend(ltReport)
            
        report.addColumnHash("Final Objective", fObjective)
    
        #--------------------
        #Write Report
        #--------------------
        if self.verbose: print "Writing report ..."
        outputFileName = resultDir + "RD_" + self.resultTag + "_" + self.modelName + "_" + self.syntheticObjectiveName + "_" + str(neighborhood) + "_" + str(iterations) +  "_" + strftime('%Y%m%dT%H%M%S') + ".txt"  
        
        writer = ReportWriter()
        writer.setFile(outputFileName)
        writer.write(report) 
        writer.closeFile()
        if self.verbose: print "Report Written [%s]" % (outputFileName)
        
        return report 
        