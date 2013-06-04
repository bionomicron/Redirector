'''
Created 2010.5.10
Updated 2013.2.07

@author: Graham Rockwell
Church Lab

User interface command line method for running redirector analysis and framework
4 major parts
1) Parsing model, control map
2) Generating control libraries
3) Generating MILP model 
4) Optimization and Progressive Target Discovery

Example usage:
>python Redirector.py -m SimpleModel1 -b R8 -s R9 --iter 3 --sn 3 -v --sconfig scip_param.set --report  --bt 0.1 
>python Redirector.py -n "Test Simple Model"
>python Redirector.py -m iAF1260 -b Biomass -s EX_C14(e) --iter 1 --sn 1 --report --bt 0.1
>python Redirecror.py -n "Test iAF1260"

'''
from core.model.LinearModelSensitivity import LinearModelVariableBoundarys
#enzymeBoundaryControl
#from core.model.LinearModelTools import ValidateModel
from core.model.ModelFactory import ModelFactory
from core.model.ConstructRegulationOptimization import ConstructRegulationOptimization
from core.model.OptimizationControlRedirector import OptimizationControlRedirector
from core.model.OptimizationControlTools import OptimizationControlFactory
from core.util.Config import ReflectionConfig
from core.util.Report import Report
from core.util.ReportWriter import ReportWriter

from optparse import OptionParser
import os, sys


def main_function():

    parser = OptionParser()
    
    parser.add_option("-v", "--verbose",
                      action="store_true", 
                      dest="verbose", 
                      default=False,
                      help="set verbose mode")
    
    parser.add_option("--masterConfig", 
                      dest="config", 
                      default="redirector.config",
                      help="master configuration file", 
                      metavar="FILE")
                  
    parser.add_option("-c", "--configfile", 
                      dest="configFileName", 
                      default="model_config.txt",
                      help="configuration file", 
                      metavar="FILE")
    
    parser.add_option("-n","--configNames",
                      type = "string",
                      dest = "configNames",
                      default = "Default",
                      help = "comma separated list of configuration settings to include",
                      )    
    
    parser.add_option("-m","--modelname",
                      type = "string", 
                      dest="modelName", 
                      default="",
                      help="name of model from configuration file", 
                      metavar="String")
    
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
                       
    parser.add_option("-o","--outputfile", 
                      dest="outputFileName",
                      default = "rd_output.csv",
                      help="name of report FILE", 
                      metavar="FILE")
    
    parser.add_option("-r", "--result_directory",
                       dest="resultDirectory", 
                       default='../../results/',
                       help = "directory where results are stored")
    
    parser.add_option("--bt", "--biologicalTarget", 
                      dest = "bioTarget", 
                      default = 0.20, 
                      type = float,
                      help = "minimum percentage of the biological objective to maintain", 
                      metavar = "float")
    
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
    
    parser.add_option("--targets",
                      dest = "targets",
                      default = '',
                      help = "Valid reaction targets"
                      )

    parser.add_option("--protectTargets",
                      action = "store_true",
                      dest = "protectTargets",
                      default = False,
                      help = "Do not allow selected targets to be altered"
                      )
    
    parser.add_option("--report",
                      action = "store_true",
                      dest = "isReport",
                      default = False,
                      help = "Write report",
                      metavar = "boolean")
    
    
    parser.add_option("--pl","--preload",
                      dest = "preload",
                      default = 0,
                      type = int,
                      help = "indicate size of pre-made library to load for faster start"
                      )
    
    parser.add_option("--ps", "--preStart",
                      dest = "preStart",
                      default = 0,
                      type = int,
                      help = "start from previous iteration state",
                      metavar = "int"
                      )
    
    parser.add_option("--aneal",
                       dest= "aneal",
                       default= 0.0,
                       help= "use maximum library then reduction")
    
    parser.add_option("--debug",
                      action="store_true",
                      dest = "debug",
                      default = False,
                      help = "turn on debug mode",
                      )
    
    parser.add_option("--gm", "--GeneMap",
                      action="store_true",
                      dest = "useGeneMap",
                      default = False,
                      help = "Use Gene Map",
                      )
    
    parser.add_option("--control",
                      dest="control",
                      default = "flat",
                      help = "control library type, flat,binary,sense")
    
    parser.add_option("--simocontrol",
                      dest="simoControl",
                      default = 1.0,
                      help = "control library type, flat,random,binary,sense")
    
    parser.add_option("--section",
                      dest="subSections",
                      default = '',
                      help = "Comma separated list of sections of the model files to use")
    
    parser.add_option("--targetSlice",
                      dest="targetSlice",
                      default = '',
                      help = "limit maximum target size, mostly for testing controls")
    
    parser.add_option("--primeBounds",
                      dest="primeBounds",             
                      default = False,
                      help = "Do min/max analysis to find boundaries of fluxes")
        
    # Parse options
    (options,args) = parser.parse_args()    
    config = ReflectionConfig()    
    config.readfp(open("Redirector.config"))

    '''
    0. configure analysis
    '''
    
    configNames = options.configNames.split(",")
    configNames.insert(0,"Redirector Model")
    for name in configNames:
        config.merge(name,"Redirector",append=True)
    
    #----------------------------------------
    # reflect options from configuration
    #----------------------------------------

    config.reflect("Redirector",options)
    config.load("Redirector", options.__dict__,override=False)
    
    #-----------------------------------------
    # Check and Build Storage Directories
    #-----------------------------------------
    
    dataDirectory = config.getValue("Redirector","dataDirectory",classType=''.__class__)
    libraryDirectory = config.getValue("Redirector","libraryDirectory",classType=''.__class__)
    resultsDirectory = config.getValue("Redirector","resultDirectory",classType=''.__class__)
    analysisDirectory = config.getValue("Redirector","analysisDirectory",classType=''.__class__)
    
    if not os.path.exists(dataDirectory):
        raise IOError("unable to find required data directory" % dataDirectory)
    if not os.path.exists(resultsDirectory):
        os.makedirs(resultsDirectory)
    if not os.path.exists(libraryDirectory):
        os.makedirs(libraryDirectory)
    if not os.path.exists(analysisDirectory):
        os.makedirs(analysisDirectory)
    
    #----------------------------
    # Parse Inputs
    #----------------------------
    
    verbose        = options.verbose 
    modelName      = options.modelName
    objectiveName  = options.bioObj
    syntheticObjectiveName = options.synthObj
    outputFileName = options.outputFileName
    
    #----------------------------------------------------
    # Initialize and set values for tools and factories
    #----------------------------------------------------
    
    
    naturalObjective = {objectiveName:-1.0}
    syntheticObjective = {syntheticObjectiveName:-1.0}
    protectedTargets = set(naturalObjective.keys()).union(syntheticObjective.keys())
    
    if verbose: print "Redirector Version 1.0"
    if verbose: print "Model names: [%s]" % (modelName)
    if verbose: print "Synthetic objective: [%s]" % (syntheticObjectiveName)
    if verbose: print "Parsing data files for [%s]" % (modelName)
    if verbose: print "Search Size [%s] Iterations [%s]" % (options.searchNeighborhood,options.iterations)
    
    '''
    I. Parse data files and configuration settings
    '''
    
    if verbose: print "----------------Loading Metabolic Models---------------"
    modelNames = modelName.split(",")        
    modelFactory = ModelFactory()
    config.reflect("Redirector",modelFactory)
    modelFactory.protectedTargets = protectedTargets
    (fluxModel,modelMatrix,reducer,geneReduction) = modelFactory.loadModel(modelNames)    
    
    if verbose: print "removing objectives from target set"
    targets = modelMatrix.targets
    targets = set(targets).difference(protectedTargets)
    if verbose: print "Targets List Size [%s]" % len(targets)
    simControl = options.simoControl
    if verbose: print "Simultanious Control %s" % simControl
           
    #-------------------------------------------------------------------
    # Pre discovery of flux bound current not used due to instability
    #-------------------------------------------------------------------
        
    primeFluxBoundaries = {}
    objectiveMinPercent = options.bioTarget
    boundarySearchSize = 1
    boundaryTargets = targets
    boundaryReportFileName = "rd_flux_boundary_M_%s_t_%s_p_%s_s_%s_analysis.csv" % (modelName,len(targets),objectiveMinPercent,boundarySearchSize)
    fluxBoundariesFile = "ControlLibraries/FluxBounds_M_%s_O_%s_T_%s_S_%s" % (modelName,objectiveName,objectiveMinPercent,syntheticObjectiveName)

    print "Prime Bounds [%s]" % options.primeBounds
    naturalFluxBounds = None
    
    if options.primeBounds:
        
        print "finding natural flux bounds"
        naturalFluxBounds = LinearModelVariableBoundarys(modelMatrix, objectiveName=objectiveName, targets=boundaryTargets, pickleFileName = fluxBoundariesFile, minObjectivePercent=objectiveMinPercent,searchSize=boundarySearchSize)
        print "finding production flux bounds"
        #syntheticFluxBounds = LinearModelVariableBoundarys(modelMatrix, objectiveName=syntheticObjectiveName, targets=boundaryTargets, pickleFileName = fluxBoundariesFile, minObjectivePercent=objectiveMinPercent,searchSize=boundarySearchSize)
        if verbose: print "Prime flux boundaries found."
        
        boundaryReport = Report()
        boundaryReport.addColumnHash(objectiveName, naturalFluxBounds)
        #boundaryReport.addColumnHash(syntheticObjectiveName, syntheticFluxBounds)
        
        writer = ReportWriter()       
        writer.setFile(boundaryReportFileName)
        writer.write(boundaryReport) 
        writer.closeFile()
        
    '''
    II. Generate Control Library(ies)
    '''
    
    if verbose:  print "===============Setting up control libraries================="
    processLibrary = OptimizationControlFactory()
    config.reflect("Redirector",processLibrary)
    
    processLibrary.geneMap = modelMatrix.controlMap
    processLibrary.fluxBounds = primeFluxBoundaries
    processLibrary.naturalObjectiveName = objectiveName
    processLibrary.syntheticObjectiveName = syntheticObjectiveName     
    
    controlLibraries = []
    controlMap = {}
    controlTags = options.control.split(",")
    
    if "sense" in options.control:
        print "====>Using sensitivity control library"
        (icontrolLibraries, controlMap) = processLibrary.generateControl(modelMatrix,targets,targetsOnly = False)
        icontrolLibraries = processLibrary._filterControl(icontrolLibraries, targets)
        controlLibraries.extend(icontrolLibraries)
            
    if "flat" in options.control.split(","):
        print "====>Using flat control library"
        (icontrolLibraries) = processLibrary.seedTestLibrary(targets, randomize=False)
        controlLibraries.extend(icontrolLibraries)
        
    if "random" in options.control.split(","): 
        print "====>Using random control library"
        (icontrolLibraries) = processLibrary.seedTestLibrary(targets, randomize=True, factor = 2.0)
        controlLibraries.extend(icontrolLibraries)
        
    if "binary" in options.control.split(","):
        print "====>Using binary control library"
        binaryRange =[-3,-2,-1,0]
        icontrolLibraries = processLibrary.generateBinaryControl(targets,binaryRange)
        controlLibraries.extend(icontrolLibraries)
    
    controlReportName = "ControlLibraries/ControlReport_%s_%s_%s_%s.csv" % (options.configNames,objectiveName,syntheticObjectiveName,len(targets))
    processLibrary.controlReport(modelMatrix, controlLibraries, controlMap, targets, controlReportName)
    
    
    print "============Control Constructed============="
    
    #--------------------------------------------------
    #configure Redirector model construction object 
    #--------------------------------------------------
    if verbose: print "==Formulating Framework=="    
    con = ConstructRegulationOptimization()
    config.reflect("Redirector",con)
    con.controlMax = options.searchNeighborhood
    con.newObjective = syntheticObjectiveName
    
    #----------------------------------------
    # Initialize redirector analysis object
    #----------------------------------------
    if verbose: print "==Running Optimization=="
    redirector = OptimizationControlRedirector()
    config.reflect("Redirector",redirector)
    
    redirector.modelFactory = modelFactory
    redirector.con = con
    redirector.controlFactory = processLibrary
    redirector.controlLibraries = controlLibraries
    redirector.primeBounds = naturalFluxBounds
    #redirector.targets = targets
    redirector.rGeneMap = modelMatrix.controlMap
    #redirector.modelName = modelName
    redirector.naturalObjectiveName = objectiveName
    redirector.syntheticObjectiveName = syntheticObjectiveName
    
    
    #--------------------------------------------------------
    # Debugging option which tests addition / removal of 
    # reactions from the objective
    #--------------------------------------------------------
    
    if "toggle0" in options.control.split(","):
        if verbose: print "-----------control test toggle 0-------------------"
        redirector.testDesignOptimization(modelMatrix)
    if "toggle" in options.control.split(","):
        if verbose: print "-----------control test toggle all-------------------"
        redirector.testDesignOptimization(modelMatrix,fullTest=True)
        
    (ltReport,oPredVal,sPredVal,finalCheckValue,fObjective) = redirector.optimizeControl(modelMatrix)
    
    finalPredictionValue = ltReport["flux"]       
    if verbose: print "Final Production: %s" % (finalPredictionValue[syntheticObjectiveName])

    #-----------------------------
    # Write Report
    #----------------------------
       
    if options.isReport:
        try:
            print "Writing Full Model Report in %s" % (options.resultDirectory)
            redirector.writeReport(ltReport,fluxModel,oPredVal,sPredVal,options.searchNeighborhood,options.iterations,fObjective,options.resultDirectory)
        except:
            print "Unable to write report: %s" % (sys.exc_info[0])
        try:
            targetReportName = "Target_Report_M_%s_N_%s_S_%s_K_%s_I_%s.txt" % (modelName,objectiveName,syntheticObjectiveName,options.searchNeighborhood,options.iterations)
            print "Writing Optimization target Report %s" % (targetReportName)
            iteration = options.preStart + options.iterations
            redirector.writeAnnotatedTargets(fObjective, modelMatrix, annotationName = "bnumber", regex="[a-zA-Z0-9\(\)]+", iteration = iteration, oPrediction = oPredVal, nPrediction = finalPredictionValue)
        except:
            print "failed to write %s" % (targetReportName)
    print "Done"
        
if __name__  == "__main__":
    main_function()
