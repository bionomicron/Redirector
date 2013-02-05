'''
Created on Dec 12, 2012
Updated Feb 05 2013
@author: Graham Rockwell
Church Lab
'''

from core.model.ModelFactory import ModelFactory
from core.model.LPSolver import LPSolver
from core.util.Config import ReflectionConfig
from core.util.ReportWriter import ReportWriter
from optparse import OptionParser
from time import strftime
import os

def main_function():

    parser = OptionParser()
    
    parser.add_option("-v", "--verbose",
                      action="store_true", 
                      dest="verbose", 
                      default=False,
                      help="set verbose mode")
    
    parser.add_option("--config", 
                      type="string",
                      metavar="FILE",
                      dest="config", 
                      default="redirector.config",
                      help="Master configuration file that provides most setting for the Flux Balanace Analysis. Extra variable settings will be ingnored.  Hence the same config can be used for multiple analysis functions", 
                      )
                  
    parser.add_option("-c", "--modelConfig",
                      type="string",
                      metavar="FILE",
                      dest="configFileName",
                      default="model_config.txt",
                      help="Configuration file which sets model files" 
                      )
    
    parser.add_option("-n","--configNames",
                      type = "string",
                      metavar = "String",
                      dest = "configNames",
                      default = "Default",
                      help = "A comma separated list of the names of configurations to use, as set out in the configuration file")    
    
    parser.add_option("-m","--modelname",
                      type = "string",
                      metavar = "String",
                      dest="modelName", 
                      default="",
                      help="Name of model(s) from the modelconfiguration file" 
                      )
    
    parser.add_option("-b","--bioObjective", 
                      type = "string",
                      dest="bioObj", 
                      default = 'Biomass',
                      help="Name / ID of biological objective reaction", 
                      metavar="String")
    
    parser.add_option("-o","--outputfile", 
                      dest="outputFileName",
                      default = None,
                      help="Name of report file to be generated", 
                      metavar="FILE")
    
    parser.add_option("-r", "--result_directory",
                      type="string",
                      metavar="directory",
                      dest="resultDirectory", 
                      default='../../results/',
                      help = "Directory where results are stored")
        
    parser.add_option("--targets",
                      type="string",
                      metavar="String",
                      dest = "targets",
                      default = '',
                      help = "List of valid reaction target, if left blank will be automatically filled in from metabolic network file")

    parser.add_option("--report",
                      action = "store_true",
                      dest = "isReport",
                      default = False,
                      help = "When this tag is used a report will be generated when the analysis is finished",
                      metavar = "boolean")
    
    
    parser.add_option("--debug",
                      action="store_true",
                      dest = "debug",
                      default = False,
                      help = "turn on debug mode")
    
    parser.add_option("--gm", "--GeneMap",
                      action="store_true",
                      dest = "useGeneMap",
                      default = False,
                      help = "Use Gene Map")
    
    parser.add_option("--section",
                      type="string",
                      metavar="String",
                      dest="subSections",
                      default = '',
                      help = "Comma separated list of sections of the model files to use")
        
    #-------------------------------
    # Parse options
    #-------------------------------
    
    (options,args) = parser.parse_args()    
    config = ReflectionConfig()    
    config.readfp(open("Redirector.config"))

    #---------------------------
    # configure preset analysis
    #---------------------------
    
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
    resultsDirectory = config.getValue("Redirector","resultDirectory",classType=''.__class__)
    analysisDirectory = config.getValue("Redirector","analysisDirectory",classType=''.__class__)
    
    if not os.path.exists(dataDirectory):
        raise IOError("unable to find required data directory" % dataDirectory)
    if not os.path.exists(resultsDirectory):
        os.makedirs(resultsDirectory)
    if not os.path.exists(analysisDirectory):
        os.makedirs(analysisDirectory)
    
    #----------------------------
    # Parse Inputs
    #----------------------------
    
    verbose        = options.verbose 
    modelName      = options.modelName
    objectiveName  = options.bioObj
    outputFileName = options.outputFileName
    
    '''
    #----------------------------------------------------
    # Initialize and set values for tools and factories
    #----------------------------------------------------
    '''
    
    naturalObjective = {objectiveName:-1.0}
    
    if verbose: print "Flux Balanace Analysis Version 1.6"
    if verbose: print "Model names: [%s]" % (modelName)
    if verbose: print "Parsing data files for [%s]" % (modelName)
    
    '''
    I. Parse data files and configuration settings
    '''
    
    if verbose: print "----------------Loading Metabolic Models---------------"
    modelNames = modelName.split(",")        
    modelFactory = ModelFactory()
    config.reflect("Redirector",modelFactory)
    (fluxModel,modelMatrix,reducer,geneReduction) = modelFactory.loadModel(modelNames)    
    
    model=modelMatrix
    if verbose: print "removing objectives from target set"
    targets = modelMatrix.targets
    if verbose: print "Targets List Size [%s]" % len(targets)
        
    lps = LPSolver()
    predictions = lps.run(model,naturalObjective)
    objectiveValue = predictions[objectiveName]
    lps.clear()
    
    if verbose: print "Optimized Objective [%s] Value [%s]" % (objectiveName,objectiveValue)

    report = fluxModel.getReport()
    report.addColumnHash("Flux", predictions)
    
    if outputFileName == None or outputFileName == '':
        outputFileName = resultsDirectory + "FBA_%s_%s.txt" % (modelName,strftime('%Y%m%dT%H%M%S'))
        
    writer = ReportWriter()
    writer.setFile(outputFileName)
    writer.write(report) 
    writer.closeFile()
    
    if verbose: print "Report Written [%s]" % (outputFileName)
    
    
if __name__  == "__main__":
    main_function()