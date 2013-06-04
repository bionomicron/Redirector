 #!/usr/bin/env python2.7
'''
@author: Graham Rockwell
Church Group Havard Genetics
Last updated 8/28/2009
@summary:
    User interface methods for optimization of recombination oligos
    Method "auto": determines which methods to run given from given input parameters
    Method "location": selected locations from genome sequence with replacement 
    Method "oligo": uses input fasta file of oligos to find optimal location and generate primers
    Method "sequence": take recombination oligo fasta file compares it to sequencing fasta file
'''

from core.reader.FlatFileParser import FlatFileParser

from core.util.Config import ReflectionConfig
from core.util.Report import Report
from core.util.ReportWriter import ReportWriter

from core.genetic.SequenceTools import SequenceTools, SequenceFactory, RecombinationOligoFactory,ControlRegionTools

from Bio import Entrez,SeqIO
from Bio.SeqRecord import SeqRecord

from optparse import OptionParser
import pickle, time


if __name__  == "__main__":
    parser = OptionParser()
    
    parser.add_option("-v", "--verbose",
                      action="store_true", 
                      dest="verbose", 
                      default=False,
                      help="print status messages to stdout")
    
    parser.add_option("-c", "--configFile",
                      dest="config", 
                      default="oligo_generation.config",
                      help="config file for oligo generation, check config for defaults",
                      metavar="FILE")
    
    parser.add_option("-n", "--configName", 
                      dest="configName", 
                      default="default",
                      help="which configuration setting to use", 
                      metavar="string")
    
    parser.add_option("-m", "--mode",
                      dest="mode", 
                      default="auto",
                      help="oligo design mode: auto, oligo, sequence", 
                      metavar="string")

    parser.add_option("-d", "--databaseFile", 
                      dest="databaseFile", 
                      default="",
                      help="file for blast database to use", 
                      metavar="FILE")
    
    parser.add_option("-b", "--blast_exe", 
                      dest="blastExe", 
                      default="/usr/bin/blastall",
                      help="location of blastall executable", 
                      metavar="FILE")

    parser.add_option("-l", "--locationFile", 
                      dest="locationFile", 
                      default="",
                      help="Tab delmited file of location information from which to create recombination oligos", 
                      metavar="FILE")

    parser.add_option("-t", "--targetFile", 
                      dest="targetFile", 
                      default="",
                      help="starting point for targets, input file of oligos to optimize, in fasta format", 
                      metavar="FILE")
    
    parser.add_option("-o","--oligoFile", 
                      dest="oligoFile",
                      help="name of report FILE", 
                      default = "recombination_oligo_file.txt",
                      metavar="FILE")
    
    parser.add_option("-s", "--sequencing",
                      dest="sequencingFile", 
                      default="",
                      help="Name of FASTA file of sequencing results to be compared to recombination oligos",
                      metavar="FILE")
    
    parser.add_option("--recCompare",
                      dest="seqCompareLog", 
                      default="recombination_comparison_log.txt",
                      help="Name of file to write sequencing of recombination analysis to",
                      metavar="FILE")
    
    parser.add_option("--flankSize", 
                      dest="flankingSize", 
                      default=45.0,
                      help="size of flanking homologus regions", 
                      metavar="FILE")
    
    parser.add_option("--primerTM", 
                      dest="primerTM", 
                      default=55.0,
                      help="melting temperature of suggested primers", 
                      metavar="FILE")
    
    parser.add_option("--minDeltaG", 
                      dest="minDeltaG", 
                      default=-13.0,
                      help="minimum folding energy of oligos", 
                      metavar="FILE")
                  
    #parse oligos
    (options,args) = parser.parse_args()
    configFileName = options.config
    configName = options.configName
    #configName = "fatty acid targets"
    
    # Parse main config and reflect into options
    config = ReflectionConfig()    
    config.readfp(open(configFileName))
    config.reflect(configName,options)
    
    #----------------------------
    #load command line
    #----------------------------
    
    verbose = options.verbose
    mode = options.mode
    
    locationFile = options.locationFile
    targetFile = options.targetFile
    oligoFileName = options.oligoFile
    
    blastDB = options.databaseFile
    blastExe = options.blastExe
    recSeqFileName = options.sequencingFile
    recCompareFileName = options.seqCompareLog
    
    
    #------------------------
    #Static Values
    #------------------------
    
    targetTag = "([^ACGTacgt]+.*[^ACGTacgt]+)"
    boundary = options.flankingSize
    #boundary = 45
    searchSize = 50
    cutOff = (None,-13.0)
    primerSize = 20
    primerTm = 58
    primerDistance = 250
    
    #--------------------------
    #Start factories
    #--------------------------
    
    fileParser = FlatFileParser()
    seqTools = SequenceTools()
    seqTools.verbose = verbose
    seqFac = SequenceFactory()
    seqFac.verbose = verbose
    recFac = RecombinationOligoFactory()
    recFac.verbose = verbose
    
    oligoRecords = None
    gRecord = None

    #----------------------------------
    # Generate new oligos from targets
    #----------------------------------
    
    if mode == "location" or (locationFile != '' and mode == "auto"):
        controlRegionTools = ControlRegionTools()
        
        targetReport = fileParser.parseToReport(locationFile,"TargetName",header=["TargetName","start","end","direction","replacement",])
        conversionMap = {"geneID":"gene","f":int(1),"r":int(-1)}
        targetFeatures = seqFac.convertReportToFeatures(targetReport,conversionMap,"start","end","direction")
    
        genomicSeqFileName = seqFac.downloadGenBank()
        time.sleep(5)
        gRecord = seqFac.parseGenBank(genomicSeqFileName)
        seqData = str(gRecord.seq).lower()
        
        oligoRecords = list()
        for targetName in targetReport.returnRowNames():
            location = float(targetReport.get(targetName, "start"))
            end = float(targetReport.get(targetName, "end"))
            dir = int(targetReport.get(targetName, "direction"))
            replace = targetReport.get(targetName, "replacement")
            seq = controlRegionTools.getSequenceRegion(seqData,location,size=1,boundary=boundary,dir=dir,useUpper=True,replace=replace,end=end)
            
            iRecord = SeqRecord(seq)
            iRecord.id = targetName
            iRecord.name = targetName
            oligoRecords.append(iRecord)
           
        print "sequences found" 
        
        oligoTargetFileHandle = open(targetFile,"w")
        SeqIO.write(oligoRecords,oligoTargetFileHandle,"fasta")
        oligoTargetFileHandle.close()
        
        print "target sequence file written"
    
    if mode == "sequence" or (mode == "auto" and targetFile != ''):
        
        # Parse Target Oliog if they allready exist.
        if oligoRecords == None or True:
            oligoRecords = SeqIO.parse(open(targetFile), "fasta")
            oligoRecords = list(oligoRecords)
        
        #--------------------------------
        #Get genomic data for primers
        #--------------------------------
        if gRecord == None:
            genomicSeqFileName = seqFac.downloadGenBank()
            time.sleep(5)
            gRecord = seqFac.parseGenBank(genomicSeqFileName)
        seqData = gRecord.seq
    
        #---------------------------
        # Analysis
        #---------------------------

        report = Report()
        
        if verbose: print "blasting sequences"
        blastedFeatures = seqTools.seqBlastToFeatures(blastDB, blastExe, targetFile, blastType = "blastn",scoreMin = 1e-5)
        if verbose: print "finished blasting locations"
        if verbose: print "generating recombination oligos"
        (targetMap,report) = recFac.generateTargetingOligos(oligoRecords, blastedFeatures, tagRE = targetTag, boundary = boundary, searchSize = searchSize, cutOff = cutOff) 
        
        seqTools = SequenceTools()
        seqTools.verbose = verbose
        primerReport = seqTools.findPrimers(seqData, targetMap, boundary=primerDistance, oligoSize=primerSize, searchSize=60, targetTm=primerTm)
        report.extend(primerReport)
        #----------------------------
        # Write Report
        #----------------------------
        writer = ReportWriter()
        writer.setFile(oligoFileName)
        writer.write(report)
        writer.closeFile()
        
    #-----------------------------------------
    # do comparison of recombination oligos
    #  to sequecing of strains
    #-----------------------------------------
    
    if mode == "sequence" or (mode == "auto" and recSeqFileName != ''):
        if verbose: print "blasting sequences of recombinations"
        #sequencingFeatures = seqTools.seqBlastToFeatures(blastDB, blastExe, seqResultFile, blastType = "blastn",scoreMin = 1e-5)
        blastRecords = seqTools.seqBlast(blastDB, blastExe, recSeqFileName, blastType = "blastn")
        if verbose: print "finsihed blasting"
        
        oligoReport = fileParser.parseToReport(oligoFileName, keyTag="Row Names", header = ["Row Names","original","genomic_start","genomic_end","genomic_strand","lagging_complement_strand","best"])
        recFac.parseRecombination(oligoReport, blastRecords, logFile=recCompareFileName, minDiff=0, scoreMin=1e-5)
    
    if verbose: print "done"
 
