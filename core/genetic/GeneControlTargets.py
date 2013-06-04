#!/usr/bin/env python2.6
'''
@author: Graham Rockwell,Church Lab, Harvard Genetics
@updated 8/28/2009
@summary:
Top level main function for finding control gene control regions (promoters and RBS).
Uses flat files of promoter regions to generate lists of sequence of upstream regions of gene optimization targets. 
These targets can then be used for metabolic production optimization.
    
'''

from core.util.ReportWriter import ReportWriter
from core.reader.FlatFileParser import FlatFileParser
from core.genetic.SequenceTools import SequenceFactory, ControlRegionTools

from Bio.SeqRecord import SeqRecord

from optparse import OptionParser
import os
    
if __name__  == "__main__":
    
    parser = OptionParser()
    
    parser.add_option("-v", "--verbose",
                      action="store_true", 
                      dest="verbose", 
                      default=False,
                      help="run in verbose mode")
                  
    parser.add_option("--geneTargets", 
                      dest="geneTargets", 
                      default="",
                      help="Tab delimited gene target file", 
                      metavar="FILE")

    parser.add_option("--promoterTargets", 
                      dest="promoterTargets", 
                      default="",
                      help="Tab delimited promoter region target file", 
                      metavar="FILE")
    
    parser.add_option("-p", "--pickleFile", 
                      dest="pickleFile", 
                      default="",
                      help="python compressed archive file", 
                      metavar="FILE")
                       
    parser.add_option("-o","--output", 
                      dest="outputFileName",
                      help="name of report FILE", 
                      default = "oligo_targets.txt",
                      metavar="FILE")
    
    parser.add_option("-c","--config", 
                      dest="configFileName",
                      help="name of configuration file", 
                      default = "",
                      metavar="FILE")
    
    parser.add_option("--slice", 
                      dest="targetSlice",
                      help="start,end of targets to use", 
                      default = "",
                      metavar="FILE")
    

    #load command line
    (options,args) = parser.parse_args()
    
    verbose = options.verbose
    outputFileName = options.outputFileName  
    geneTargetFileName = options.geneTargets
    promoterTargetFileName = options.promoterTargets
    #sequenceFile = "../../data/Genome/EColi_K12_MG1655.fasta"
    
    if geneTargetFileName == '':
        geneTargetFileName = "../../data/Genome/transcription start promec.csv"
    if promoterTargetFileName == '':
        promoterTargetFileName = "../../data/Genome/ecoli_regulonDB_promoters.txt"
    
    targetReader = FlatFileParser()
    seqFactory = SequenceFactory()
    controlRegionTools = ControlRegionTools()
    
    seqFactory.verbose = verbose
    seqFile = seqFactory.downloadGenBank()
    gRecord = seqFactory.parseGenBank(seqFile)
    genes = seqFactory.extractCDSncRNAFeatures(gRecord)
    
    targetRecord = SeqRecord(id="promoters",seq=gRecord.seq,features = [])
    
    if os.path.exists(geneTargetFileName):
        targetReport = targetReader.parseToReport(geneTargetFileName,"geneID",header=["geneID","start","direction"])
        conversionMap = {"geneID":"gene","F":int(1),"R":int(-1)}
        targetFeatures = seqFactory.convertReportToFeatures(targetReport,conversionMap,"start","start","direction")
        targetRecord.features.extend(targetFeatures)
    
    if os.path.exists(promoterTargetFileName):
        targetReport2 = targetReader.parseToReport(promoterTargetFileName,"PromoterID",header=["PromoterID","Promoter Name","Strand","(+1)","SigmaFactor","Sequence (+1 upper case)","evidence"])
        conversionMap2 = {"promoterName":"gene","forward":int(1),"reverse":int(-1)}
        targetFeatures2 = seqFactory.convertReportToFeatures(targetReport2,conversionMap2,"(+1)","(+1)","Strand")
        targetRecord.features.extend(targetFeatures2)
    
    slice = options.targetSlice
    if slice != '':
        (sStart,sEnd) = slice.split(",")
        genes = genes[sStart:sEnd]

    sdata = str(gRecord.seq).lower()
    #controlRegions = controlRegionTools.localFeatures(genes,targetRecord)   
    seqReport = controlRegionTools.getControlReport(genes,targetRecord,sdata,boundary = 100)
    
    if verbose: print "writing report"
    
    writer = ReportWriter()
    writer.setFile(outputFileName)
    writer.write(seqReport) 
    writer.closeFile()
    
    
    print "done"
    