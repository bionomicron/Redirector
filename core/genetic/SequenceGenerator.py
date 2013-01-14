#!/usr/bin/env python2.6

#---------------------------------
#Graham Rockwell
#Church Group Harvard Genetics
#Last updated 8/28/2009
#---------------------------------

from core.reader.FlatFileParser import FlatFileParser
from core.reader.SequenceReader import SequenceReader

from core.writer.Report import Report
from core.writer.ReportWriter import ReportWriter

from core.genetic.SequenceTools import SequenceFactory, ControlRegionTools

from Bio.SeqRecord import SeqRecord

from optparse import OptionParser
import pickle


    
def main_function():
    '''
    Uses datasets of promoter regions to generate list of
    optimization targets for 
    diverstiy of genetric control for metabolic production optimization.
    '''
    
    parser = OptionParser()
    
    parser.add_option("-v", "--verbose",
                      action="store_true", 
                      dest="verbose", 
                      default=False,
                      help="run in verbose mode")
                  
    parser.add_option("-s", "--sequencefile", 
                      dest="sequenceFile", 
                      default="",
                      help="FASTA file of squence", 
                      metavar="FILE")
    
    parser.add_option("-t", "--targetfiles", 
                      dest="targetFile", 
                      default="",
                      help="comma seperated list of tab delimited target files, listed in configuration", 
                      metavar="FILE")
    
    parser.add_option("-p", "--pickleFile", 
                      dest="pickleFile", 
                      default="",
                      help="python compressed archive file", 
                      metavar="FILE")
                       
    parser.add_option("-o","--outputfile", 
                      dest="outputFileName",
                      help="name of report FILE", 
                      default = "output.txt",
                      metavar="FILE")
    
    parser.add_option("-c","--config", 
                      dest="configFileName",
                      help="name of configuration file", 
                      default = "",
                      metavar="FILE")

    (options,args) = parser.parse_args()
    
    #load command line
    
    verbose = options.verbose
    sequenceFile = options.sequenceFile
    outputFileName = options.outputFileName  
    pickleFileName = options.pickleFile
    configFileName = options.configFileName
    
    sequenceFile = "../../data/Genome/EColi_K12_MG1655.fasta"
    targetFile = "../../data/Genome/transcription start promec.csv"
    targetFile2 = "../../data/Genome/ecoli_regulonDB_promoters.txt"
    
    sequenceReader = SequenceReader()
    targetReader = FlatFileParser()
    seqFactory = SequenceFactory()
    controlRegionTools = ControlRegionTools()

    targetReport = targetReader.parseToReport(targetFile,"geneID",header=["geneID","start","direction"])
    conversionMap = {"geneID":"gene","F":int(1),"R":int(-1)}
    targetFeatures = seqFactory.convertReportToFeatures(targetReport,conversionMap,"start","start","direction")
    
    #targetReport2 = targetReader.parseToReport(targetFile2,"PromoterID",header=["PromoterID","Promoter Name","Strand","(+1)","SigmaFactor","Sequence (+1 upper case)","evidence"])
    #conversionMap2 = {"promoterName":"gene","forward":int(1),"reverse":int(-1)}
    #targetFeatures2 = seqFactory.convertReportToFeatures(targetReport2,conversionMap2,"(+1)","(+1)","Strand")
    
    #tf = []
    #for f in targetFeatures2:
        #if f.location.start.position != 0:
            #tf.append(f)
    #targetFeatures2 = tf
    
    seqFactory.verbose = verbose
    seqFile = seqFactory.downloadGenBank()
    gRecord = seqFactory.parseGenBank(seqFile)
    genes = seqFactory.extractCDSncRNAFeatures(gRecord)
    
    targetRecord = SeqRecord(id="promoters",seq=gRecord.seq,features = targetFeatures)
    #targetRecord.features.extend(targetFeatures2)
    
    #genes = genes[1:100]

    sdata = str(gRecord.seq).lower()
    controlRegions = controlRegionTools.localFeatures(genes,targetRecord)   
    seqReport = controlRegionTools.getControlReport(genes,targetRecord,sdata,boundary = 100)
    
    if verbose: print "writing report"
    
    writer = ReportWriter()
    writer.setFile("oligo_targets_all.txt")
    writer.write(seqReport) 
    writer.closeFile()
    
    
    print "done"
    
if __name__  == "__main__":
    main_function()
    