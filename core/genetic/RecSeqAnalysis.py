'''
Created on Jun 10, 2013

@author: Graham Rockwell
@summary:
    Method for parsing and analyzing sequencing data from recombination engineered strains.
'''
from core.reader.FlatFileParser import FlatFileParser

from core.util.Config import ReflectionConfig
from core.util.Report import Report
from core.util.ReportWriter import ReportWriter

from core.genetic.SequenceTools import SequenceTools, SequenceFactory, RecombinationOligoFactory,ControlRegionTools

from Bio import Entrez,SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

from optparse import OptionParser
import pickle, time

def writeToLog(stringValue,logFileH,verbose=False):
    try:
        if verbose: print stringValue
        logFileH.write(stringValue+"\n")
    except:
        print "failed to write to log file [%s]" % (logFileH)
    return None


def replaceSeqTarget(seq,newSeq,loc):
    prefix = seq[:loc]
    post = seq[loc+len(newSeq):]
    result = prefix + newSeq + post
    return result

def findTargets(targetReport,feature,minQuality=100,logFileName="Target_anlaysis_log.txt"):
    start = feature.location.start.position
    end = feature.location.end.position
    locations = targetReport.getRowNames()
    matches = filter(lambda x: start<=float(x)<=end,locations)
    
    querySeq = feature.qualifiers["query"]
    subjectSeq = feature.qualifiers["subject"]
    refSeq = Seq("_"*len(subjectSeq))
    readSeq = Seq("_"*len(subjectSeq))
    qualityValues = []
    
    dCount = 0
        
    for loc in matches:
        floc = float(loc)
        seqStart = int(floc-start)
        
        target = targetReport.getRow(loc)
        ref = target["REF"]
        alts = target["ALT"].split(",")
        quality = target["QUAL"]
        if float(quality) < minQuality:
            continue
        else:
            dCount += 1
        qualityValues.append([seqStart,quality])
        
        refSeq = replaceSeqTarget(refSeq,ref,seqStart)
        for alt in alts:
            readSeq = replaceSeqTarget(readSeq,alt,seqStart)
    
    result = {}
    if dCount > 0:
        qualityValues.sort()
        
        logFileH = open(logFileName,"w")
        s = "Feature [%s] (%s <-> %s) ===> %s" % (feature.id,start,end,qualityValues)
        writeToLog(s, logFileH, True)
        s = "%s [subject]" % (subjectSeq)
        writeToLog(s, logFileH, True)
        s = "%s [query]" % (querySeq)
        writeToLog(s, logFileH, True)
        s= "%s [vcf]" % (refSeq)
        writeToLog(s, logFileH, True)
        s = "%s [alt]" % (readSeq)
        writeToLog(s, logFileH, True)
        logFileH.close()
    
    result["query"] = querySeq
    result["subject"] = subjectSeq
    result["refSeq"] = refSeq
    result["alt-read"] = readSeq
    
    return result

def annotateAlignment(targetMap,featuresArray,logFileName="log.txt"):
    result = {}
    for features in featuresArray:
        for feature in features:
            id = feature.id
            targets = findTargets(targetMap, feature, minQuality= 100, logFileName = logFileName)
            result[id] = targets
    
    return None
            
            
def parseVCFFile(fileName):
    parser = FlatFileParser(delimiter='\t', comment='##', emptyField='NA')
    header={"#CHROM":"CHROM","POS":"POS","ID":"ID","REF":"REF","ALT":"ALT","QUAL":"QUAL","FILTER":"FILTER","INFO":"INFO","FORMAT":"FORMAT","unknown":"unknown"}
    keyTag = "POS"
    record = parser.parseToReport(fileName, keyTag, header, unique=True)
    
    return record
    

if __name__  == "__main__":
    parser = OptionParser()
    
    parser.add_option("-v", "--verbose",
                      action="store_true", 
                      dest="verbose", 
                      default=False,
                      help="print status messages to stdout")
    
    parser.add_option("-c", "--configFile",
                      dest="config", 
                      default="sequence_analysis.config",
                      help="config file for oligo generation, check config for defaults",
                      metavar="FILE")
    
    parser.add_option("-n", "--configName", 
                      dest="configName", 
                      default="default",
                      help="which configuration setting to use", 
                      metavar="string")
    
    parser.add_option("-g", "--genome_seq",
                      dest="genomeSeq", 
                      default="",
                      help="fasta formated file of recombination target regions",
                      metavar="FILE")
   
    parser.add_option("-t", "--target_seq",
                      dest="targetSeq", 
                      default="",
                      help="fasta formated file of recombination target regions",
                      metavar="FILE")
   
    parser.add_option("-r", "--recombination_seq",
                      dest="recSeq", 
                      default="",
                      help="VCF file format of recombination sequencing results",
                      metavar="FILE") 
    
    #Parse command line arguments, parse configuration, reflect into options
    (options,args) = parser.parse_args()
    configFileName = options.config
    configName = options.configName
    config = ReflectionConfig()    
    config.readfp(open(configFileName))
    config.reflect(configName,options)
    
    #Set variables
    verbose = options.verbose
    blastDB = config.get(configName, "databasefile")
    blastExe = config.get(configName, "blastExe")
    readFile = options.recSeq
    targetFile = options.targetSeq
    
    
    seqTools = SequenceTools()
    seqTools.verbose = verbose
    
    seqFac = SequenceFactory()
    config.reflect(configName,seqFac)
    
    recFac = RecombinationOligoFactory()
    recFac.verbose = verbose
    
    #Parse data files
    #gRecord = seqFac.getGenBankSequence(genbankID=None, filename=None)
    #seqData = str(gRecord.seq).lower()
    
    readRecord = parseVCFFile(readFile)
    targetRecords= SeqIO.parse(open(targetFile), "fasta")
    targetRecords = list(targetRecords)

    print len(targetRecords)
    for t in targetRecords:
        #print t
        pass
    
    if verbose: print "blasting sequences"
    blastedFeatures = seqTools.seqBlastToFeatures(blastDB, blastExe, targetFile, blastType = "blastn",scoreMin = 1e-5)
    
    if verbose: print "finished blasting locations"
    alignmentReport = annotateAlignment(readRecord, blastedFeatures)
        
    print "Done"