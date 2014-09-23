'''
Created on Jun 10, 2013
Updated Sept 22, 2014

@author: Graham Rockwell
@summary:
    Method for working through sequencing of recombination strains.
    Focused on comparing recombination oligos used with sequencing results
    Searching multiple file and running analysis in parallel.
    
    !Currently broken and under reconstruction!
'''

from core.reader.FlatFileParser import FlatFileParser
from core.genetic.SequenceTools import BlastTools

from core.util.Config import ReflectionConfig
from core.util.Report import Report
from core.util.ReportWriter import ReportWriter
from util.DataReport import DataReport, DataReportIO

from core.genetic.SequenceTools import SequenceTools, SequenceFactory, RecombinationOligoFactory,ControlRegionTools

from Bio import Entrez,SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

from optparse import OptionParser
import pickle, time, re, csv

def writeToLog(stringValue,logFileH,verbose=False):
    if logFileH == None:
        return None
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

def findTargets(targetReport,feature,minQuality=100,logFile=None):
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
        
        s = "Feature [%s] (%s <-> %s) ===> %s" % (feature.id,start,end,qualityValues)
        writeToLog(s, logFile, True)
        s = "%s [subject]" % (subjectSeq)
        writeToLog(s, logFile, True)
        s = "%s [query]" % (querySeq)
        writeToLog(s, logFile, True)
        s= "%s [vcf]" % (refSeq)
        writeToLog(s, logFile, True)
        s = "%s [alt]" % (readSeq)
        writeToLog(s, logFile, True)
    
    result["query"] = querySeq
    result["subject"] = subjectSeq
    result["refSeq"] = refSeq
    result["alt-read"] = readSeq
    
    return result

def annotateAlignment(targetMap,featuresArray,logFile=None):
    result = {}
    for features in featuresArray:
        for feature in features:
            id = feature.id
            targets = findTargets(targetMap, feature, minQuality= 100, logFile = logFile)
            result[id] = targets
    
    return None
               
def alignVCF(targetFile, vcfFile, targetType = "fasta", verbose = True):
    
    readRecord = parseVCFFile(vcfFile)
    
    targetRecords= SeqIO.parse(open(targetFile), targetType)
    targetRecords = list(targetRecords)

    if verbose: print "Processing [%s] records" % (len(targetRecords))
    
    if verbose: print "blasting sequences"
    blastedFeatures = seqTools.seqBlastToFeatures(blastDB, blastExe, targetFile, blastType = "blastn",scoreMin = 1e-5)
    
    if verbose: print "finished blasting locations"
    readLogName = "read_log_%s.txt" % (readTag)
    logFile = open(readLogName,"w")
    alignmentReport = annotateAlignment(readRecord, blastedFeatures,logFile=logFile)
    logFile.close()

class ProcessVCF:
    '''
    @summary:
    Tools set for dealing with variant call format (vcf) Files.
    Specifically this tools set does the following
    1. Find expanding regions of variants with a designated distance across all samples
    2. Combine variants in each sample into a report object to be used for machine learning analysis of genotype to phenotype.
    3. Create a distance matrix using given variants between given strains. 
    
    '''
    
    def __init__(self):
        self.blastTools = BlastTools()
        self.verbose = False
        
    def writeToLog(self,stringValue,logFileH,verbose=False):
        if logFileH == None:
            return None
        try:
            if verbose: print stringValue
            logFileH.write(stringValue+"\n")
        except:
            print "failed to write to log file [%s]" % (logFileH)
        return None

    def replaceSeqTarget(self,seq,newSeq,loc):
        prefix = seq[:loc]
        post = seq[loc+len(newSeq):]
        result = prefix + newSeq + post
        return result

    def findTargets(self,targetMap,feature,minQuality=10,logFile=None):
        '''
        @var 
        @summary:
        Find variant calls in proximity to given feature.
        '''
        start = feature.location.start.position
        end = feature.location.end.position
        locations = targetMap.keys()
        matches = filter(lambda x: start<=float(x)<=end,locations)
        
        querySeq = feature.qualifiers["query"]
        subjectSeq = feature.qualifiers["subject"]
        refSeq = Seq("_"*len(subjectSeq))
        readSeq = Seq("_"*len(subjectSeq))
        qualityValues = []
        chrom = ''        
        
        dCount = 0    
        for loc in matches:
            floc = float(loc)
            seqStart = int(floc-start)
            
            target = targetMap[loc]
            chrom = target["CHROM"]
            ref = target["REF"]
            alts = target["ALT"].split(",")
            quality = target["QUAL"]
            qualityValues.append((seqStart,quality))
            
            if float(quality) < minQuality:
                print "(failed) Feature [%s] (%s <-> %s) ===> %s" % (feature.id,start,end,qualityValues)
            else:
                dCount += 1
            
            refSeq = self.replaceSeqTarget(refSeq,ref,seqStart)
            for alt in alts:
                readSeq = self.replaceSeqTarget(readSeq,alt,seqStart)
        
        result = {}
        
        result["name"] = feature.id
        result["chrom"] = chrom
        result["start"] = start
        result["end"] = end
        result["quality"] = qualityValues 
        result["query"] = querySeq
        result["subject"] = subjectSeq
        result["refSeq"] = refSeq
        result["alt-read"] = readSeq
        
        return result

    def annotateAlignment(self,targetMap,featuresArray,idTag):
        readLogName = "read_log_%s.txt" % (idTag)
        logFP = open(readLogName,"w")
        result = {}
        
        for features in featuresArray:
            for feature in features:
                id = feature.id
                targets = self.findTargets(targetMap, feature, minQuality= 10, logFile = logFP)
                result[id] = targets
        return result

    def alignVCF(self,targetFile,vcfFile,idTag):
        '''
        @return: ReportObect
        @summary:
        Check alignment to annotated Genomic sequence using BLAST.
        
        '''
        readRecord = parseVCFFile(vcfFile)
    
        targetRecords= SeqIO.parse(open(targetFile), "fasta")
        targetRecords = list(targetRecords)
    
        print "Processing [%s] records" % (len(targetRecords))
        
        if verbose: print "blasting sequences"
        
        self.blastTools.verbose = verbose
        blastedFeatures = self.blastTools.seqBlastToFeatures(blastDB, blastExe, targetFile, blastType = "blastn",scoreMin = 1e-5)
        
        if verbose: print "finished blasting locations"
        alignmentReport = self.annotateAlignment(readRecord, blastedFeatures,idTag)
        
        return alignmentReport

def parseRecombinationOligoFile(filename):
    '''
    Read oligo file, of type used for making recombinant strains.
    
    '''
    header = ["Row Names", "original", "hits", "genomic_start", "genomic_end", "genomic_strand", "sequencing location"]
    fp = open (filename,'rb')
    reader = csv.DictReader(filter(lambda row: row[0]!='#', fp),delimiter="\t",fieldnames=header)
    result = {}
    for l in reader:
        key = l["Row Names"]
        result[key] = l
    return result

def parseVCF(filename, comment_tag="##", delimiter="\t"):
    '''  
    @var filename: name of vcf file
    @filename: String
    @result:[{vcf}]
    
    @summary:
    Parse VCF file to list of dictionary objects
    Split alternative sequences into separate entries
    '''
    
    header = ["#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","unknown"]
    
    fp = open(filename,'rb')
    reader = csv.DictReader(filter(lambda row: row[0:2]!='##', fp),delimiter="\t",fieldnames=header)
    reader.next()
    
    result = []
    for d in reader:
        d["CHROM"] = d["#CHROM"]
        del d["#CHROM"]
        d["QUAL"] = float(d["QUAL"])
        
        alt = d["ALT"].split(",")
        for a in alt:
            d["ALT"] = a
            result.append(d)
        
    return result

def parseVCFFile(filename):
    '''
    @result{position:[{vcf}]
    '''
    reader = parseVCF(filename)
    result = {}
    for d in reader:
        key = d["POS"]
        result[key] = d
    
    return result

def oligoAlignmentAnalysis(oligoSeqFile,vcfFile,sampleCode,blast_tools):
    
    processVCF = ProcessVCF()
    processVCF.blastTools = blast_tools
    vcfAlignment = processVCF.alignVCF(oligoSeqFile,vcfFile,idTag=sampleCode)
        
    alignmentReport = DataReport()
    alignmentLogFileName = workFolder + "%s.log" % (sampleCode)
    logFile = open(alignmentLogFileName,"w")
        
    for (targetKey,alignment) in vcfAlignment.items():
        if verbose: print "alignment [%s]" % (targetKey)
            
        name = alignment["name"]
        irefname = alignment["chrom"]
        quality = alignment["quality"]
        start = alignment["start"]
        end = alignment["end"]
        subject = alignment["subject"]
        query = alignment["query"]
        refSeq = "%s" % alignment["refSeq"]
        altRead = "%s" % alignment["alt-read"]
            
        if verbose: print "[%s]:%s:%s" % (irefname,start,end)
        if irefname == '':
            irefname = refname
            
        readList = []
        alreads = samFile.fetch(irefname,start,end)
        alreads = list(alreads)
        alreadCount = len(alreads)
        
        if len(altRead.replace("_","")) != 0:
            rTag = "%s_count[%s]_quality[%s]" % (targetKey,alreadCount,quality)
            print rTag
            logFile.write(">%s\n" % rTag)
            logFile.write(subject+"[subject]\n") 
            logFile.write(query+"[query]\n")
            logFile.write(refSeq+"[refSeq]\n")
            logFile.write(altRead+"[alt-read]\n")   
        
    logFile.close()
    
    return (alignmentReport,readList)


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
    readRegx = 'Sample\.(.*)\.txt'
    rregx = re.search(readRegx,readFile)
    
    #temporary debugging variables
    #oligoReport = "20120709_FattyAcid_sequencing_colonies_primers.csv"
    #oligoSeqFile = "20120709_FA_recombination_oligos_v2.fasta"  #! Make this something taken from the configuration file
    #oligoSeqFile = options.workFolder + oligoSeqFile    
    
    #bTools = BlastTools()
    #bTools.verbose = verbose
    #bTools.blastDB = blastDB
    #bTools.blastExe = blastExe
    
    if  rregx != None:
        readTag = rregx.group(1).replace(".","_")
    else:
        readTag = "read_analysis"
    
    readRecord = parseVCFFile(readFile)
    
    targetRecords= SeqIO.parse(open(targetFile), "fasta")
    targetRecords = list(targetRecords)

    print "Processing [%s] records" % (len(targetRecords))
    
    if verbose: print "blasting sequences"
    seqTools = SequenceTools()
    seqTools.verbose = verbose
    blastedFeatures = seqTools.seqBlastToFeatures(blastDB, blastExe, targetFile, blastType = "blastn",scoreMin = 1e-5)
    
    if verbose: print "finished blasting locations"
    readLogName = "read_log_%s.txt" % (readTag)
    logFile = open(readLogName,"w")
    alignmentReport = annotateAlignment(readRecord, blastedFeatures,logFile=logFile)
    logFile.close()
    
    print "Done"