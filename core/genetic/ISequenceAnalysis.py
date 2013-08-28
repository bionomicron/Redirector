'''
Created on Aug 11, 2013

@author: Graham Rockwell
@summary:
    Objects and methods for working analysis of sequencing data.
    Particular tools for analysis of recombination strains.
'''

from core.util.Config import ReflectionConfig
from core.reader.FlatFileParser import FlatFileParser 
from core.genetic.SequenceTools import BlastTools
from core.util.Report import Report
from core.util.ReportWriter import ReportWriter

from Bio.Seq import Seq
from Bio import Entrez,SeqIO
from Bio.SeqRecord import SeqRecord

import os,re, csv
import subprocess
from optparse import OptionParser
import difflib

def parseSampleNameTag(name1,name2,regex='.*_([ATCG]+)_.*',tag="Sample_%s"):
    difference = difflib.SequenceMatcher()
    
    difference.set_seqs(name1, name2)
    match = difference.get_matching_blocks()
    if len(match) != 0:
        m = match[0]
        idTag = name1[m[0]:m[0]+m[2]]
    if idTag == '':
        idTag = name1
    
    match = re.match(regex,idTag)
    if match == None:
        print "No match on regex [%s]" %idTag
        return idTag
    
    rtag = match.group(1)
    result = tag % rtag
    return result

def parseRecombinationOligoFile(filename):
    header = ["Row Names", "original", "hits", "genomic_start", "genomic_end", "genomic_strand", "sequencing location"]
    fp = open (filename,'rb')
    reader = csv.DictReader(filter(lambda row: row[0]!='#', fp),delimiter="\t",fieldnames=header)
    result = {}
    for l in reader:
        key = l["Row Names"]
        result[key] = l
    return result
    

def parseRecOligoFile(filename):
    parser = FlatFileParser(delimiter='\t', comment='##', emptyField='NA')
    header={"Row Names":"Row Names", "original":"original", "hits":"hits", "genomic_start":"genomic_start", "genomic_end":"genomic_end", "genomic_strand":"genomic_strand","sequencing location":"sequencing location"}
    keyTag = "Row Names"
    record = parser.parseToReport(filename, keyTag, header, unique=True)
    return record

def parseVCFFile(fileName):
    header = ["#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","unknown"]
    fp = open(fileName,'rb')
    reader = csv.DictReader(filter(lambda row: row[0:2]!='##', fp),delimiter="\t",fieldnames=header)
    result = {}
    for d in reader:
        key = d["POS"]
        result[key] = d
    
    #parser = FlatFileParser(delimiter='\t', comment='##', emptyField='NA')
    #header={"#CHROM":"CHROM","POS":"POS","ID":"ID","REF":"REF","ALT":"ALT","QUAL":"QUAL","FILTER":"FILTER","INFO":"INFO","FORMAT":"FORMAT","unknown":"unknown"}
    #keyTag = "POS"
    #record = parser.parseToReport(fileName, keyTag, header, unique=True)
    
    return result

class SeqReadTools:
    
    def __init__(self):
        # Location of this script. We may find other paths relative to this.
        self.verbose = False
        self.PWD = os.path.dirname(os.path.realpath(__file__ ))
        self.workingDir = "./"
        self.refGenomeDir = "./"
        self.refGenomeTag = ""
        
    def samProcessRead(self,r1,r2,idTag = None):
        difference = difflib.SequenceMatcher()
        
        if not os.path.exists(workFolder):
            if self.verbose: print "making work folder %s" % (workFolder)
            os.mkdir(workFolder)
    
        if idTag == None:
            difference.set_seqs(r1, r2)
            match = difference.get_matching_blocks()
            if len(match) != 0:
                m = match[0]
                idTag = r1[m[0]:m[0]+m[2]]
            if idTag == '':
                idTag = r1
            
        drefTag = self.refGenomeDir + self.refGenomeTag
        dIdTag = self.workingDir + idTag
        
        if not os.path.exists("%s.1.bt2" % drefTag):
            print "building reference index"
            call = "bowtie2-build %s %s > bowtie2_ref_index_log.txt" % (drefTag,drefTag)
            if self.verbose: print "executing command: [%s]" % call
            subprocess.call(call,shell=True)
        else:
            if self.verbose: print "building reference index found %s" % (drefTag)
        
        if self.verbose: print "aligning sequences to reference"
        call = "bowtie2 -q -p 4 -k 3 --local %s -1 %s -2 %s -S %s.sam" % (drefTag,r1,r2,dIdTag)
        if self.verbose: print "executing command: [%s]" % call
        subprocess.call(call,shell=True)
        
        if self.verbose: print "creating bam file"
        call = "samtools view -S -b %s.sam -o %s.bam" % (dIdTag,dIdTag)
        if self.verbose: print "executing command: [%s]" % call
        subprocess.call(call,shell=True)
        
        if self.verbose: print "sorting bam file"
        call = "samtools sort %s.bam %s_s" % (dIdTag,dIdTag)
        if self.verbose: print "executing command: [%s]" % call
        subprocess.call(call,shell=True)
        
        if self.verbose: print "indexing bam file"
        call = "samtools index %s_s.bam" % (dIdTag)
        if self.verbose: print "executing command: [%s]" % call
        subprocess.call(call,shell=True)
        
        resultFile = "%s_s.bam" % (dIdTag)
        
        return (idTag,dIdTag,resultFile)
    
    def vcfProcess(self,bamFile,idTag):
        drefTag = self.refGenomeDir + self.refGenomeTag
        
        if self.verbose: print "processing vcf file"
        call = "freebayes --fasta-reference %s %s -v %s.vcf" % (drefTag,bamFile,idTag)
        if self.verbose: print "executing command: [%s]" % call
        subprocess.call(call,shell=True)
        
        resultFile = "%s.vcf" % idTag
        
        return resultFile

class ProcessVCF:
    
    def __init__(self):
        self.blastTools = BlastTools()
        
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
        start = feature.location.start.position
        end = feature.location.end.position
        locations = targetMap.keys()
        matches = filter(lambda x: start<=float(x)<=end,locations)
        
        querySeq = feature.qualifiers["query"]
        subjectSeq = feature.qualifiers["subject"]
        refSeq = Seq("_"*len(subjectSeq))
        readSeq = Seq("_"*len(subjectSeq))
        qualityValues = []
        
        print "Feature [%s] => Matches [%s]" %(feature.id,len(matches))
        dCount = 0    
        for loc in matches:
            floc = float(loc)
            seqStart = int(floc-start)
            
            target = targetMap[loc]
            ref = target["REF"]
            alts = target["ALT"].split(",")
            quality = target["QUAL"]
            
            if float(quality) < minQuality:
                print "(failed) Feature [%s] (%s <-> %s) ===> %s" % (feature.id,start,end,qualityValues)
            else:
                dCount += 1
            qualityValues.append([seqStart,quality])
            
            refSeq = self.replaceSeqTarget(refSeq,ref,seqStart)
            for alt in alts:
                readSeq = self.replaceSeqTarget(readSeq,alt,seqStart)
        
        result = {}
        if dCount > 0:
            qualityValues.sort()
            
            s = "Feature [%s] (%s <-> %s) ===> %s" % (feature.id,start,end,qualityValues)
            self.writeToLog(s, logFile, True)
            s = "%s [subject]" % (subjectSeq)
            self.writeToLog(s, logFile, True)
            s = "%s [query]" % (querySeq)
            self.writeToLog(s, logFile, True)
            s= "%s [vcf]" % (refSeq)
            self.writeToLog(s, logFile, True)
            s = "%s [alt]" % (readSeq)
            self.writeToLog(s, logFile, True)
        
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
        readRecord = parseVCFFile(vcfFile)
    
        targetRecords= SeqIO.parse(open(targetFile), "fasta")
        targetRecords = list(targetRecords)
    
        print "Processing [%s] records" % (len(targetRecords))
        
        if verbose: print "blasting sequences"
        
        self.blastTools.verbose = verbose
        blastedFeatures = self.blastTools.seqBlastToFeatures(blastDB, blastExe, targetFile, blastType = "blastn",scoreMin = 1e-5)
        
        if verbose: print "finished blasting locations"
        alignmentReport = self.annotateAlignment(readRecord, blastedFeatures)
        
        return alignmentReport
        
        
def filterReadFile(fileName,output,headerRegx,dataRegx):
    f = open(fileName,"r")
    o = open(output,"w")
    header = ''
    for line in f:
        hMatch = re.match(line,headerRegx)
        dMatch = re.match(dataRegx)
        if header != '' & dMatch:
            o.write(header)
            o.write(line)
            header=''
        elif hMatch:
            header = line
        else:
            header = ''
    
    f.close()
    o.close()
            
    return True    


if __name__ == '__main__':
    
    parser = OptionParser()
    
    parser.add_option("-v", "--verbose",
                      action="store_true", 
                      dest="verbose", 
                      default=False,
                      help="set verbose mode")
    
    parser.add_option("-c","--config", 
                      dest="config",
                      default = "sequence_analysis.config",
                      help="name of sequencing analysis configuration file", 
                      metavar="FILE")
    
    parser.add_option("-n", "--configName", 
                      dest="configName", 
                      default="default",
                      help="which configuration setting to use", 
                      metavar="string")
    
    parser.add_option("-m", "--mode", 
                      dest="mode", 
                      default='',
                      metavar="String",
                      help="comma separated list of analysis modes")
    
    parser.add_option("--wd", 
                      dest="workFolder", 
                      default="./",
                      metavar="Directory",
                      help="Directory to store intermediate files")
    
    parser.add_option("-d","--genomic_seq_tag", 
                      dest="genomicSeqTag", 
                      default="NC_000913",
                      metavar="File",
                      help="Tag for genome sequence data")

    parser.add_option("--r1", 
                      dest="readTag1", 
                      default=None,
                      metavar="File",
                      help="readTag for first paired end read")
    
    parser.add_option("--r2", 
                      dest="readTag2", 
                      default=None,
                      metavar="File",
                      help="readTag for first paired end read")


    (options,args) = parser.parse_args()
    configFileName = options.config
    configName = options.configName
    
    config = ReflectionConfig()    
    config.readfp(open(configFileName))
    config.reflect(configName,options)
    
    #Set variables
    verbose = options.verbose
    mode = options.mode.split(",")
    blastDB = config.get(configName, "databasefile")
    blastExe = config.get(configName, "blastExe")
    
    workFolder = options.workFolder
    refTag = options.genomicSeqTag
    readFile1 = options.readTag1
    readFile2 = options.readTag2
    
    #Testing hard coded values
    verbose = True
    mode = ["vcf"]
    refTag = "NC_000913_2"
    oligoReport = "20120709_FattyAcid_sequencing_colonies_primers.csv"
    oligoSeqFile = "20120709_FA_recombination_oligos.fasta"
    oligoSeqFile = options.workFolder + oligoSeqFile
    
    #readFile1 = "LIB007074_GEN00011829_TCCAGT_L001_R1.fastq.bz2"
    #readFile1 = "Sample_TCCAGT_R1.fastq"
    #readFile1 = "TCCAGT_R1_2.fastq"
    #readFile2 = "LIB007074_GEN00011829_TCCAGT_L001_R2.fastq.bz2"
    
    bTools = BlastTools()
    bTools.verbose = verbose
    bTools.blastDB = blastDB
    bTools.blastExe = blastExe
    
    sampleCode = parseSampleNameTag(readFile1, readFile2, regex='.*_([ATCG]+)_.*',tag="Sample_%s")
    print "Using Sample Code [%s]" % (sampleCode)
    
    if "hit_count" in mode:
        oligoRecords = parseRecOligoFile(oligoReport)
        targetMap = {}
        for rowName in oligoRecords.getRowNames():
            start = oligoRecords.getElement(rowName, "genomic_start")
            end = oligoRecords.getElement(rowName, "genomic_end")
            middle = (float(start)+float(end))/2.0
            location = oligoRecords.getElement(rowName, "sequencing location")
            targetMap[rowName] = (float(start),float(end),float(location))
            print "[%s] => ([%s] <-> [%s]) [%s]" % (rowName,start,end,location)
                   
        targetCount = bTools.blastAlignSequences(targetMap, readFile1,readType="fastq")
    
        blastReport = Report()
        blastReport.addColumnHash("blast_hit_count", targetCount)
    
        writer = ReportWriter()
        writer.setFile("Blast_hit_count.csv")
        writer.write(blastReport)
    
    if "vcf" in mode:
        
        vcfFile = "%s.vcf" % (sampleCode)
        sTools = SeqReadTools()
        sTools.verbose = True
        sTools.refGenomeTag = refTag
        sTools.workingDir = workFolder
        sTools.refGenomeDir = workFolder
        sTools.refGenomeTag = refTag
        #(id,dID,bamFile) = sTools.samProcessRead(readFile1,readFile2,idTag=sampleCode)
        #vcfFile = sTools.vcfProcess(bamFile, idTag = id)
        
        processVCF = ProcessVCF()
        processVCF.blastTools = bTools
        vcfAlignment = processVCF.alignVCF(oligoSeqFile,vcfFile,idTag=sampleCode)
        
        #vcfReportName = "%s_vcf_report.csv" % (sampleCode)
        #writer = ReportWriter()
        #writer.setFile(vcfReportName)
        #writer.write(vcfAlignmentReport)
    
    print "done"   
    #print "finding VCF information"
    #call = "freebayes --fasta-reference %s %s_s.bam -v %s_%s.vcf" % (refTag,dreadTag,refTag,readTag)
    #print call
    #subprocess.call(call,shell=True)
    

