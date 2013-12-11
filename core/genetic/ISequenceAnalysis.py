'''
Created on Aug 11, 2013

@author: Graham Rockwell
@summary:
    Objects and methods for working analysis of sequencing data.
    Particular tools for analysis of recombination strains.
'''

from util.Config import ReflectionConfig 
from core.genetic.SequenceTools import BlastTools
from util.Report import Report
from util.ReportWriter import ReportWriter
from util.FlatFileParser import FlatFileParser

from Bio.Seq import Seq
from Bio import Entrez,SeqIO
from Bio.SeqRecord import SeqRecord

import os,re, csv, pickle
import subprocess
from optparse import OptionParser
import difflib
import pysam

class VarianceTools:
    '''
    Tools for processing high throughput sequencing assembly and analysis.
    '''
    
    def __init__(self):
        self.verbose = False
        
        # Location pof this script. We may find other paths relative to this.
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
        #if True:
            print "building reference index"
            call = "bowtie2-build %s %s > bowtie2_ref_index_log.txt" % (drefTag,drefTag)
            if self.verbose: print "executing command: [%s]" % call
            subprocess.call(call,shell=True)
        else:
            if self.verbose: print "building reference index found %s" % (drefTag)
        
        if self.verbose: print ">aligning sequences to reference"
        call = "bowtie2 -q -p 4 -k 3 --local %s -1 %s -2 %s -S %s.sam" % (drefTag,r1,r2,dIdTag)
        if self.verbose: print "executing command: [%s]" % call
        subprocess.call(call,shell=True)
        
        if self.verbose: print ">creating bam file"
        call = "samtools view -S -b -o %s.bam %s.sam" % (dIdTag,dIdTag)
        if self.verbose: print "executing command: [%s]" % call
        subprocess.call(call,shell=True)
        
        if self.verbose: print ">sorting bam file"
        call = "samtools sort %s.bam %s_s" % (dIdTag,dIdTag)
        if self.verbose: print "executing command: [%s]" % call
        subprocess.call(call,shell=True)
        
        if self.verbose: print ">indexing bam file"
        call = "samtools index %s_s.bam" % (dIdTag)
        if self.verbose: print "executing command: [%s]" % call
        subprocess.call(call,shell=True)
        
        resultFile = "%s_s.bam" % (dIdTag)
        
        return (idTag,dIdTag,resultFile)
    
    def vcfProcess(self,bamFile,idTag):
        drefTag = self.refGenomeDir + self.refGenomeTag
        resultFile = self.workingDir + "%s.vcf" % idTag
        
        if self.verbose: print "processing vcf file"
        call = "freebayes --fasta-reference %s %s -v %s" % (drefTag,bamFile,resultFile)
        if self.verbose: print "executing command: [%s]" % call
        subprocess.call(call,shell=True)
        
        return resultFile
    
class ProcessVCF:
    
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

def parseSampleNameTag(name1,name2,regex='.*_([ATCG\-]+)_.*',tag="Sample_%s"):
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
        print "No match on regex [%s] in  [%s]" %(regex, idTag)
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

def parseVCF(filename):
    '''
    @var filename: name of vcf file
    @filename: String
    @result:[{vcf}]
    @summary:
    parse VCF file to list of dicts
    split alt sequences into seperate entries
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
    reader = parseVCF(filename)
    result = {}
    for d in reader:
        key = d["POS"]
        result[key] = d
    
    return result
    
def joinVariants(vcfData,vcfCollection={},strainID='None',joinDistance=50):
    
    for vcf in vcfData:
        vcf["Strain"] = strainID
        vcf["ID"] = (strainID,vcf["POS"],vcf["ALT"])
        
        start = float(vcf["POS"])
        end = start + len(vcf["REF"])
        fstart = start
        fend = end
        
        print "[%s] (%s,%s)" % (vcf["Strain"],start,end)
        f = lambda x: (fstart-joinDistance) <= x <= (fend + joinDistance)
        g = lambda x: f(x[0]) or f(x[1])
        matches = filter(g,vcfCollection.keys())
        
        result = {vcf["ID"]:vcf}
        
        while len(matches) > 0:
            for loc in matches:
                (istart,iend) = loc
                if istart < fstart: fstart = istart
                if iend < fend: fend = iend
                vcfValues = vcfCollection[loc]
                if vcfValues != None:
                    result.update(vcfValues)
                del vcfCollection[loc]
            matches = filter(g,vcfCollection.keys())
        iLoc = (fstart,fend)
        vcfCollection[iLoc] = result                
                
    return vcfCollection

def joinStrainVariants(vcfMasterReport,vcfCollection=None,joinDistance=50):
    print "joining variants"
    strainIDs = vcfMasterReport.keys()
    vcfCollection = {}
    for strainID in strainIDs:
        vcfData = vcfMasterReport[strainID]
        joinVariants(vcfData,vcfCollection=vcfCollection,strainID=strainID,joinDistance=joinDistance)
    return vcfCollection

def vcfCollectionReport(vcfCollection):
    result = Report()
    
    for (loc,vcfs) in vcfCollection.items():
        (start,end) = loc
        if vcfs == None: continue
        
        result.add(loc,"position",start)
        
        row = {}
        for (id,vcf) in vcfs.items():
            strain = vcf["Strain"]
            if strain not in row.keys():
                row[strain] = ''
            row[strain] = row[strain] + "(%s,%s,%s)" % (vcf["POS"],vcf["READS"],vcf["ALT"])
        
        result.add(loc,"count",len(row)) 
        for (strain,variantString) in row.items():
            result.add(loc,strain,variantString)
          
    return result
        
        
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

def buildVCF(workFolder,sampleCode,mode):
        
        bamFile = workFolder + "%s_s.bam" % (sampleCode)
        vcfFile = workFolder + "%s.vcf" % (sampleCode)
        
        if not os.path.exists(bamFile) or "rebuild" in mode:
            print "Building bam file %s" % (bamFile)
            (id,dID,bamFile) = sTools.samProcessRead(readFile1,readFile2,idTag=sampleCode)
        else:
            print "%s found" % (bamFile)
                        
        if not os.path.exists(vcfFile) or "rebuild" in mode:
        #if True:
            print "Building vcf file %s" % (vcfFile)    
            vcfFile = sTools.vcfProcess(bamFile, idTag = sampleCode)
        else:
            print "%s found" % (vcfFile)
        
        return(bamFile,vcfFile)

def vcfAnalysis(vcfData, samFile, strainID = '', refName='', minQuality = 30, minCount = 1):
    vReport = Report()
    
    for vcf in vcfData:
        vcf["STRAIN"] = strainID
        irefname = vcf["CHROM"]
        start = int(vcf["POS"])
        ref = vcf["REF"]
        alt = vcf["ALT"]
        quality = vcf["QUAL"]
        end = start + len(vcf)
        vcf["END"] = end
        
        alreads = samFile.fetch(irefname,start,end)
        alreads = list(alreads)
        rCount = len(alreads)
        
        vcf["READS"] = rCount
        
        if float(quality < minQuality) or rCount < minCount:
            continue
        
        if verbose: print "[%s]:(%s:%s)x[%s]=[%s] : [%s] ==> [%s]" % (irefname,start,end,rCount,quality,ref,alt)
        if irefname == '':
            irefname = refname
        
        vID = "%s_%s" % (strainID,start)
        vReport.add(vID,"strainID", strainID)
        vReport.add(vID,"start", start)
        vReport.add(vID,"end", end)
        vReport.add(vID,"chrom",irefname)
        vReport.add(vID,"count",rCount)
        vReport.add(vID,"qual",quality)
        vReport.add(vID,"ref",ref)
        vReport.add(vID,"alt",alt)
        vReport.add(vID,"bamFile",bamFile)
        vReport.add(vID,"bamFile",samFile)
    
    return (vcfData, vReport)

def oligoAlignmentAnalysis(oligoSeqFile,vcfFile,sampleCode):
    
    processVCF = ProcessVCF()
    processVCF.blastTools = bTools
    vcfAlignment = processVCF.alignVCF(oligoSeqFile,vcfFile,idTag=sampleCode)
        
    alignmentReport = Report()
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
                      default='vcf',
                      metavar="String",
                      help="comma separated list of analysis modes")
    
    parser.add_option("-w", 
                      dest="workFolder", 
                      default="./",
                      metavar="Directory",
                      help="Directory to store intermediate files")
    
    parser.add_option("-g", 
                      dest="genomeFolder", 
                      default="./",
                      metavar="Directory",
                      help="Directory to for genomic files")
    
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
    genomeFolder = options.genomeFolder
    refTag = options.genomicSeqTag
    readFile1 = options.readTag1
    readFile2 = options.readTag2
    
    #Quality values, may be taken from configuration later
    minQuality = 30.00
    minCount = 5 
    
    #Testing hard coded values
    refTag = "NC_000913_2"
    #oligoReport = "20120709_FattyAcid_sequencing_colonies_primers.csv"
    #oligoSeqFile = "20120709_FA_recombination_oligos_v2.fasta"  #! Make this something taken from the configuration file
    #oligoSeqFile = options.workFolder + oligoSeqFile    
    
    bTools = BlastTools()
    bTools.verbose = verbose
    bTools.blastDB = blastDB
    bTools.blastExe = blastExe
    
    print "ISeuqence Analysis Version 1.2"
    print "Running Mode %s" % (mode)
    
    idCode = parseSampleNameTag(readFile1, readFile2, regex='.*_([ATCG\-]+)_.*',tag="%s")
    sampleCode = "Sample_%s" % idCode
    if verbose: print "Using Sample Code [%s]" % (sampleCode)
    
    currentPath = os.getcwd()
    pathRegex = '.*GEN([0-9].+)'
    match = re.match(pathRegex,currentPath)
    
    if match == None:
        print "No match on regular expression [%s] in  [%s]" %(pathRegex, currentPath)
        pathCode = "Zero"
    else:
        pathCode = match.group(1)
    if verbose: print "Path Code [%s]" % (pathCode)
    
    #File Names
    varianceRepositoryFile = workFolder + "Variant_Repository.bac"
    
    masterReportFile = workFolder + "Master_Report.csv"
    mReportFile = workFolder + "Master_Report.bac"
    varianceReportFile = workFolder + "Variance_Report.csv"
    mVarianceLog = workFolder + "Master_Var.bac"
    
    #Strain Identifier
    strainID = "%s_%s" % (sampleCode,pathCode)  
    
    #Variance Tools.
    sTools = VarianceTools()
    sTools.verbose = True
    sTools.refGenomeTag = refTag
    sTools.workingDir = workFolder
    sTools.refGenomeDir = genomeFolder
    sTools.refGenomeTag = refTag
     
    #========================
    # Construct variant data
    #======================== 
        
    if os.path.exists(mReportFile):
        mReportFh = open(mReportFile,"r")
        masterReport = pickle.load(mReportFh)
        mReportFh.close()
    else:
        masterReport = {}
        
    if os.path.exists(varianceRepositoryFile):
        fh = open(varianceRepositoryFile,"r")
        varRepository = pickle.load(fh)
        fh.close()
    else:
        varRepository = {}        
        
    #if 'vcf' not in mode:
    #if strainID in masterReport.getColumnNames()
    #if ('vcf' not in mode) or (strainID in masterReport.keys() ):
    if False:
        print "VCF already in Master Report"
    else:
        print "Process VCF to Report for sample code [%s]" % (sampleCode)
        
        (bamFile,vcfFile) = buildVCF(workFolder, sampleCode, mode)
        
        print "Opening bam file [%s]" % (bamFile)
        samFile = pysam.Samfile(bamFile,"rb")
        refname = samFile.getrname(0)
            
        vcfData = parseVCF(vcfFile)
        
        print "Building Variant information from VCF report"
        strainID = "%s_%s" % (sampleCode,pathCode)     
        (vcfData, vReport) = vcfAnalysis(vcfData=vcfData, samFile=samFile, strainID = strainID, refName=refname, minQuality=minQuality, minCount=minCount)
        #masterReport.extend(vReport)
        varRepository[strainID] = vcfData
        masterReport[strainID] = vReport
        
 
        #Update data file of reported variants
        print "Building VCF master pickle %s" % (mReportFile)
        mReportFh = open(mReportFile,"w")
        pickle.dump(masterReport,mReportFh)
        mReportFh.close()
        
        #Update data file of reported variants
        print "Building VCF master pickle %s" % (varianceRepositoryFile)
        FH = open(varianceRepositoryFile,"w")
        pickle.dump(varRepository,FH)
        FH.close()
          
    #if 'join' in mode:
        
        if os.path.exists(mVarianceLog):
            tempFH = open(mVarianceLog,"r")
            varCollection = pickle.load(tempFH)
            tempFH.close()
        else:
            varCollection = {}   
    
        #Join variants that are close together
        strainID = "%s_%s" % (sampleCode,pathCode)        
        #variantCollection = joinVariants(vcfData, joinDistance=100, vcfCollection=varCollection,strainID = strainID)
        varCollection = joinStrainVariants(varRepository, varCollection, joinDistance=100)
        varGroupReport = vcfCollectionReport(varCollection)
            
        #Update data file of joined variants
        print "Building Variance Collection %s" % (mVarianceLog)
        varFH = open(mVarianceLog,"w")
        pickle.dump(varCollection,varFH)
        varFH.close()
        
    #Write flat file of vcf master report and joined variant matrix
    if "vReport" in mode:
            print "Building VCF Master Report"
            writer = ReportWriter()
            writer.setFile(masterReportFile)
            writer.write(masterReport)
            writer.closeFile()
            
            print "Building Variance Array"
            writer = ReportWriter()
            writer.setFile(varianceReportFile)
            writer.write(varGroupReport)
            writer.closeFile()
    
    print "Done"