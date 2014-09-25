'''
Created on Aug 11, 2013
Updated Sep 22, 2014

@author: Graham Rockwell
@summary:
    Objects and methods for working analysis of sequencing data.
    Particular tools for analysis of recombination strains.
'''

#Redirector methods and objects
from util.Config import ReflectionConfig 
from util.DataReport import DataReport, DataReportIO
from core.genetic.SequenceTools import SequenceFactory

#Bio Python Objects
from Bio.Seq import Seq
from Bio import Entrez,SeqIO
from Bio.SeqRecord import SeqRecord
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist

#System objects and methods
import os, re, csv, pickle
import subprocess, difflib

#Python samtools
import pysam

#Command line parsing
from optparse import OptionParser

class SequenceAssemblyTools:
    '''
    @summary:
    Tools for processing high throughput sequencing assembly and analysis.
    '''
    
    def __init__(self):
        self.verbose = False
        self.PWD = os.path.dirname(os.path.realpath(__file__ ))
        self.workingDir = "./"
        self.refGenomeDir = "./"
        self.refGenomeTag = ""
        self.reg_genome = self.refGenomeDir + self.refGenomeTag
        self.analysis_tools_path= ""
        
    def samProcessRead(self,r1,r2,idTag = None, rebuild=True):
        '''
        @summary:
        Takes in paired sequencing files and process them into index bam files
        '''
        
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
            
        reference_tag = self.refGenomeDir + self.refGenomeTag
        dIdTag = self.workingDir + idTag
        
        #Check Path
        print os.environ["PATH"]
        if self.analysis_tools_path not in os.environ["PATH"]:
            os.environ["PATH"] += self.analysis_tools_path
        
        print "Analysis tools must be in current path [%s]" % os.path.expandvars("$PATH")
        if not os.path.exists("%s.1.bt2" % reference_tag):
        #if True:
            print "building reference index"
            call = "bowtie2-build %s %s > bowtie2_ref_index_log.txt" % (reference_tag,reference_tag)
            if self.verbose: print ">%s" % call
            subprocess.call(call,shell=True)
        else:
            if self.verbose: print "building reference index found %s" % (reference_tag)
        
        if self.verbose: print ">aligning sequences to reference"
        #try BWA script.
        if os.path.exists("%s.sam" % dIdTag) and not rebuild:
            if self.verbose: print "found %s.sam" % dIdTag
        else: 
            print "building %s.sam" % dIdTag
            call = "bowtie2 -q -p 4 -k 3 --local %s -1 %s -2 %s -S %s.sam" % (reference_tag,r1,r2,dIdTag)
            if self.verbose: print ">%s" % call
            subprocess.call(call,shell=True)
        
        if os.path.exists("%s.bam" % dIdTag) and not rebuild:
            print "found %s.bam" % dIdTag
        else:
            print "Building %s.bam" % dIdTag
            call = "samtools view -S -b -o %s.bam %s.sam" % (dIdTag,dIdTag)
            if self.verbose: print ">%s" % call
            subprocess.call(call,shell=True)
        
        if os.path.exists("%s_s.bam" % dIdTag) and not rebuild:
            print "found %s.bam" % dIdTag
        else:
            print "Building %s_s.bam" % dIdTag
            call = "samtools sort %s.bam %s_s" % (dIdTag,dIdTag)
            if self.verbose: print ">%s" % call
            subprocess.call(call,shell=True)
        
            if self.verbose: print "Indexing bam file"
            call = "samtools index %s_s.bam" % (dIdTag)
            if self.verbose: print ">%s" % call
            subprocess.call(call,shell=True)
        
        if os.path.exists("%s.fq" % dIdTag) and not rebuild:
            print "found %s.fq" % dIdTag
        else:
            if self.verbose: print "building fasta sequence"
            call = "samtools mpileup -uf %s %s_s.bam | bcftools view -cg - | vcfutils.pl vcf2fq > %s.fq" % (reference_tag,dIdTag,dIdTag)
            if self.verbose: print ">%s" % call
            subprocess.call(call,shell=True)
            
        bam_file = "%s_s.bam" % (dIdTag)
        fasta_file = "%s.fq" % (dIdTag)
        
        return (idTag,dIdTag,bam_file,fasta_file)
    
    def vcfProcess(self,bamFile,idTag):
        '''
        @summary:
        Creates variant call format (vcf) from index bam file.
        Currently uses freebayes program, this may not be the best.
        Currently focused on E. coli sequencing data.
        '''
        
        reference_tag = self.refGenomeDir + self.refGenomeTag
        resultFile = self.workingDir + "%s.vcf" % idTag
        
        if self.verbose: print "processing vcf file"
        call = "freebayes --fasta-reference %s %s -v %s" % (reference_tag,bamFile,resultFile)
        if self.verbose: print "executing command: [%s]" % call
        subprocess.call(call,shell=True)
        
        return resultFile
    
    def buildAlignment(self,workFolder,sampleCode,mode):
        '''
        @summary:
        Pipeline to process paired sequencing files into index bam files and vcf (variant call format) file.
        '''
        
        bam_file = workFolder + "%s_s.bam" % (sampleCode)
        vcf_file = workFolder + "%s.vcf" % (sampleCode)
        fasta_file = workFolder + "%s.fq" % (sampleCode)
            
        if not os.path.exists(bam_file) or "rebuild2" in mode:
            print "Building bam file %s" % (bam_file)
            (id,dID,bam_file,fasta_file) = self.samProcessRead(readFile1,readFile2,idTag=sampleCode)
        else:
            print "%s found" % (bam_file)
                            
        if not os.path.exists(vcf_file) or "rebuild2" in mode:
        #if True:
            print "Building vcf file %s" % (vcf_file)    
            vcf_file = sTools.vcfProcess(bam_file, idTag = sampleCode)
        else:
            print "%s found" % (vcf_file)
            
        return(bam_file,vcf_file,fasta_file)

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
        d["START"] = start = float(d["POS"])
        d["END"] = start + len(d["REF"])
        
        alt = d["ALT"].split(",")
        for a in alt:
            d["ALT"] = a
            result.append(d)
        
    return result

def parseVcfToMap(filename):
    '''
    @result{position:[{vcf}]
    '''
    reader = parseVCF(filename)
    result = {}
    for d in reader:
        key = d["POS"]
        result[key] = d
    
    return result

class VariantRegionAnalysis:
    '''
    Set of methods for finding regions of variants in a set of samples.
    Then combining variants for each sample in those regions.
    Useful for finding areas of high recombination / mutation activity.
    Also used for preparing matrix of variance regions for machine learning.
    '''
    
    def __init__(self):
        pass

    def getVariantRegions(self,vcfData,strainID,verbose=False):
        '''
        @summary:
        Collect region (start, end) sequence data from list of variants.
        '''
        regions = []
        
        for vcf in vcfData:
            vcf["Strain"] = strainID
            vcf["ID"] = (strainID,vcf["POS"],vcf["ALT"])
            
            start = float(vcf["POS"])
            end = start + len(vcf["REF"])
            fstart = start
            fend = end
            
            if verbose: print "[%s] (%s,%s)" % (vcf["Strain"],start,end)
            regions.append((fstart,fend))
            
        return regions

    def _filterRegions(self,regions,start,end,join_distance):
        '''
        @summary:
        return list of regions with the join distance from the given start and end locations
        '''
        f = lambda x: (start-join_distance) <= x <= (end + join_distance)
        g = lambda x: f(x[0]) or f(x[1])
        matches = filter(g,regions)
        return matches

    def joinVariantRegions(self,regions,join_distance=50,verbose=False):
        '''
        @return: [regions]
        @summary:
        Combine regions with in join_distance
        '''
        regions.sort()
        
        for region in regions:
            fstart = region[0]
            fend = region[1]
            
            matches = self._filterRegions(regions,fstart,fend,join_distance)
            
            if len(matches) == 1: continue
            
            #if verbose: print "region matches [%s]" % (len(matches))
            
            while len(matches) > 0:
                locations = []
                for loc in matches:
                    locations.extend(loc)
                    regions.remove(loc)
                fstart = min(locations)
                fend = max(locations)
                matches = self._filterRegions(regions,fstart,fend,join_distance)
                
            iLoc = (fstart,fend)
            regions.append(iLoc)
            
            #if verbose: print "final region (%s,%s)" % (fstart,fend)
            
        return regions

    def findVariantRegions(self,vcfData,regions=[],strainID='',joinDistance=50,verbose=False):
        '''
        @summary:
        Build list of regions in which variants are found.
        Used to build a set of Genome comparable regions.
        Build regions by join close by locations of variants
        '''
        for vcf in vcfData:
            vcf["Strain"] = strainID
            vcf["ID"] = (strainID,vcf["POS"],vcf["ALT"])
            
            start = float(vcf["POS"])
            end = start + len(vcf["REF"])
            fstart = start
            fend = end
            
            if verbose: print "[%s] (%s,%s)" % (vcf["Strain"],start,end)
            f = lambda x: (fstart-joinDistance) <= x <= (fend + joinDistance)
            g = lambda x: f(x[0]) or f(x[1])
            matches = filter(g,regions)
            
            while len(matches) > 0:
                for loc in matches:
                    (istart,iend) = loc
                    if istart < fstart: fstart = istart
                    if iend > fend: fend = iend
                    regions.remove(loc)
                matches = filter(g,regions)
            iLoc = (fstart,fend)
            regions.append(iLoc)  
            
        return regions              
            
    def joinVariants(self,vcfData,regions,samData,verbose=False):
        '''
        @summary:
        Take selected regions of a selected Genome.
        Groups variants in selected regions
        Return region grouped variants.
        '''
        irefname = samData.getrname(0)
        vcfRegions = {}
        
        
        for region in regions:
            start = region[0]
            end = region[1]
        
            alreads = samData.fetch(irefname,start,end)
            alreads = list(alreads)
            rCount = len(alreads)
        
            f = lambda x: (float(start) <= float(x["POS"]) <= float(end))
            
            matches = filter(f,vcfData)
            l = len(matches)
            if l > 1:
                #print "[%s],(%s)" % (l,matches)
                pass
            vcfRegion = {"Count":rCount}
            vcfRegion["vcfData"] = matches
                
            vcfRegions[region] = vcfRegion
            
        return vcfRegions
            
    def joinStrainVariants(self,varianceData,joinDistance=50,verbose=False):
        '''
        @return: {[strain_id]:{region:[variants]}}
        @summary:
        Merge near by variants into variant regions for better analysis
        '''
        
        strainIDs = varianceData.keys()
        vcfCollection = {}
        regions = []
        
        if verbose: print "Populating variance region list with initial variants"
        for strainID in strainIDs:
            vcfData = varianceData[strainID]["vcfData"] 
            #regions = findVariantRegions(vcfData, regions, strainID, joinDistance, verbose)
            iregions = self.getVariantRegions(vcfData, strainID, verbose)
            regions.extend(iregions)
        
        if verbose: print "Joining variants in selected variance regions with in %s"  % (joinDistance)    
        regions = self.joinVariantRegions(regions, joinDistance, verbose=True)
        
        if verbose: print "Using regions to join sets of close variants into report"
        for strainID in strainIDs:
            print "Joining variants for [%s]" % (strainID)
            vcfData = varianceData[strainID]["vcfData"]
            bamFile = varianceData[strainID]["sourceFile"]
            samData = pysam.Samfile(bamFile,"rb")
            vcfRegions = self.joinVariants(vcfData,regions,samData)
            vcfCollection[strainID] = vcfRegions
            
        print "Variant report complete"
        return vcfCollection
    
    def vcfCollectionReport(self,vcfCollection,vcf,reorder=True,useCount=True,fill_blank="NA"):
        '''
        @summary:
        Build report object form variant calls groups by regions into "variance collection".
        Report is made to be easily written to delimited matrix / spread sheet format.
        '''
        result = DataReport()
        coverageReport = DataReport()
        data = vcfCollection.items()
        data.sort()
        
        strainIDs = vcfCollection.keys()
        #sort strain IDs
        
        #if reorder:
        if True:
            print "reordering report"
            columnNames = strainIDs
            cmap = {}
            regex = ".*_([0-9]+).*"
            for cName in columnNames:
                match = re.match(regex,cName)
                if match != None:
                    key = match.group(1)
                else:
                    key = cName
                cmap[key] = cName
            strainKeys = cmap.keys()
            strainKeys.sort()
         
        #for (strainID,vcfRegions) in vcfCollection.items():
        for strainKey in strainKeys:
            strainID = cmap[strainKey]
            vcfRegions = vcfCollection[strainID]
            
            print "Collecting [%s] regions" % (strainID)
            for (loc,vcfRegion) in vcfRegions.items():
                
                count = vcfRegion["Count"]
                vcfData = vcfRegion['vcfData']
                item = "[%s]:" % (count)
                coverageReport.add(loc,strainID,count)
                    
                if useCount:
                    if result.get(loc,"Region_Count") == None:
                        #result.add(loc,"Region_Count",0)
                        rCount = 0
                    else:
                        rCount = result.get(loc,"Region_Count")
                        rCount = float(rCount)
                        
                    if len(vcfData) != 0:
                        rCount += 1
                        #result.add(loc,"Region_Count",rCount) 
                
                for vcf in vcfData:
                    item = item + "(%s,%s,%s):" % (vcf["POS"],vcf["READS"],vcf["ALT"])
                if len(vcfData) != 0:
                    result.add(loc,strainID,item)
                else:
                    result.add(loc,strainID,fill_blank)
    
        return (result,coverageReport)

def parseSampleNameTag(name1,name2,regex='.*_([ATCG\-]+)_.*',tag="Sample_%s"):
    '''
    @summary:
    Process sample file names to find important identificaiton codes
    '''
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
               
def processVariantData(vcf_data, bam_file, strainID = '', refName='', minQuality = 30, minCount = 1,minSize=0):
    '''
    @return [vcf]
    @summary:
    Process variant data link combine with related sequencing reference.
    Also creates a report of all variant data.
    '''
    #vReport = Report() #Testing new report object
    vReport = DataReport()
    vcfDataResult = []
    
    print "Opening bam file [%s]" % (bam_file)
    samFile = pysam.Samfile(bam_file,"rb")
    refname = samFile.getrname(0)
    
    for vcf in vcf_data:
        vcf["STRAIN"] = strainID
        irefname = vcf["CHROM"]
        start = int(vcf["POS"])
        ref = vcf["REF"]
        alt = vcf["ALT"]
        quality = vcf["QUAL"]
        end = start + len(ref)
        vcf["END"] = end
        
        alreads = samFile.fetch(irefname,start,end)
        alreads = list(alreads)
        rCount = len(alreads)
        
        vcf["READS"] = rCount
        vcf["READFILE"] = bam_file
        
        if float(quality < minQuality) or rCount < minCount or len(alt) < minSize:
            #print "vcf [%s] quality [%s] count [%s] size [%s] has failed to pass" % (vcf,quality,rCount,len(alt))
            continue
        else:
            vcfDataResult.append (vcf)
        
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
        vReport.add(vID,"bamFile",bam_file)
        #vReport.add(vID,"samFile",samFile)
    
    return (vcfDataResult, vReport)
    
def variantComparison(ref_seq,variant1,variant2):
    r1 = variant1["REF"]
    a1 = variant1["ALT"]
    s1 = int(variant1["POS"])-1
    e1 = s1 + len(r1)
    
    r2 = variant2["REF"]
    a2 = variant2["ALT"]
    s2 = int(variant2["POS"])-1
    e2 = s2 + len(r2)
    
    s = min(s1,s2)
    e = max(e1,e2)
    
    seq1 = ref_seq[s:s1] + a1 + ref_seq[e1:e]
    seq2 = ref_seq[s:s2] + a2 + ref_seq[e2:e]
    
    print "[%s] vs [%s]" % (seq1,seq2)
    

    matrix = matlist.blosum62
    gap_open = -10
    gap_extend = -0.5
    alns = pairwise2.align.globalds(seq1, seq2, matrix, gap_open, gap_extend)
    #score = alns[3]
    max_score = 0
    for align in alns:
        score = align[2]
        if score >= max_score:
            max_score = score
    print max_score
    
    return max_score
    

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
    
    parser.add_option("--nameTag", 
                      dest="nameTag", 
                      default='',
                      metavar="String",
                      help="Name tag to track analysis by")


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
    minSize = 2
    
    #Testing hard coded values
    refTag = "NC_000913_2"

    
    print "ISeuqence Analysis Version 1.2"
    print "Running Mode %s" % (mode)
    
    analysisName = options.nameTag
    
    #File Names
    varianceRepositoryFile = workFolder + "Variant_%s_Repository.bac" % (analysisName)
    mReportFile = workFolder + "Master_%s_Report.bac" % (analysisName)
    mVarianceLog = workFolder + "Master_%s_Var.bac" % (analysisName)
    
    masterReportFile = workFolder + "Master_%s_Report.txt" % (analysisName)
    varianceReportFile = workFolder + "Variance_%s_Report.txt" % (analysisName)
    coverageReportFile = workFolder + "Coverage_%s_Report.txt" % (analysisName)
     
    #========================
    # Construct variant data
    #======================== 
        
    if os.path.exists(mReportFile) and 'rebuild' not in mode:
        print "Loading master repository [%s]" % (mReportFile)
        mReportFh = open(mReportFile,"r")
        masterReport = pickle.load(mReportFh)
        mReportFh.close()
    else:
        print "==> Starting new master report [%s]" % (mReportFile)
        masterReport = DataReport()
        
    if os.path.exists(varianceRepositoryFile) and 'rebuild' not in mode:
        print "Loading variance repository [%s]" % (varianceRepositoryFile)
        fh = open(varianceRepositoryFile,"r")
        varRepository = pickle.load(fh)
        fh.close()
    else:
        print "==> Staring new variant repository"
        varRepository = {}    
        
    '''
    Running in 'build' mode processes sequencing data into 
    index bam files and vcf files.
    '''
    if 'build' in mode:
        
        idCode = parseSampleNameTag(readFile1, readFile2, regex='.*_([ATCG\-]+)_.*',tag="%s")
        sampleCode = "Sample_%s" % idCode
        if verbose: print "Using Sample Code [%s]" % (sampleCode)
        
        currentPath = os.getcwd()
        pathRegex = '.*GEN([0-9].+)'
        match = re.match(pathRegex,currentPath)
        
        if match == None:
            print "No match on regular expression [%s] in  [%s]" %(pathRegex, currentPath)
            pathCode = "000"
        else:
            pathCode = match.group(1)
        if verbose: print "Path Code [%s]" % (pathCode)
        
        #Strain Identifier
        #strainID = "%s_%s" % (sampleCode,pathCode)  
        
        #Variance Tools.
        sTools = SequenceAssemblyTools()
        sTools.verbose = True
        sTools.workingDir = workFolder
        sTools.refGenomeDir = genomeFolder
        sTools.refGenomeTag = refTag
        
        print "Process VCF to Report for sample code [%s]" % (sampleCode)
        
        (bam_file,vcf_file,fasta_file) = sTools.buildAlignment(workFolder, sampleCode, mode)
        
    
    '''
    Process sam and variant calls into unified repository file for processing.
    Consider changing mode to "process"
    '''
    if "build" in mode:
        
        print "Opening bam file [%s]" % (bam_file)
        samFile = pysam.Samfile(bam_file,"rb")
        refname = samFile.getrname(0)
            
        vcf_data = parseVCF(vcf_file)
        
        #----------------------------------------------------
        # Build data record object for variant information
        #----------------------------------------------------
        print "Building Variant information from VCF report"
        strainID = "%s_%s" % (sampleCode,pathCode)     
        (vcf_data_post, vReport) = processVariantData(vcf_data=vcf_data, bam_file=bam_file, strainID = strainID, refName=refname, minQuality=minQuality, minCount=minCount, minSize=minSize)
        
        print "Adding new entries [%s] to master report [%s]" % (len(vReport._data),len(masterReport._data))
        masterReport.extend(vReport)        
        if strainID not in varRepository.keys(): varRepository[strainID] = {}
        varRepository[strainID]["vcfData"] = vcf_data_post
        varRepository[strainID]["sourceFile"] = bam_file
        varRepository[strainID]["seq_file"] = fasta_file

        #Update data file of reported variants
        print "Building Master Report pickle %s" % (mReportFile)
        mReportFh = open(mReportFile,"w")
        pickle.dump(masterReport,mReportFh)
        mReportFh.close()
        
        #Update repository of variants
        print "Building VCF Repository pickle %s" % (varianceRepositoryFile)
        FH = open(varianceRepositoryFile,"w")
        pickle.dump(varRepository,FH)
        FH.close()
        
        if "report" in mode:
            print "Building VCF Report [%s]" % (masterReportFile)
            reportIO = DataReportIO()
            reportIO.writeReport(masterReport, masterReportFile, delimiter='\t')
          
    '''
    Running in "join mode combines variants by regions into matrix of strains and regions.
    '''
    if 'join' in mode:
        
        region_analysis = VariantRegionAnalysis()
            
        print "Joining strain variants"
        varCollection = region_analysis.joinStrainVariants(varRepository, joinDistance=20)
        print"Building collection report"
        (varianceReport, coverageReport) = region_analysis.vcfCollectionReport(varCollection,varRepository)
            
        #Update data file of joined variants
        print "Building Variance Collection %s" % (mVarianceLog)
        varFH = open(mVarianceLog,"w")
        pickle.dump(varCollection,varFH)
        varFH.close()
        
        #Write flat file of VCF master report and joined variant matrix
        if "report" in mode:
            reportIO = DataReportIO()
            print "Writing Variance Regions Report [%s]" % (varianceReportFile)
            reportIO.writeReport(varianceReport, varianceReportFile, delimiter='\t')
            print "Writing Variance Coverage Report [%s]" % (coverageReportFile)
            reportIO.writeReport(coverageReport, coverageReportFile, delimiter='\t')
            
    '''
    Mode for constructing difference matrix and clustering
    !update to only do a pop / or once through of sequence ids.
    '''    
    if 'difference' in mode:
        seqFac = SequenceFactory()
        #config.reflect(configName, seqFac)
        #seqFac.verbose = verbose
        #gRecord = seqFac.getGenBankSequence(genbankID=None, filename=None)
        #seqData = str(gRecord.seq).lower()
        reference_genome = genomeFolder + refTag
        reference_record = SeqIO.parse(open(reference_genome), "fasta")
        reference_record = list(reference_record)[0]
        
        result = {}
        
        
        strain_ids = varRepository.keys()
        strain_ids.sort()
        for strain_id_1 in strain_ids:
            strain_ids.remove(strain_id_1)
            vcf_data_1 = varRepository[strain_id_1]["vcfData"]
            
            for strain_id_2 in strain_ids:
                vcf_data_2 = varRepository[strain_id_2]["vcfData"]
                total_score = 0
                    
                for var in vcf_data_1:         
                    strain_id_check = var["STRAIN"]
                    ref = var["REF"]
                    alt = var["ALT"]
                    quality = var["QUAL"]
                    start = int(var["POS"])-1
                    end = start + len(var["REF"])
          
                    ref_seq = str(reference_record.seq[start:end])
                    print "Var1 (%s,%s): genome [%s] ref[%s] => alt[%s]" % (start,end,ref_seq,ref,alt)
                    
                    var_match = filter(lambda var: start <= int(var["POS"]) <= end or start <= (int(var["POS"]) + len(var["REF"])) <= end, vcf_data_2)
                    
                    for var_target in var_match:
                        print "Var2 (%s,%s): genome [%s] ref[%s] => alt[%s]" % (int(var_target["POS"])-1,int(var_target["POS"])+len(var_target["REF"])-1,ref_seq,var_target["REF"],var_target["ALT"])
                        score = variantComparison(ref_seq, var, var_target)
                        total_score += score
                print "(%s,%s) = %s" % (strain_id_1,strain_id_2,total_score)
                result[(strain_id_1,strain_id_2)] = total_score
                
    print "Done"