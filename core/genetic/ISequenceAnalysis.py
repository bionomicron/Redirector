'''
Created on Aug 11, 2013
Updated Sep 22, 2014

@author: Graham Rockwell
@summary:
    Objects and methods for working analysis of sequencing data.
    Particular tools for analysis of recombination strains.
'''

#Redirector utility methods and objects
from util.Config import ReflectionConfig 
from util.DataReport import DataReport, DataReportIO

#Bio Python Objects
from Bio.Seq import Seq
from Bio import Entrez,SeqIO
from Bio.SeqRecord import SeqRecord

#System tools
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
            
        drefTag = self.refGenomeDir + self.refGenomeTag
        dIdTag = self.workingDir + idTag
        
        #Check Path
        print os.environ["PATH"]
        if self.analysis_tools_path not in os.environ["PATH"]:
            os.environ["PATH"] += self.analysis_tools_path
        
        print "Analysis tools must be in current path [%s]" % os.path.expandvars("$PATH")
        if not os.path.exists("%s.1.bt2" % drefTag):
        #if True:
            print "building reference index"
            call = "bowtie2-build %s %s > bowtie2_ref_index_log.txt" % (drefTag,drefTag)
            if self.verbose: print "executing command: [%s]" % call
            subprocess.call(call,shell=True)
        else:
            if self.verbose: print "building reference index found %s" % (drefTag)
        
        if self.verbose: print ">aligning sequences to reference"
        #try BWA script.
        if os.path.exists("%s.sam" % dIdTag):
            if self.verbose: print "found %s.sam" % dIdTag
        else: 
            print "building %s.sam" % dIdTag
            call = "bowtie2 -q -p 4 -k 3 --local %s -1 %s -2 %s -S %s.sam" % (drefTag,r1,r2,dIdTag)
            if self.verbose: print "executing command: [%s]" % call
            subprocess.call(call,shell=True)
        
        if os.path.exists("%s.bam" % dIdTag):
            print "found %s.bam" % dIdTag
        else:
            print "building %s.bam" % dIdTag
            call = "samtools view -S -b -o %s.bam %s.sam" % (dIdTag,dIdTag)
            if self.verbose: print "executing command: [%s]" % call
            subprocess.call(call,shell=True)
        
        if os.path.exists("%s_s.bam" % dIdTag) and not rebuild:
            print "found %s.bam" % dIdTag
        else:
            print "building %s_s.bam" % dIdTag
            call = "samtools sort %s.bam %s_s" % (dIdTag,dIdTag)
            if self.verbose: print "executing command: [%s]" % call
            subprocess.call(call,shell=True)
        
        if self.verbose: print ">indexing bam file"
        call = "samtools index %s_s.bam" % (dIdTag)
        if self .verbose: print "executing command: [%s]" % call
        subprocess.call(call,shell=True)
        
        resultFile = "%s_s.bam" % (dIdTag)
        
        return (idTag,dIdTag,resultFile)
    
    def vcfProcess(self,bamFile,idTag):
        '''
        @summary:
        Creates variant call format (vcf) from index bam file.
        Currently uses freebayes program, this may not be the best.
        Currently focused on E. coli sequencing data.
        '''
        
        drefTag = self.refGenomeDir + self.refGenomeTag
        resultFile = self.workingDir + "%s.vcf" % idTag
        
        if self.verbose: print "processing vcf file"
        call = "freebayes --fasta-reference %s %s -v %s" % (drefTag,bamFile,resultFile)
        if self.verbose: print "executing command: [%s]" % call
        subprocess.call(call,shell=True)
        
        return resultFile
    
    def buildAlignment(self,workFolder,sampleCode,mode):
        '''
        @summary:
        Pipeline to process paired sequencing files into index bam files and vcf (variant call format) file.
        '''
        
        bamFile = workFolder + "%s_s.bam" % (sampleCode)
        vcfFile = workFolder + "%s.vcf" % (sampleCode)
            
        if not os.path.exists(bamFile) or "rebuild2" in mode:
            print "Building bam file %s" % (bamFile)
            (id,dID,bamFile) = self.samProcessRead(readFile1,readFile2,idTag=sampleCode)
        else:
            print "%s found" % (bamFile)
                            
        if not os.path.exists(vcfFile) or "rebuild2" in mode:
        #if True:
            print "Building vcf file %s" % (vcfFile)    
            vcfFile = sTools.vcfProcess(bamFile, idTag = sampleCode)
        else:
            print "%s found" % (vcfFile)
            
        return(bamFile,vcfFile)

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
               
def processVariantData(vcfData, bamFileName, strainID = '', refName='', minQuality = 30, minCount = 1,minSize=0):
    '''
    @return [vcf]
    @summary:
    Process variant data link combine with related sequencing reference.
    Also creates a report of all variant data.
    '''
    #vReport = Report() #Testing new report object
    vReport = DataReport()
    vcfDataResult = []
    
    print "Opening bam file [%s]" % (bamFileName)
    samFile = pysam.Samfile(bamFile,"rb")
    refname = samFile.getrname(0)
    
    for vcf in vcfData:
        vcf["STRAIN"] = strainID
        irefname = vcf["CHROM"]
        start = int(vcf["POS"])
        ref = vcf["REF"]
        alt = vcf["ALT"]
        quality = vcf["QUAL"]
        end = start + len(alt)
        vcf["END"] = end
        
        alreads = samFile.fetch(irefname,start,end)
        alreads = list(alreads)
        rCount = len(alreads)
        
        vcf["READS"] = rCount
        vcf["READFILE"] = bamFileName
        
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
        vReport.add(vID,"bamFile",bamFile)
        #vReport.add(vID,"samFile",samFile)
    
    return (vcfDataResult, vReport)
    

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
        print "=====>Starting new master report"
        masterReport = DataReport()
        
    if os.path.exists(varianceRepositoryFile) and 'rebuild' not in mode:
        print "Loading variance repository [%s]" % (varianceRepositoryFile)
        fh = open(varianceRepositoryFile,"r")
        varRepository = pickle.load(fh)
        fh.close()
    else:
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
        sTools.refGenomeTag = refTag
        sTools.workingDir = workFolder
        sTools.refGenomeDir = genomeFolder
        sTools.refGenomeTag = refTag
        
        print "Process VCF to Report for sample code [%s]" % (sampleCode)
        
        (bam_file,vcf_file) = sTools.buildAlignment(workFolder, sampleCode, mode)
        
        print "Opening bam file [%s]" % (bam_file)
        samFile = pysam.Samfile(bam_file,"rb")
        refname = samFile.getrname(0)
            
        vcf_data = parseVCF(vcf_file)
        
        #----------------------------------------------------
        # Build data record object for variant information
        #----------------------------------------------------
        print "Building Variant information from VCF report"
        strainID = "%s_%s" % (sampleCode,pathCode)     
        (vcf_data_post, vReport) = processVariantData(vcfData=vcf_data, bamFileName=bam_file, strainID = strainID, refName=refname, minQuality=minQuality, minCount=minCount, minSize=minSize)
        print "Adding new entries [%s] to master report [%s]" % (len(vReport._data),len(masterReport._data))
        masterReport.extend(vReport)
        
        #Fill repository
        if strainID not in varRepository.keys():
            varRepository[strainID] = {}
            
        varRepository[strainID]["vcfData"] = vcf_data_post
        varRepository[strainID]["sourceFile"] = bam_file

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
    join mode combines variants by regions into matrix of strains and regions.
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
        
    print "Done"