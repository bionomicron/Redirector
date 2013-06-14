'''
@author: Graham Rockwell
@change: 2013 06 06
@summary: 
Sequence Tool Box
Tool kit for 
1) Generating recombination oligos
2) Finding sequence melting temperatures
3) Blasting sequences against target genomes or other sequences

Many functions based in part or whole, with individual revisions, 
on work from Sriram's further based on Billy's code
'''

import os, re, math, subprocess
import time, copy

from Bio import Entrez,SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature,FeatureLocation
from Bio.Blast import NCBIStandalone, NCBIXML
from Bio.SeqRecord import SeqRecord
from Bio.Blast.Applications import BlastallCommandline

from core.util.Report import Report

class SequenceTools:
    '''
    Group of tools for finding sequence properties and generally usefull information
    
    Functions:
    * Calculate secondary structure energy for sequence
    * Blast sequence and return record list
    * Find oligoTM
    * Optimize a group of oligos for a target TM
    * Generate a list of primers optimzied for a TM, optimized for boundary distance.
    
    Dependencies
    * blast
    * hybrid-ss-min (used to be mfold)
    
    '''
    
    def __init__(self):
        self.verbose = False
        
    def calcSecondaryStructure(self, seq):
        '''
        Caculate secondary structure using
        Uses what used to be mfold.
        @param seq: input sequence
        @type seq: string or seq object
        @rtype: float
        '''
        
        cmd = "hybrid-ss-min --NA=DNA -q " + str(seq)
        p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
        p.wait()
        dg = p.stdout.read()
        ans = float(dg)
        return ans
    
    def seqBlast(self, blastDB, blastExe, seqFile, blastType = "blastn", scoreMin = 1e-3, logFile = None):
        '''
        command line blast
        blastall -d database -i query -p blastn -o blastout
        '''
        
        if not os.path.exists(os.path.expanduser(seqFile)):
            print "(ignore) %s file not found" %(seqFile)
        
        if not os.path.exists(os.path.expanduser(blastDB+".nsq")):
            print "(ignore) %s file not found" % (blastDB)
            
        (resultHandle,errorHandle) = NCBIStandalone.blastall(blastExe,blastType,blastDB,seqFile)       
        time.sleep(5)
        blastRecords = NCBIXML.parse(resultHandle)

        blastRecords = list(blastRecords)
        resultHandle.close()
        errorHandle.close()

        return blastRecords
    
    def blastSelection(self, blastRecords, start=float("-inf"), end=float("inf"), include=True, scoreMin=1e-3):
        '''
         Convert blast records from within region into list of alignments
        '''
        
        result = []
        iBlastRecords = copy.deepcopy(blastRecords)
        for r in iBlastRecords:
            ialignments = []
            for alignment in r.alignments:
                ihsps = []
                
                for hsp in alignment.hsps:

                    if hsp.expect < scoreMin:
                        
                        istart = hsp.sbjct_start
                        iend = hsp.sbjct_end
                        
                        if start<istart and iend<end:
                            ihsps.append(hsp)
                        elif include and (istart < end and start < iend):
                            ihsps.append(hsp)
                            
                if len(ihsps) > 0:
                    alignment.hsps=ihsps
                    ialignments.append(alignment)
                    
            r.alignments = ialignments
            result.append(r)
            
        return result
            
    def seqBlastToFeatures(self, blastDB, blastExe, seqFile, blastType = "blastn",scoreMin = 1e-3, logFile = None):
        '''
        Blast sequence file against blast database 
        parse files into records.
        This function may not work as well on very large blast comparisons because 
        it does a full read of the result for the conversion to features.
        '''
        blastRecords = self.seqBlast( blastDB, blastExe, seqFile, blastType = "blastn", scoreMin = 1e-3, logFile = None)
        
        result = []
        index = 0
        for r in blastRecords:
            recordFeatures = []
            for alignment in r.alignments:
                name = alignment.title
                query = r.query
                for hsp in alignment.hsps:
                    if hsp.expect < scoreMin:
                        (ts,ss) = hsp.frame
                        strand = ss
                        start = hsp.sbjct_start
                        end = hsp.sbjct_end
                        location = FeatureLocation(start,end)
                        feature = SeqFeature(id=query,location=location,strand=strand)
                        aMatch = hsp.query + "\n" + hsp.match + "\n" + hsp.sbjct
                        feature.qualifiers["query"] = hsp.query
                        feature.qualifiers["subject"] = hsp.sbjct
                        feature.qualifiers["alignment"] = aMatch
                        recordFeatures.append(feature)
            result.append(recordFeatures)
            index = index + 1
            
        return result

    def oligoTm(self, seqobj, C_primer= 50, C_Mg= 1.5, C_MonovalentIon= 50, C_dNTP=0.4):
        """
        Takes either a SeqRecord object, a Seq object, or a string
        and computes the melting temperature based on the NN model (yes?).
        This is based largely on Kun's code
        
        CHECK THE NN PARAMETERS
        From Uri Laserson
        
        The defaults are from Kun's caculator settings online.
        TMs returned are slightly different than other caculators on line,
         this may be a result of molecular levels that are not accounted for in those
         calculators.
        
        @type seqobj: String, Seq object, or SeqRecord
        @var seqObj:  The sequence for which the TM is being produced
        @type C_primer: float
        @type C_mg: float
        @type C_MonovalentIon: float
        @rtype: float
        @return: Returns the melting temperature (TM) for the input sequence 
        """
        
        if isinstance(seqobj,SeqRecord):
            seq = seqobj.seq.tostring().upper()
        elif isinstance(seqobj,Seq):
            seq = seqobj.tostring().upper()
        elif isinstance(seqobj,str):
            seq = seqobj.upper()
        
        # set the default tm parameters
        C_primer = 50 # 10nM intracellular concentration (1 molecule in 1e-15L
        C_Mg = 1.5 #10mM intracellular concentration
        C_MonovalentIon = 50  
        C_dNTP = 0.8
        percentage_DMSO = 0
        percentage_annealed = 50
        
        percentage_annealed = percentage_annealed/100.0
        percentage_DMSO = percentage_DMSO/100.0
        #Some constants
        R = 1.987
        deltaH =  {"AA": -7.9,  "TT": -7.9, "AT": float((-1)*7.2), "TA": -7.2, "CA": -8.5, "TG": -8.5, "GT": -8.4, "AC": -8.4,"CT": -7.8, "AG": -7.8, "GA": -8.2, "TC": -8.2,"CG": -10.6,"GC": -9.8, "GG": -8.0, "CC": -8.0, "A": 2.3, "T": 2.3, "G": 0.1, "C": 0.1}
        deltaS = { "AA": -22.2, "TT": -22.2, "AT": -20.4, "TA": -21.3, "CA": -22.7, "TG": -22.7, "GT": -22.4, "AC": -22.4, "CT": -21.0, "AG": -21.0, "GA": -22.0, "TC": -22.0,"CG": -27.2, "GC": -24.4, "GG": -19.9, "CC":-19.9, "A": 4.1, "T": 4.1, "G": -2.8, "C": -2.8}
        
        C_SodiumEquivalent = C_MonovalentIon + 120 * math.sqrt(C_Mg-C_dNTP)
        seqLength = len(seq)
        dH = deltaH[str(seq[0])] + deltaH[str(seq[len(seq)-1])]
        dS = deltaS[seq[0]] + deltaS[seq[len(seq)-1]]
        for     i in range(0, seqLength - 1):
            dH = dH + deltaH[str(seq[i:i+2])]
            dS = dS +  deltaS[seq[i:i+2]]
        dS = dS + 0.368 * seqLength * math.log(C_SodiumEquivalent/1000.0)
        val = math.log(C_primer*(1-percentage_annealed)/percentage_annealed)
        Tm =(dH * 1000) / (dS + R * (math.log(C_primer*(1-percentage_annealed)/percentage_annealed)-21.4164)) - 273.15 - 0.75*percentage_DMSO
        return Tm
    
    def oligolistTm(self, sequences):
        '''
        finds TMs for list of sequences
        @param sequences: list of sequences for which to calculate TMs
        @type sequences: list of strings or sequence objects.
        '''
        Tm = []
        for i in sequences:
            Tm.append(oligoTm(i))
        return Tm
    
    def scanOligoTm(self,seq,start,size,searchSize,targetTM,delta = 0.1):
        '''
        Search though possible oligo locations around targte area for oligos with desired TM
        preference is given to tm similarity then to close ness to start sight
        with a margin of TM difference of delta, default 0.1
        '''
        result = {}
        minTmDiff = 500
        minAdjust= 500
        rTm = None
        rAdjust = None
        rSeq = None
        
        for adjust in range(-searchSize,searchSize):
            iStart = start+adjust
            iEnd = start+size+adjust
            iSeq = seq[iStart:iEnd]
            iTM = self.oligoTm(iSeq)
            tmDiff = abs(iTM - targetTM)
            if tmDiff < minTmDiff or (abs(adjust) < abs(minAdjust) and abs(tmDiff - minTmDiff) < delta):
                minTmDiff = tmDiff
                minAdjust = adjust
                rTm = iTM
                rAdjust = adjust
                rSeq = iSeq
        result = (rTm,rAdjust,rSeq)
        
        return result
    
    def optimizeOligosForTm(self,seqMap,start,size,searchSize, targetTM):
        '''
        optimize a list of oligos for a target TM
        '''
        seqTms = {}
        for k in seqMap.keys:
            s = seqMap[k]
            v = self.scanOligoTm(s,start,size,searchSize,targetTM)
            seqTms[k] = v
            
        return seqTms
    
    def findPrimers(self,seq,targetMap,boundary,oligoSize,searchSize,targetTm):
        '''
        @input targetMap: a list of genomic locations with sequence names
        
        find a list of sequencing primers for a list of target locations
        
        usage:
        set boundary and oligo size
        read in report, list of locations
        run find sequence primers.
        write report of result
        '''
        result = Report()
 
        for k in targetMap.keys():
            targetLocation = targetMap[k]
            
            upLocation = targetLocation - boundary
            upStart = targetLocation - boundary - searchSize - oligoSize
            upEnd = targetLocation - boundary + searchSize
            upSeq = seq[upStart:upEnd]
            
            downLocation = targetLocation + boundary
            downStart = targetLocation + boundary - searchSize
            downEnd = targetLocation + boundary + searchSize + oligoSize
            downSeq = seq[downStart:downEnd]
            downSeq = downSeq.reverse_complement()
            
            start = searchSize
            (uTm,uAdjust,oUpSeq) = self.scanOligoTm(upSeq,start,oligoSize,searchSize,targetTm)
            (dTm,dAdjust,oDownSeq) = self.scanOligoTm(downSeq,start,oligoSize,searchSize,targetTm)
            dAdjust = -dAdjust
            
            ucLocation = upLocation + uAdjust
            dcLocation = downLocation + dAdjust
            
            result.add(k,"sequencing location",targetLocation)
            
            result.add(k,"foward primer",oUpSeq)
            result.add(k,"foward adjust",uAdjust)
            result.add(k,"foward location",ucLocation)
            result.add(k,"foward TM",uTm)
            
            result.add(k,"reverse primer",oDownSeq)
            result.add(k,"reverse adjust",dAdjust)
            result.add(k,"reverse location",dcLocation + dAdjust)
            result.add(k,"reverse TM",dTm)
            
        return result        

class FeatureTools:
    '''
    Tools set for discovering relationships between features.
    Finds overlapping sequences 
    '''
    
    def __init__(self):
        pass
    
    def testFeatureOrder(self,featureList):
        '''
        Finds features that are overlapping and produces a list of updated positions
        Other algorithms assume feature list is ordered.  I think this is true for E. coli, though I don't know now that there are ncRNAs.  
        Also does not check if there are double overlaps between genes... though don't think so for E. coli
        '''
        startLocations=[]
        #skips first gene, assumes no overlap from other side of circular dna... OK for E. coli for now
        for i in range(0,len(featureList)):
            #find upstream homology regions
            startLocations.append(featureList[i].location.start.position)
        for i in range(1,len(startLocations)):
            if startLocations[i-1]>startLocations[i]:
                print "Found a feature at position " +str(startLocatios[i-1])+" in list that is not completely in order"
    
    def findFeatureOverlaps(self,featureList):
        '''
        Finds features that are overlapping and produces a list of updated positions
        Assumes genbank feature list is ordered, and no more than a single overlap.  I think this is true for E. coli, though I don't know now that there are ncRNAs.
        '''
        startLocations=[]
        endLocations=[]
        updatedStartLocations=[]
        updatedEndLocations=[]
        overlapStart=[]
        overlapEnd=[]
        #skips overlap between first and last gene, assumes no overlap from other side of circular dna... OK for now with E. coli
        for i in range(0,len(featureList)):
            #find upstream homology regions
            startLocations.append(featureList[i].location.start.position)
            endLocations.append(featureList[i].location.end.position)
            updatedStartLocations.append(featureList[i].location.start.position)
            updatedEndLocations.append(featureList[i].location.end.position)
    
        for i in range(0,len(startLocations)): 
            overlapStart.append(False)
            overlapEnd.append(False)
            if not i==0:
                basesBetweenGenesStart = startLocations[i]-endLocations[i-1]
            else: 
                basesBetweenGenesStart = 1 # arbitrary to pass screening, just needs to be positive
            
            if not i==(len(startLocations)-1):            
                basesBetweenGenesEnd = startLocations[i+1]-endLocations[i]
            else:
                basesBetweenGenesEnd = 1 # arbitrary to pass screening, just needs to be positive
            
            #update start locations of gene i
            if str(featureList[i].type) == "CDS": # maintain reading frame of CDS
                j=3 # remove start codon
                if basesBetweenGenesStart<-2: # correct for overlap maintaining reading frame of CDS
                    j = ((basesBetweenGenesStart * -1)//3 + 1) * 3
                    overlapStart[i]=True
                updatedStartLocations[i] = startLocations[i] + j
            elif basesBetweenGenesStart<1 and str(featureList[i].type) == "ncRNA" : # no need to maintain ORF for ncRNA
                updatedStartLocations[i] = endLocations[i-1] + 1
                overlapStart[i]=True
    
            #updated end locations of gene i
            if str(featureList[i].type) == "CDS": # maintain reading frame of CDS
                j = 3 # reomve stop codon
                if basesBetweenGenesEnd<-2: # correct for overlap maintaining reading frame of CDS
                    j = ((basesBetweenGenesEnd * -1)//3 + 1) * 3
                    overlapEnd[i]=True
                updatedEndLocations[i] = endLocations[i] - j
            elif basesBetweenGenesStart<1 and str(featureList[i].type) == "ncRNA" : # no need to maintain ORF for ncRNA
                updatedEndLocations[i] = endLocations[i+1] - 1
                overlapEnd[i]=True
      
        return updatedStartLocations, updatedEndLocations, overlapStart, overlapEnd
    
    def testFeatureOverlaps(self,featureList):
        '''
        Debugging method to find features that are overlapping 
        Assumes feature list is ordered, and no more than a single overlap.  I think this is true for E. coli, though I don't know now that there are ncRNAs.
        '''
        startLocations=[]
        endLocations=[]
        #skips first gene, assumes no overlap from other side of circular dna... OK for E. coli for now
        for i in range(0,len(featureList)):
            #find upstream homology regions
            startLocations.append(featureList[i].location.start.position)
            endLocations.append(featureList[i].location.end.position)
        for i in range(1,len(startLocations)):
            basesBetweenGenes = startLocations[i]-endLocations[i-1]
            if basesBetweenGenes<0:
                print "Gene " + str(i)+" is overlapping by " + str(-1*basesBetweenGenes) +" bp."
            # checks if there are very large overlaps and prints out some diagnostics
            if basesBetweenGenes<-500:
                print str(i-1) + "weirdo"
                print featureList[i-1]
                print str(i) + "weirdo"
                print featureList[i]
                
    
class SequenceFactory:
    '''
    Factory 
    Generates genome sequence files.
    Pulls down genome sequence files from database
    or serves local genome sequence file into Seq object.
    '''
    
    def __init__(self,dataDirectory=''):
        self.dataDirectory = dataDirectory
        self.genbankID = ''
        self.email = "automated"
        self.verbose = False
        
    def _genbankFileName(self,genbankID):
        return self.dataDirectory + "gi_" + genbankID + ".gbk"
    
    def downloadGenBank(self,genbankID=None,filename = None):
        '''
        Downloads GenBank file given an GenBank ID and an email address to give NCBI
        Defaults to GenBank record of E. coli K12 MG1655
        '''
        if genbankID == None:
            genbankID = self.genbankID
        if filename == None:
            filename = self.dataDirectory + "gi_" + genbankID + ".gbk"
        
        if not os.path.isfile(filename):
            if self.verbose: print "Downloading..."
            net_handle = Entrez.efetch(db="nucleotide", rettype="genebank", id=genbankID)
            out_handle = open(filename, "w")
            out_handle.write(net_handle.read())
            out_handle.close()
            net_handle.close()
            if self.verbose: print "Saved"
            return filename
        else:
            if self.verbose: print "Getting GenBank record from local disk"
            return filename
    

    def parseGenBank(self, filename = None):
        '''
        parse genbank file into SeqRecord object
        default is E. Coli genome
        '''
        
        if filename == None:
            filename = self.dataDirectory + "gi_" + self.genbankID + ".gbk"
            
        if self.verbose: print "Parsing Whole Genome..."
        record = SeqIO.read(open(filename,"rU"), "genbank")
        if self.verbose: print "Parsing Genome Complete"
        return record
    
    def getGenBankSequence(self,genbankID=None,filename=None):
        genbankFile = self.downloadGenBank(genbankID, filename)
        if not os.path.exists(genbankFile):
            time.sleep(5)
        genbankSequence = self.parseGenBank(genbankFile)
        return genbankSequence
    

    def printFeature(self,feature):
        '''
        debugging method to print a single feature in a smaller format
        '''
        out = "id: %s\t" % feature.id
        out += "type: %s\t" % feature.type 
        out += "location: %s\t" % feature.location 
        out += "strand: %s" % feature.strand
        print out
     

    def printFeatures(self,record):
        '''
        debugging method to print all features in a seqrecord in a smaller format
        '''
        print "Parsing Individual Genes..."
        for feature in record.features:
            self.printFeature(feature)
        print "Parsing Genes Complete"
    
    def extractCDSncRNAFeatures(self,record):
        '''
        extracts the CDS and ncRNA features and returns them in a list format; does not include pseudogenes
        '''        
        cdsncRNAs = []
        if self.verbose: print "Extracting CDS and ncRNA features"
        for feature in record.features:
            if str(feature.type) == "CDS" or str(feature.type) == "ncRNA":
                #Test to make sure not a pseudo-gene
                if not feature.qualifiers.has_key("pseudo"):
                    cdsncRNAs.append(feature)
        
        return cdsncRNAs
    
    def _parseConversion(self, value, conversionMap, default = 0):
        '''
        utility function to parse values from a report
        '''
        if value not in conversionMap.keys():
            result = default
        else:
            result = int(conversionMap[value])
        return result
                
                
    def convertReportToFeatures(self,report,conversionMap,startTag,endTag,strandTag,qualTag="gene"):
        '''
        parse 2d hash into feature list
        '''
        result = []
        for rowName in report.returnRowNames():
            start = float(report.getElement(rowName,startTag))
            end = float(report.getElement(rowName,endTag))
            location = FeatureLocation(start,end)
            strand = self._parseConversion(report.getElement(rowName,strandTag),conversionMap)
            feature = SeqFeature(id = rowName,location=location,strand=strand)
            feature.qualifiers[qualTag] = [rowName]
            for colName in [rowName,report.returnColumnNames()]:
                if colName in conversionMap.keys():
                    newName = conversionMap[colName]
                    value = report[colName]
                    feature.qualifiers[newName] = [value]
            result.append(feature)
        return result  
        
        
class ControlRegionTools:
    '''
    Tools for performing recombinations on control regions (Promoter or RBS)
    Generates control region targets for recombination
    Usage:
     Input list of genes and get list of control regions for each gene with upper case motif locations.
    '''
    
    def __init__(self):
        self.verbose = False
    
    def getSequenceRegion(self,seq,location,size=1,boundary=0,dir=1.0,useUpper=True,replace='',end=None):
        '''
        Produce a sequence from the target location.
        Put indicated region in uppercase if so noted.
        Shows control regions and surrounding regions.
        Used as starting point for generating control region diversity oligos.
        @Location: start of target region
        @end: end of target region, used if replacement size != target region
        '''        
        seq = Seq(seq)
        
        if replace != '':
            size = len(replace)
        
        #find start and end    
        if end == None:
            if dir == 1:
                start = int(location)
                end = int(location - size)
            if dir == -1:
                start = int(location - size)
                end = int(location)
        else:
            start = int(location)
            end = int(end)
            size = int(end-start)
        
        #set boundary locations
        lowBound = int(start-boundary)
        highBound = int(end+boundary)
        
        #get region 
        t = seq[lowBound:highBound]
            
        if start < lowBound:
            print "sequence retreval error, start %s less than lower bound %s" % (start,lowBound)
        if end > highBound:
            print "sequence retreval error, end %s less than higher bound %s" % (end,highBound)
        
        a = seq[lowBound:start]
        b = seq[start:end]
        if replace != '':
            b = replace
        c = seq[end:highBound]

        if dir == 1:
            if useUpper: b = str(b).upper()
            result = a + b + c
        if dir == -1:
            b = b.reverse_complement()
            if useUpper: b = str(b).upper()
            result = c.reverse_complement() + b + a.reverse_complement() 
        
        return result
    
    def generateReplacementFeatures(self,features,sequence):
        pass
        #currently a place holder.
    
    def localFeatures(self,features,targetRecord,range = 1000, closestOnly = False):
        '''
        Find the local features to a target feature and create hash out of them.
        '''
        result = {}
        targetFeatures = targetRecord.features
        
        for feature in features:
            subTRecord = None
            name = feature.qualifiers["gene"][0]
            start = feature.location.start.position
            end = feature.location.end.position
            strand = feature.strand
            if strand == 1:
                subTRecord = targetRecord[start-range:start]
                mod = start-range
            if strand == -1:
                subTRecord = targetRecord[end:end + range]
                mod = end
            localFeatures = []
            for tFeature in subTRecord.features:
                tName = tFeature.qualifiers["gene"][0]
                tStart = tFeature.location.start.position + mod
                tEnd = tFeature.location.start.position + mod
                tFeature.location.start.position = tStart
                tFeature.location.end.position = tEnd

                if strand == 1: arrow = "->"
                if strand == -1: arrow = "<-"

                if tFeature.strand == strand:
                    #if verbose: print "%s %s %s [%s] %s %s %s " % (str(start),arrow, str(end),str(strand),str(tStart),arrow, str(tEnd))
                    localFeatures.append(tFeature)
            result[name] = localFeatures
        
        return result
        
    
    def getControlReport(self,genes,targetRecord,sequence,boundary=0,range=1000):
        '''
        Create report of control regions for listed genes
        '''
        result = Report()
        seqProp = RecombinationOligoFactory()
        
        sdata = str(sequence).lower() 
        promoters = {}
        locations = {}
        arrow= ["<-","-","->"]
        
        if self.verbose: print "finding local features"
        
        localFeatures = self.localFeatures(genes,targetRecord,range=range)
    
        if self.verbose: print "features found, creating report"
        
        for feature in genes:
            name = feature.qualifiers["gene"][0]
            locTag = feature.qualifiers["locus_tag"][0]
            start = feature.location.start.position
            end = feature.location.end.position
            strand = feature.strand
            
            result.add(name,"locus_tag",locTag)                                                     
            result.add(name,"gene_start",start)                                                     
            result.add(name,"gene_end",end)                                                     
            result.add(name,"gene_strand",strand)                                                     
            
            if strand == 1:
                loc = start
            if strand == -1:
                loc = end
                
            rbsSeq = self.getSequenceRegion(sdata, loc, 3, boundary, strand)
            oligoStrandRbs = seqProp.strandChooser(feature)
                
            result.add(name,"rbs_start",loc)                                                     
            result.add(name,"rbs_region",rbsSeq)
            result.add(name,"rbs_Oligo_Strand",str(oligoStrandRbs))
            
            iLocalFeatures = localFeatures[name]
                                                       
            count = 0
            
            for tFeature in iLocalFeatures:
                count = count + 1
                tName = tFeature.qualifiers["gene"][0]
                tStart = tFeature.location.start.position
                tEnd = tFeature.location.end.position
                tStrand = tFeature.strand
                oligoStrand = seqProp.strandChooser(tFeature)
                iArrow = arrow[tStrand + 1]
                if tStrand == strand:
                    tSize = -1*tStrand
                    pSeq = self.getSequenceRegion(sdata, tStart, tSize, boundary, strand)
                tag = "%s:[%s %s %s] t[%s] = %s" % (tName,tStart,iArrow,tEnd,oligoStrand,pSeq)
                colName = "promoter_%s" % (count)
                result.add(name,colName,tag)
            
        return result
    
        
class RecombinationOligoFactory:
    '''
    Factory for generating oligos for recombination
    '''
    
    def __init__(self):
        '''
        Factory for recombation oligos
        defaults based on ecoli genome
        '''
        self.sequenceTools = SequenceTools()
        self.verbose = False
        self.stars = 4
        
        #defaults for replichore
        #replichore start region, middle of oriC region in E. coli MG1655  (3924035+3924478)/2
        self.repstart = 3924256
        #replicore end region, middle of E. coli dif sites  (1588774+1588801)/2
        self.repend = 1588787
        
        #defaults for recombination olgio target discovery
        self.oligoSearchSize = 45
        self.structureScoreCutOff = -12.5
        
        #defaults for primer optimization
        self.primerDistance = 250
        self.primerTm = 55
        self.primerSize = 20
        self.primerSearchSize = 60
    
    def strandChooser(self, thefeature):
        '''
        Choosed target strand for oligo recombination,
        based on origin of replication and terminus locations. 
        Returns which strand to target for any particular sequence feature 
        Default settings are for MG1655
        
        Fuction:
        If target is in Replichore 1: 
            Replication is continious on + strand (genomic: start 5'-> end 3')
            Target should be complementry to + strand 
            There for target should be on - strand (genomic start 3' <- end 5')
        If target is on Replichore 2: 
            want to target - strand, so return + strand because it is complementary continious replication strand
        '''

        repstart = self.repstart
        repend = self.repend
        
        featurestart = thefeature.location.start.position
        featureend = thefeature.location.end.position
        featuremiddle = (featurestart+featureend)/2
        
        
        if repstart<featuremiddle or featuremiddle<repend:
            #Replicore 1
            return -1
        else:
            #Replicore 2
            return 1
        
    
    def retrievePrimerHRregions(self, therecord, thefeature, updatedstart, updatedend, hrlength, numofprimers):
        '''
        unsure function
        '''
        updatedstart = updatedstart-hrlength
        updatedend = updatedend+hrlength
        myseq = therecord.seq[updatedstart:updatedend]
        strand = self.strandChooser(thefeature)
        #reverse complement if necessary
        if strand==1:
            #reverse complement
            myseq = myseq.reverse_complement()
        #construct possible forward subsequences:
        forwardhr = []
        reversehr = []
        for i in range(numofprimers):
            start = i*3
            end = i*3+hrlength
            if strand==1:
                reversehr.append(myseq[start:end])
            else:
                forwardhr.append(myseq[start:end])
        #construct possible reverse subsequences:
        seqlength = len(myseq)
        for i in range(numofprimers):
            start = (seqlength-hrlength-i*3)
            end = (seqlength-(i*3))
            if strand==1:
                forwardhr.append(myseq[start:end])
            else:
                reversehr.append(myseq[start:end])
        #return the two lists
        
        return forwardhr, reversehr
    
    def optimizeSecondaryStructure(self,seq,start,end,searchSize,cutOff,ratio):
        '''
        Find oligo with best secondary structure in particular location.
        ratio determines the weigh of the distance to the prefered location vs the 2ndary structure energy
        '''
        seqStack = []
        seqBest = []
        for adjust in range(-searchSize,searchSize):
            iStart = start + adjust
            iEnd = end + adjust
            subSeq = seq[iStart:iEnd]
            foldScore = self.sequenceTools.calcSecondaryStructure(subSeq)
            (low,high) = cutOff
            if high != None and foldScore > high:
                score= abs(adjust)
                seqStack.append((score,abs(adjust),foldScore,subSeq))
            elif low == None or foldScore > low:
                score = abs(adjust) + abs(foldScore)*ratio
                seqStack.append((score,abs(adjust),foldScore,subSeq))
        return seqStack
    
    def findSequenceLocation(self, record, seq):
        '''
        find the location of a given sequence in a the sequence of a record.
        '''
        s =record.seq
        if seq in s:
            strand = 1
            start = s.index(seq)
            end = start + len(seq)
            return SeqFeature(start = start, end = end, strand = strand)
        rc = record.seq.reverse_complement()
        if seq in rc:
            strand = -1
            start  = rc.index(seq)
            end = start + len(seq)
            return SeqFeature(start = start, end = end, strand = strand)
        return None
    
    def addStars(self,seq,number):
        '''
        add * character to sequence
        indicating protection to 5' end of recombination oligos
        '''
        result = ""
        
        s1 = 0
        s2 = number
        s3 = len(seq)-number
        s4 = len(seq)
        
        for i in range(s1,s2):
            result += seq[i] + "*"
            
        result += seq[s2:s3]
        
        for i in range(s3,s4):
            result += "" + seq[i]
            
        return result
    
    def optimizeOligo(self,start,end,sequence,boundary,searchSize,cutOff):
        '''
        finds optimal subsection of sequence for recombination oligo
        '''
        s = str(r.seq)
        sx = sequence[start:end]
        

        lowSearchBound = start - searchSize - 5
        highSearchBound = end + searchSize + 5
            
        if lowSearchBound < 0:
            lowSearchBoud = 0
        if highSearchBound > len(s):
            highSearchBoud = len(s)
                
        if genomicStrand == laggingComplementStrand:            
            s = r.seq[lowSearchBound:highSearchBound]
        else:
            s = r.seq[lowSearchBound:highSearchBound].reverse_complement()
    
        searchEnd = (end-start+searchSize+5)
        searchStart = searchSize + 5
            
        try:
            testSeqs = self.optimizeSecondaryStructure(s, searchStart, searchEnd, searchSize, cutOff, ratio = 2)
            testSeqs.sort()
            if len(testSeqs) == 0:
                (adjust,foldScore,bestSeq) = ("na","na","na")
            else:
                (score,adjust,foldScore,bestSeq) = testSeqs.pop(0)
        except:
            testSeqs = []
            (adjust,foldScore,bestSeq) = ("na","na",sx)
        
        return (adjust,foldScore,bestSeq)
    
    def parseAlignments(self,records,featureLocations):
        result = Report()
        index = 0
        
        logFile = open("oligoLog.txt","w")
        targetMap = {}
        
        for r in records:
            id = r.id

            features = featureLocations[index]
            hits =  len(features)
            genomicStart = features[0].location.start.position
            genomicEnd = features[0].location.end.position
            genomicStrand = features[0].strand
            laggingComplementStrand = self.strandChooser(features[0])
            
            if "alignment" in features[0].qualifiers.keys():
                aMatch = features[0].qualifiers["alignment"]
            else: 
                aMatch = ''
            
            targetMap[id] = (genomicStart + genomicEnd)/2
            
            logFile.write(id+"\n")
            logFile.write(aMatch+"\n")
            #if self.verbose: print aMatch
            
            s = str(r.seq)
            originalString = s
            
            result.add(id,"original", originalString) 
            result.add(id,"hits", hits) 
            result.add(id,"genomic_start", genomicStart) 
            result.add(id,"genomic_end", genomicEnd) 
            result.add(id,"genomic_strand", genomicStrand)
            result.add(id,"match", aMatch)
        
        return result
        
        
    def generateTargetingOligos(self, records, featureLocations, tagRE, boundary, searchSize, cutOff):
        '''
        Generate a list of oligos for recombination in target locations
        and return report with the targets and sequencing oligos
        
        All oligos printed 5' -> 3'
        Control upstream will be to the left if strands are preserved and to the right 
        if strands are switched when matching lagging complement
        
        oligos are selected discovered as an optimized subsection of the presented sequences 
        
        @records: sequences from which to select oligos
        @featureLocations: a list of locations that place the features in a genome
        @tagRE: regular expression for finding taged sequence with in feature sequences.
        @bounary: oligo flanking region size
        @searchSize: distance in base pairs to search for optimal oligo
        @cutOff: limit of viable fold change energy for chose oligos.
        '''
        
        result = Report()
        index = 0
        
        logFile = open("oligoLog.txt","w")
        targetMap = {}
        sRegions = []
        
        for r in records:
            id = r.id

            features = featureLocations[index]
            hits =  len(features)
            genomicStart = features[0].location.start.position
            genomicEnd = features[0].location.end.position
            genomicStrand = features[0].strand
            laggingComplementStrand = self.strandChooser(features[0])
            
            if "alignment" in features[0].qualifiers.keys():
                aMatch = features[0].qualifiers["alignment"]
            else: 
                aMatch = ''
            
            targetMap[id] = (genomicStart + genomicEnd)/2
            
            logFile.write(id+"\n")
            logFile.write(aMatch+"\n")
            #if self.verbose: print aMatch
            
            s = str(r.seq)
            originalString = s
            
            #Find larget section using special targeting tag
            matchTag = re.search(tagRE,s)
            if matchTag == None:
                tagLoc = len(s)/2
            else:
                targetTag = matchTag.group(0)
                tagLoc = s.index(targetTag) + len(targetTag)/2
                
            if self.verbose: print "Target [%s] location [%s]" % (targetTag,tagLoc)

            start = int(tagLoc - boundary)
            end = int(tagLoc + boundary)
            
            #!May need a little touch up
            if start < 0:
                start = 0
                end = int(boundary*2)
            if end > len(s):
                end = len(s)
            
            if self.verbose: print"region %s -> %s of %s" % (start,end,len(s))

            lowSearchBound = start - searchSize - 5
            highSearchBound = end + searchSize + 5
            
            if lowSearchBound < 0:
                lowSearchBound = 0
            if highSearchBound > len(s):
                highSearchBound = len(s)
                
            if genomicStrand == laggingComplementStrand:            
                s = r.seq[lowSearchBound:highSearchBound]
                sx = r.seq[start:end]
            else:
                s = r.seq[lowSearchBound:highSearchBound].reverse_complement()
                sx = r.seq[start:end].reverse_complement()
            
            searchEnd = (end-start+searchSize+5)
            searchStart = searchSize + 5
            
            try:
                testSeqs = self.optimizeSecondaryStructure(s, searchStart, searchEnd, searchSize, cutOff, ratio = 2)
                testSeqs.sort()
                if len(testSeqs) == 0:
                    (score,adjust,foldScore,bestSeq) = ("na","na","na","na")
                else:
                    (score,adjust,foldScore,bestSeq) = testSeqs.pop(0)
            except:
                testSeqs = []    
                (score,adjust,foldScore,bestSeq) = ("na","na","na",sx)
                print "failed to exicute secondary structure test"
                #print "[%s]" % (s)
            
            #bestSeq = self.addStars(bestSeq,self.stars)
            if len(bestSeq) < boundary*2:
                print "Short Sequence"
            if self.verbose: print "%s best %s S:%s [%s] (%s)" % (len(testSeqs), adjust, score, foldScore, len(bestSeq))
            if self.verbose: print "[%s]" % (bestSeq)
            
            result.add(id,"original", originalString) 
            result.add(id,"hits", hits) 
            result.add(id,"genomic_start", genomicStart) 
            result.add(id,"genomic_end", genomicEnd) 
            result.add(id,"genomic_strand", genomicStrand) 
            result.add(id,"lagging_complement_strand", laggingComplementStrand) 
            result.add(id,"best", bestSeq) 
            result.add(id,"fold score", foldScore) 
            result.add(id,"off center", adjust) 
            
            #append to list of sequence regions
            sRegion = r
            sRegion.seq = sx
            sRegions.append(sRegion)
            
            index = index + 1
        
        logFile.close()
        
        return (targetMap,result,sRegions)
    
    def _selectFeatures(self,featureList,start,end):
        result = []
        for features in featureList:
            for feature in features:
                fstart = feature.location.start.position
                fend = feature.location.end.position
                if (start < fstart < end) | (start< fend < end):
                    result.append(feature)
                elif (fstart<start<fend) | (fstart<end<fend):
                    result.append(feature)
        return result
    
    def _matchSequenceLength(self,sequence,start,end,targetStart,targetEnd,bufferChar="X"):
        iseq = str(sequence)
        iseq = iseq.replace('*','')
        sdiff = int(start-targetStart)
        ediff = int(targetEnd - end)
        
        if sdiff > 0:
            for i in range(0,sdiff):
                iseq = bufferChar + iseq
        
        if ediff > 0:
            for i in range(0,ediff):
                iseq = iseq + bufferChar
        
        if sdiff < 0:
            iseq = iseq[-sdiff:]
        
        if ediff < 0:
            iseq = iseq[:ediff]
            
        return iseq
            
    
    def parseRecombination(self, oligoReport, blastRecords, logFile="recombination_log.txt",minDiff=0,scoreMin = 1e-3):
        seqTools = SequenceTools()
        output = open(logFile,"w")
        
        rowNames = oligoReport.returnRowNames()
        
        for rName in rowNames:
            
            oligo = oligoReport["original"][rName]
            ostart = float(oligoReport["genomic_start"][rName])
            oend = float(oligoReport["genomic_end"][rName])
            strand = oligoReport["genomic_strand"][rName]
            
            location = (ostart + oend) / 2
            start = location - 50
            end = location + 50
            
            output.write("\n--Recombinant Name--\n" + rName + "\n")
            output.write(oligo + "\n")
            output.write("[%s,%s]\n" % (str(start),str(end)))
            
            sBlastRecords = seqTools.blastSelection(blastRecords,start=start,end=end,scoreMin=scoreMin)

            for r in sBlastRecords:            
                name = r.query
                for alignment in r.alignments:
                    output.write("==" + name + "\n")
                    for hsp in alignment.hsps:
                        (ts,ss) = hsp.frame
                        targetStrand = ss
                        targetStart = hsp.sbjct_start
                        targetEnd = hsp.sbjct_end
                        o = Seq(oligo)
                        
                        if int(strand) != int(ts):
                            #pass
                            o = o.reverse_complement()
                            
                        output.write("(oligo %s subject %s target %s)\n" % (strand,ss,ts))
                        output.write("(%s %s)\n" % (str(targetStart),str(targetEnd)))
                        iOligo = self._matchSequenceLength(o,ostart,oend,start,end)
                        
                        iquery = self._matchSequenceLength(hsp.query,targetStart,targetEnd,start,end)
                        imatch = self._matchSequenceLength(hsp.match,targetStart,targetEnd,start,end)
                        isbjct = self._matchSequenceLength(hsp.sbjct,targetStart,targetEnd,start,end)
                        
                        aMatch = iOligo + "\n" + iquery + "\n" + imatch + "\n" + isbjct
                        
                        output.write(str(aMatch)+"\n\n")
                        location = FeatureLocation(start,end)
            
        output.close()
        return None
            
    def analyseRecombindationResults(self,oligoReport,featureList,outputFile="oligo_recomb_analysis.txt"):
        
        result = Report()
        result.extend(oligoReport)
        rowNames = oligoReport.returnRowNames()
        recordAlignment = SeqRecord(Seq(""))
        recordAlignment.features = featureList
        for rName in rowNames:
            output.write(rName + "\n")
            
            start = float(oligoReport["genomic_start"][rName])
            end = float(oligoReport["genomic_end"][rName])
            oligo = oligoReport["best"][rName]
            
            keyString = "\n(%s,%s) %s\n\n" % (start,end,oligo)
            output.write(keyString)
            
            subFeatures = self._selectFeatures(featureList,start,end)
            index = 0
            for feature in subFeatures:
                matchName = "match_" + str(index)
                id = feature.id
                qValue = feature.qualifiers["query"]
                sValue = feature.qualifiers["subject"]
                matchString = feature.qualifiers["alignment"]
                output.write(id + "\n")
                output.write(matchString + "\n\n")
                result.add(rName,matchName,qValue)
                index += 1
        output.close()
        return result
        
        
    
