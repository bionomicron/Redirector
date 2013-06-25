"""
Initial analysis scripts.
"""

import os,re
import subprocess
from optparse import OptionParser

class BowTieTools:
    
    def __init__(self):
        # Location of this script. We may find other paths relative to this.
        PWD = os.path.dirname(os.path.realpath(__file__ ))

        # Main root of data.
        DATA_ROOT = '/home/nick/Data/cccb.dfci.harvard.edu/frd/Nicholas_Guido_C23YAACXX/demux'
        
        
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

# Location of this script. We may find other paths relative to this.
PWD = os.path.dirname(os.path.realpath(__file__ ))

# Main root of data.
DATA_ROOT = '/home/nick/Data/cccb.dfci.harvard.edu/frd/Nicholas_Guido_C23YAACXX/demux'

# Root of where to put output.
OUTPUT_ROOT = os.path.join(PWD, 'data/output')

# Illumina fastq files.
# NOTE: We assume we are running on gmc.ant.
FASTQ_ROOT = ('/home/nick/Data/cccb.dfci.harvard.edu/frd/Nicholas_Guido_C23YAACXX/demux')

FASTQ_S1 = 'Sample.ACAGTG.R1.fastq.gz'
FASTQ_S2 = 'Sample.ACAGTG.R2.fastq.gz'
#FASTQ_S3 = 'GEN009334_S3_L001_R1_001.fastq.gz'
#FASTQ_S4 = 'GEN009335_S4_L001_R1_001.fastq.gz'
#FASTQ_S5 = 'GEN009336_S5_L001_R1_001.fastq.gz'
#FASTQ_UNDETERMINED = 'Undetermined_S0_L001_R1_001.fastq.gz'

ALL_FASTQ = [
    (FASTQ_S1,FASTQ_S2),
#    FASTQ_S3,
#    FASTQ_S4,
#    FASTQ_S5,
#    FASTQ_UNDETERMINED
]


def run(reference_seqs, analysis_output_dir):
    # Create the bowtie2 index for the reference genome if hasn't been
    # created yet. We check for the existence of a particular file.
    bowtie2_index_test_file = os.path.splitext(reference_seqs)[0] + '.1.bt2'
    if not os.path.exists(bowtie2_index_test_file):
        bowtie2_build_index(reference_seqs)

    # Create the output dir if it doesn't exist.
    if not os.path.exists(analysis_output_dir):
        os.mkdir(analysis_output_dir)

    # Run bowtie2 with a guess at the best parameters and bucket the
    # results.
    for fastq_filename_list in ALL_FASTQ:
        current_fastq1 = os.path.join(FASTQ_ROOT, fastq_filename_list[0])
	current_fastq2 = os.path.join(FASTQ_ROOT, fastq_filename_list[1])
        output_sam_filename = os.path.splitext(os.path.splitext(
                fastq_filename_list[0])[0])[0] + '.bowtie2_align.sam'
        output_sam_path = os.path.join(analysis_output_dir, output_sam_filename)
        align_with_bowtie2(reference_seqs, current_fastq1, current_fastq2, output_sam_path)

        # Bucket the results.
        print '...Bucketing ' + os.path.split(output_sam_path)[1]
        bucket_alignments(output_sam_path)

    # 3. Possibly tweak the parameters to figure out if there is a better solution.


###############################################################################
# Bowtie2 Indexing
###############################################################################

def bowtie2_build_index(ref_genome_location):
    """Build the index for bowtie2."""
    index_prefix = _get_bowtie2_index_prefix(ref_genome_location)

    # Create the bowtie2 index.
    subprocess.check_call([
        'bowtie2-build',
        ref_genome_location,
        index_prefix
    ])


def _get_bowtie2_index_prefix(ref_genome_location):
    return os.path.splitext(ref_genome_location)[0]


###############################################################################
# Bowtie2 Alignment
###############################################################################

def align_with_bowtie2(ref_fasta, input_reads_1_fq, input_reads_2_fq, output_sam):
    """
	Perform alignment using bowtie2.
    """
    index_prefix = _get_bowtie2_index_prefix(ref_fasta)

    stderr_log_path = output_sam + '.log'

    with open(output_sam, 'w') as output_fh:
        with open(stderr_log_path, 'w') as stderr_log_fh:
            THREADS = '4'
            args = [
                'bowtie2',
                '--local',
                '--fast-local',
                '-q',
                '-p', THREADS,
                #'--phred33',
                #'-k', '3', # Search for multiple distinct alignments
                index_prefix,
                '-1', input_reads_1_fq,
		'-2', input_reads_2_fq,
            ]
            subprocess.call(args, stdout=output_fh, stderr=stderr_log_fh)


###############################################################################
# Analyzing alignments: bucketing, etc.
###############################################################################

def bucket_alignments(sam_file_path):
    """
	Sort the alignments into buckets according to the genomes they aligned against.
    """
    args = "grep -v '^@' " + sam_file_path + ' | cut -f 3 | sort | uniq -c'

    buckets_filename = sam_file_path + '.buckets'

    with open(buckets_filename, 'w') as fh:
        subprocess.call(args, shell=True, stdout=fh)

from optparse import OptionParser

if __name__ == '__main__':
    '''
    REFERENCE_SEQS_0 = os.path.join(DATA_ROOT,
            'FattyAcid_sequencing_targets.fasta')
    ANALYSIS_OUTPUT_DIR_0 = os.path.join(OUTPUT_ROOT,
            'analysis_output_2013_06_03')
    run(REFERENCE_SEQS_0, ANALYSIS_OUTPUT_DIR_0)
    '''
    
    parser = OptionParser()
    
    parser.add_option("-v", "--verbose",
                      action="store_true", 
                      dest="verbose", 
                      default=False,
                      help="set verbose mode")
    
    parser.add_option("-c","--config", 
                      dest="configFile",
                      default = "seq_builder.config",
                      help="name of sequencing analysis configuration file", 
                      metavar="FILE")
    
    parser.add_option("-m", "--Mode", 
                      dest="mode", 
                      default=None,
                      metavar="String",
                      help="Sequence Analysis Mode (consensus,report)")
    
    parser.add_option("--working_dir", 
                      dest="workFolder", 
                      default="./seq_work/",
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
    #verbose = options.verbose
    #mode = options.mode
    #config = options.config
    
    workFolder = options.workFolder
    refTag = options.genomicSeqTag
    readTag1 = options.readTag1
    readTag2 = options.readTag2
    
    refTag = "NC_000913_2"
    readTag1 = "Sample.GATCAG.R1"
    readTag2 = "Sample.GATCAG.R2"
    
    
    if not os.path.exists(workFolder):
        print "making work folder %s" % (workFolder)
        os.mkdir(workFolder)
    
    readTag = readTag1
    drefTag = workFolder + refTag
    dreadTag = workFolder + readTag1
    dreadTag1 = workFolder + readTag1
    dreadTag2 = workFolder + readTag2
    
    print "building reference index"
    call = "bowtie2-build %s %s > bowtie2_ref_index_log.txt" % (refTag,drefTag)
    print call
    #subprocess.call(call,shell=True)
    
    print "aligning sequences to reference"
    call = "bowtie2 -q -p 4 --mp 1,1- --end-to-end %s -1 %s.fastq.gz -2 %s.fastq.gz -S %s.sam" % (drefTag,readTag1,readTag2,dreadTag1)
    print call
    #subprocess.call(call,shell=True)
    
    print "creating bam file"
    call = "samtools view -S -b %s.sam -o %s.bam" % (dreadTag,dreadTag)
    print call
    #subprocess.call(call,shell=True)
    
    print "sorting bam file"
    call = "samtools sort %s.bam %s_s" % (dreadTag,dreadTag)
    print call
    #subprocess.call(call,shell=True)
    
    print "indexing bam file"
    call = "samtools index %s_s.bam" % (dreadTag)
    print call
    #subprocess.call(call,shell=True)
    
    print "finding VCF information"
    call = "freebayes --fasta-reference %s %s_s.bam -v %s_%s.vcf" % (refTag,dreadTag,refTag,readTag)
    print call
    subprocess.call(call,shell=True)
    

