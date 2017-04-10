from debug_MACRO import *
import pdb

if 0:
    Python_command = '/Library/Frameworks/Python.framework/Versions/Current/bin/python'
    BowtiePath = '/Users/Soheil/Documents/Bioinformatics/bowtie-0.12.8'
    Bowtie2Path = '/Users/Soheil/Documents/Bioinformatics/bowtie2-2.1.0'
    TopHatPath = '/Users/Soheil/Documents/Bioinformatics/tophat-2.0.9.OSX_x86_64'
    SamPath = '/Users/Soheil/Documents/Bioinformatics/samtools-0.1.19'
    BCFPath = '/Users/Soheil/Documents/Bioinformatics/samtools-0.1.19/bcftools/bcftools'
    SimulPath = '/Users/Soheil/Documents/Bioinformatics/Simulation/RNASeqReadSimulator-master/src'
    FA2FQ_address = '/Users/Soheil/Documents/Bioinformatics/Simulation/FA2FQ/fasta_to_fastq.pl'
    igvPath = '/Users/Soheil/Documents/Bioinformatics/IGV_2.2.13'
    rsemPath = '/Users/Soheil/Documents/Bioinformatics/rsem-1.2.7'
    GATKPath = '/Users/Soheil/Documents/Bioinformatics/GenomeAnalysisTK-2.5-2-gf57256b'
    PicardPath = '/Users/Soheil/Documents/Bioinformatics/picard-tools-1.102/picard-tools-1.102'
    # Default_Ref_Path = '/Users/Soheil/Documents/Bioinformatics/Simulation/Chr15-new'
    Default_Ref_Path = '/Users/Soheil/Dropbox/Transcriptome/SNP-Calling-Summer15/data'
else:    
    #modified paths -- server side
    isServerSide = True #False
    if isServerSide: #infolab server
        #"""
        Python_command = '/home/shunfu1/anaconda2/bin/python' #try: which python
        BowtiePath = '/home/shunfu1/bowtie-1.1.2//'
        Bowtie2Path = '/usr/bin/' #'bowtie2'
        TopHatPath = '/usr/local/bin/'
        SamPath = '/home/shunfu1/samtools/bin//'
        #BCFPath = '/Users/Soheil/Documents/Bioinformatics/samtools-0.1.19/bcftools/bcftools'
        SimulPath = '/home/shunfu1/software/RNASeqReadSimulator/src/'
        FA2FQ_address = '/home/shunfu1/SNP_Calling_Summer15/tools/FA2FQ/fasta_to_fastq.pl'
        #igvPath = '/Users/Soheil/Documents/Bioinformatics/IGV_2.2.13'
        rsemPath = '/home/shunfu1/SNP_Calling_Summer15/tools/rsem-1.2.22/'
        STAR_path = '/home/shunfu1/STAR-2.5.2a/bin/Linux_x86_64_static/' #used by rsem
        #GATKPath = '/Users/Soheil/Documents/Bioinformatics/GenomeAnalysisTK-2.5-2-gf57256b'
        #PicardPath = '/Users/Soheil/Documents/Bioinformatics/picard-tools-1.102/picard-tools-1.102'
        Default_Proj_Path = '/home/shunfu1/SNP_Calling_Summer15/'
        Default_Ref_Path_Root = '/data1/shunfu1/SNPCalling/'
        Default_Ref_Path = '/data1/shunfu1/SNPCalling/data/'
        Default_Genome_Sources = '/data1/shunfu1/SNPCalling/genome/'
        #"""
    else:
        pdb.set_trace()