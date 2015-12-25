from debug_MACRO import *
import pdb

#pdb.set_trace()

if en_debug==0:
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
    isServerSide = False #False
    ServerNum = 1
    
    if isServerSide and ServerNum==1:
        #"""
        Python_command = '/usr/bin/python' #try: which python
        BowtiePath = '/usr/local/bin/'
        Bowtie2Path = '/usr/local/bin/' #'bowtie2'
        TopHatPath = '/usr/local/bin/'
        SamPath = '/usr/local/bin/' #'/Users/Soheil/Documents/Bioinformatics/samtools-0.1.19'
        #BCFPath = '/Users/Soheil/Documents/Bioinformatics/samtools-0.1.19/bcftools/bcftools'
        SimulPath = '/home/sreeramkannan/singleCell/SNP_Calling_Summer15/tools/RNASeqReadSimulator-master/src'
        FA2FQ_address = '/home/sreeramkannan/singleCell/SNP_Calling_Summer15/tools/FA2FQ/fasta_to_fastq.pl'
        #igvPath = '/Users/Soheil/Documents/Bioinformatics/IGV_2.2.13'
        #rsemPath = '/usr/local/bin/'
        rsemPath = '/home/sreeramkannan/singleCell/SNP_Calling_Summer15/tools/rsem-1.2.22/'
        #GATKPath = '/Users/Soheil/Documents/Bioinformatics/GenomeAnalysisTK-2.5-2-gf57256b'
        #PicardPath = '/Users/Soheil/Documents/Bioinformatics/picard-tools-1.102/picard-tools-1.102'
        # Default_Ref_Path = '/Users/Soheil/Documents/Bioinformatics/Simulation/Chr15-new'
        Default_Proj_Path = '/home/sreeramkannan/singleCell/SNP_Calling_Summer15/'
        #Default_Ref_Path = '/home/sreeramkannan/singleCell/SNP_Calling_Summer15/data_0814_stat/'
        #Default_Ref_Path = '/home/sreeramkannan/singleCell/SNP_Calling_Summer15/data_0814/'
        #Default_Ref_Path = '/home/sreeramkannan/singleCell/SNP_Calling_Summer15/data_SNP1k_Reads10M/'
        #Default_Ref_Path = '/home/sreeramkannan/singleCell/SNP_Calling_Summer15/data_SNP20_Reads100K/'
        Default_Ref_Path = '/home/sreeramkannan/singleCell/SNP_Calling_Summer15/data_SNP20_Reads100K_diffMPExp/'
        #Default_Ref_Path = '/home/sreeramkannan/singleCell/SNP_Calling_Summer15/data_SNP20_Reads100K_debug_line_cover/'
        #"""  
    elif isServerSide and ServerNum==2:
        #"""
        Python_command = '/usr/bin/python' #try: which python
        BowtiePath = '/usr/bin/'
        Bowtie2Path = '/usr/bin/' #'bowtie2'
        TopHatPath = '/usr/local/bin/'
        SamPath = '/usr/bin/' #'/Users/Soheil/Documents/Bioinformatics/samtools-0.1.19'
        #BCFPath = '/Users/Soheil/Documents/Bioinformatics/samtools-0.1.19/bcftools/bcftools'
        SimulPath = '/home/sreeramkannan/singleCell/SNP_Calling_Summer15/tools/RNASeqReadSimulator-master/src'
        FA2FQ_address = '/home/sreeramkannan/singleCell/SNP_Calling_Summer15/tools/FA2FQ/fasta_to_fastq.pl'
        #igvPath = '/Users/Soheil/Documents/Bioinformatics/IGV_2.2.13'
        rsemPath = '/usr/local/bin/'
        #GATKPath = '/Users/Soheil/Documents/Bioinformatics/GenomeAnalysisTK-2.5-2-gf57256b'
        #PicardPath = '/Users/Soheil/Documents/Bioinformatics/picard-tools-1.102/picard-tools-1.102'
        # Default_Ref_Path = '/Users/Soheil/Documents/Bioinformatics/Simulation/Chr15-new'
        Default_Proj_Path = '/home/sreeramkannan/singleCell/SNP_Calling_Summer15/'
        #Default_Ref_Path = '/home/sreeramkannan/singleCell/SNP_Calling_Summer15/data_0814_stat/'
        #Default_Ref_Path = '/home/sreeramkannan/singleCell/SNP_Calling_Summer15/data_0814/'
        Default_Ref_Path = '/home/sreeramkannan/singleCell/SNP_Calling_Summer15/data_SNP1k_Reads10M/'
        #Default_Ref_Path = '/home/sreeramkannan/singleCell/SNP_Calling_Summer15/data_SNP20_Reads100K_para/'
        #Default_Ref_Path = '/home/sreeramkannan/singleCell/SNP_Calling_Summer15/data_0814_copy/'
        #"""  
    else:
        #modified paths -- note PC linux side
        Python_command = '/home/olivo/anaconda/bin/python' #try: which python
        #Python_command = '/usr/bin/python' 
        BowtiePath = '/home/olivo/bio_tools/bowtie-1.1.2/'
        Bowtie2Path = '/home/olivo/bio_tools/bowtie2-2.2.4/' #'bowtie2'
        TopHatPath = '/home/olivo/bio_tools/tophat-2.1.0.Linux_x86_64/'
        SamPath = '/home/olivo/bio_tools/samtools-1.2/' #'/Users/Soheil/Documents/Bioinformatics/samtools-0.1.19'
        #BCFPath = '/Users/Soheil/Documents/Bioinformatics/samtools-0.1.19/bcftools/bcftools'
        SimulPath = '/home/olivo/bio_tools/RNASeqReadSimulator-master/src/'
        FA2FQ_address = '/home/olivo/bio_tools/FA2FQ/fasta_to_fastq.pl'
        #igvPath = '/Users/Soheil/Documents/Bioinformatics/IGV_2.2.13'
        rsemPath = '/home/olivo/bio_tools/rsem-1.2.22/'
        #GATKPath = '/Users/Soheil/Documents/Bioinformatics/GenomeAnalysisTK-2.5-2-gf57256b'
        #PicardPath = '/Users/Soheil/Documents/Bioinformatics/picard-tools-1.102/picard-tools-1.102'
        # Default_Ref_Path = '/Users/Soheil/Documents/Bioinformatics/Simulation/Chr15-new'
        #Default_Ref_Path = '/home/olivo/Desktop/SNP-Calling-Summer15/data_0806_modi/'
        #Default_Proj_Path = '/home/olivo/Desktop/SNP-Calling-Summer15/'
        #Default_Ref_Path = '/home/olivo/Desktop/SNP-Calling-Summer15/data_0814_stat/' #en_debug_0814
        Default_Ref_Path = '/home/olivo/Desktop/SNP-Calling-Summer15/data_SNP20_Reads100K/' #data_0814_stat/snp_analysis/use_data_0814/' #en_debug_0814
        Default_Ref_Path = '/home/olivo/Documents/SNP-Calling/' 
        #"""
