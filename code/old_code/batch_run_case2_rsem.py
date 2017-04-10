"""
- case 2 pipeline (gatk)
    
- main() is non-parallel version

"""

from Address import *
from Synthesis import *
from ReadProcess import *
from count_read_lambda import *
from final_caller19thjuly_m import *
from debug_MACRO import *
from sim_stat import *
from snp_res_statistics import *
#from find_pos import *
from para_operations import *
import pdb
import subprocess

def do_gen_exp(BED_address, EXP_fn_m, EXP_fn_p, flag=True):
    if flag==True:
        [BED_sorted_address, exp_address_m] = BED2ExpressionLevel(BED_address, EXP_fn_m) # generate a random expression level file
        [BED_sorted_address, exp_address_p] = BED2ExpressionLevel(BED_address, EXP_fn_p)
    else:
        exp_address_m = Default_Ref_Path + EXP_fn_m #for test purpose
        exp_address_p = Default_Ref_Path + EXP_fn_p
        BED_sorted_address = Default_Ref_Path + '/hg19_chr15-UCSC-sorted.bed'
    return [exp_address_m, exp_address_p, BED_sorted_address]

def do_gen_cov(exp_address_m,
               exp_address_p,
               BED_sorted_address,
               COV_fn_m,
               COV_fn_p,
               Stat_m,
               Stat_p,
               flag=True):
    if flag==True:
        #old cov
        coverage_address_m = ExpressionLevel2Coverage(BED_sorted_address, exp_address_m, COV_fn_m, Stat_m) # Find the exonic parts with coverage not less than a ceryain value
        coverage_address_p = ExpressionLevel2Coverage(BED_sorted_address, exp_address_p, COV_fn_p, Stat_p)
    else:
        coverage_address_m = Default_Ref_Path + COV_fn_m #for test purpose
        coverage_address_p = Default_Ref_Path + COV_fn_p
        
    return [coverage_address_m, coverage_address_p]

def do_gen_fil_cov(Stat_m,
                   Stat_p,
                   COV_fn_m,
                   COV_fn_p,
                   cov_level=90,
                   flag=True):
    qt = [cov_level] #snp appears in high coverage region (%)
    if flag==True:        
        Stat_m.set_qt(qt) 
        Stat_m.get_acc_cov_hd(COV_fn_m)
        Stat_m.set_acc_cov_qt()
        line_cover_threshold = Stat_m.acc_cover_qt[0]
        COV_fn_filtered_m = COV_fn_m[:-4] + '_qt'+repr(qt[0])+'.txt'
        coverage_address_filtered_m = Default_Ref_Path + COV_fn_filtered_m
        #use filtered coverage.txt to make snps generated in high exp-level region
        FilterCoverage(coverage_address_m, coverage_address_filtered_m, line_cover_threshold) 
        
        Stat_p.set_qt(qt) 
        Stat_p.get_acc_cov_hd(COV_fn_p)
        Stat_p.set_acc_cov_qt()
        line_cover_threshold = Stat_p.acc_cover_qt[0]
        COV_fn_filtered_p = COV_fn_p[:-4] + '_qt'+repr(qt[0])+'.txt'
        coverage_address_filtered_p = Default_Ref_Path + COV_fn_filtered_p
        #use filtered coverage.txt to make snps generated in high exp-level region
        FilterCoverage(coverage_address_p, coverage_address_filtered_p, line_cover_threshold)        
    else:        
        COV_fn_filtered_m = COV_fn_m[:-4] + '_qt'+repr(qt[0])+'.txt'
        coverage_address_filtered_m = Default_Ref_Path + COV_fn_filtered_m
        
        COV_fn_filtered_p = COV_fn_p[:-4] + '_qt'+repr(qt[0])+'.txt'
        coverage_address_filtered_p = Default_Ref_Path + COV_fn_filtered_p
    return [coverage_address_filtered_m, coverage_address_filtered_p]

def do_gen_tar(ref_address,
               coverage_address_filtered_m,
               coverage_address_filtered_p,
               Num_SNP,
               tar_address_m,
               tar_address_p,
               SNP_address_m,
               SNP_address_p,
               flag=True):
    if flag==True:   
        GenTarget(ref_address, coverage_address_filtered_m, Num_SNP, tar_address_m, SNP_address_m ) # Generate 2 random target2 (for paternal and maternal) sequence from ref_address by insering Num_SNP SNPs in the exonic positions
        GenTarget(ref_address, coverage_address_filtered_p, Num_SNP, tar_address_p, SNP_address_p )
    return

def do_gen_reads(BED_address,
                 tar_address_m,
                 tar_address_p,
                 exp_address_m,
                 exp_address_p,
                 N,
                 L,
                 error_rate,
                 flag=True):
    if flag==True:
        [readBED_address_m, readFA_address_m, readFQ_address_m] = ReadGeneration(tar_address_m,
                                                                                 BED_address,
                                                                                 exp_address_m,
                                                                                 N,
                                                                                 L,
                                                                                 error_rate)
        [readBED_address_p, readFA_address_p, readFQ_address_p] = ReadGeneration(tar_address_p,
                                                                                 BED_address,
                                                                                 exp_address_p,
                                                                                 N,
                                                                                 L,
                                                                                 error_rate)
        readFQ_address = Default_Ref_Path + '/Tar_read_l'+ repr(L) +'.fastq'  #for test purpose
        subprocess.call('cat '  + readFQ_address_m + ' ' + readFQ_address_p + ' > ' + readFQ_address, shell=True )
    else:
        readBED_address_m = Default_Ref_Path + '/Tar_m_read_l'+ repr(L) +'.bed'
        readFA_address_m = Default_Ref_Path + '/Tar_m_read_l'+ repr(L) +'.fasta'
        readFQ_address_m = Default_Ref_Path + '/Tar_m_read_l'+ repr(L) +'.fastq' #for test purpose
        
        readBED_address_p = Default_Ref_Path + '/Tar_p_read_l'+ repr(L) +'.bed'
        readFA_address_p = Default_Ref_Path + '/Tar_p_read_l'+ repr(L) +'.fasta'
        readFQ_address_p = Default_Ref_Path + '/Tar_p_read_l'+ repr(L) +'.fastq' #for test purpose
        
        readFQ_address = Default_Ref_Path + '/Tar_read_l'+ repr(L) +'.fastq'  #for test purpose
    return [readFQ_address]

def main(case = 2,
         para_comp = False,
         num_p = 1,
         Threshold_num_reads=1):
    
    print('Default data folder is %s'%Default_Ref_Path)
    
    Stat_m = sim_stat('/sim_stat_dmp_m.txt')
    Stat_p = sim_stat('/sim_stat_dmp_p.txt')
    
    ref_address = Default_Ref_Path + '/Chr15.fa'
    BED_address = Default_Ref_Path + '/hg19_chr15-UCSC.bed'
    
    EXP_fn_m = '/exp_m.txt' #for true exp
    EXP_fn_p = '/exp_p.txt'
    
    COV_fn_m = '/coverage_m.txt'
    COV_fn_p = '/coverage_p.txt'
    
    Num_SNP = 20 #1000 #20 #for true snp
    
    tar_address_m = Default_Ref_Path + '/Tar_m.txt'
    tar_address_p = Default_Ref_Path + '/Tar_p.txt'

    SNP_address_m = Default_Ref_Path + '/SNP_m.txt' 
    SNP_address_p = Default_Ref_Path + '/SNP_p.txt'
    
    N=100000 #10000000 #100000 #Number of Reads
    L=100       #Read Length
    error_rate = 0
    
    [exp_address_m, exp_address_p, BED_sorted_address] = do_gen_exp(BED_address,
                                                                    EXP_fn_m,
                                                                    EXP_fn_p,
                                                                    flag=False)
    
    [coverage_address_m, coverage_address_p] = do_gen_cov( exp_address_m,
                                                           exp_address_p,
                                                           BED_sorted_address,
                                                           COV_fn_m,
                                                           COV_fn_p,
                                                           Stat_m,
                                                           Stat_p,
                                                           flag=False)
                                                           
    [coverage_address_filtered_m, coverage_address_filtered_p] = do_gen_fil_cov(Stat_m,
                                                                                Stat_p,
                                                                                COV_fn_m,
                                                                                COV_fn_p,
                                                                                cov_level=90,
                                                                                flag=False)
    do_gen_tar(ref_address,
               coverage_address_filtered_m,
               coverage_address_filtered_p,
               Num_SNP,
               tar_address_m,
               tar_address_p,
               SNP_address_m,
               SNP_address_p,
               flag=False)

    [readFQ_address] = do_gen_reads( BED_address,
                                     tar_address_m,
                                     tar_address_p,
                                     exp_address_m,
                                     exp_address_p,
                                     N,
                                     L,
                                     error_rate,
                                     flag=False)
    
    ####
    [Res_addr] = do_gatk_best_practice()
    #Res_addr = '/data1/shunfu1/SNPCalling/data/data_GATK/GATK_out/'+"/raw_variants.vcf"
    ####
    
    #pdb.set_trace()
    #vcf file converted here
    snp_res_stat = do_snp_res_statistics(Res_addr, SNP_address_m, SNP_address_p, Res_addr[:-4]+'_[snp_res_stat].txt')
    
    return

def do_gatk_best_practice():
    
    Ref_dir=Default_Proj_Path
    #"/home/sreeramkannan/singleCell/SNP_Calling_Summer15/"
    GATK_addr=Ref_dir+"/tools/GATK/GenomeAnalysisTK.jar"
    #"/home/sreeramkannan/singleCell/SNP_Calling_Summer15/tools/GATK/GenomeAnalysisTK.jar"
    Picard_addr=Ref_dir+"/tools/picard-tools-1.138/picard.jar"
    #"/home/sreeramkannan/singleCell/SNP_Calling_Summer15/tools/picard-tools-1.138/picard.jar"
    
    #pdb.set_trace()
    data_GATK_addr=Default_Ref_Path+"/data_GATK_useGTF/"#"/data_GATK_rsem/"    
    cmd = "mkdir -p "+data_GATK_addr
    subprocess.call(cmd, shell=True)
        
    fa_addr=Default_Ref_Path+"/Chr15.fa"
    Genome_addr = data_GATK_addr+"/Chr15.fa"
    cmd = "cp "+fa_addr+" "+ Genome_addr
    subprocess.call(cmd, shell=True)

    #pdb.set_trace()
    gtf_addr = Default_Ref_Path+"/hg19_chr15-UCSC.gtf"
    
    Read_addr=Default_Ref_Path+"/Tar_read_l100.fastq"
    Genome_dict_addr=data_GATK_addr+"/Chr15.dict"
    
    #"""
    #1. mapping to reference
    genomeDir=data_GATK_addr+"/Chr15/"
    cmd = "mkdir -p "+genomeDir
    subprocess.call(cmd, shell=True)
    cmd =   "STAR --runThreadN 20 "+\
            "--runMode genomeGenerate "+\
            "--genomeDir "+genomeDir+" "+\
            "--genomeFastaFiles "+Genome_addr+" "+\
            "--sjdbGTFfile "+gtf_addr
    #pdb.set_trace()
    subprocess.call(cmd, shell=True)
    
    #2.
    runDir=data_GATK_addr+"/1pass/"
    cmd = "mkdir -p "+runDir
    subprocess.call(cmd, shell=True)
    #cmd = "cd "+runDir
    #subprocess.call(cmd, shell=True)
    cmd =   "STAR --runThreadN 20 --genomeDir "+genomeDir+\
            " --readFilesIn "+Read_addr+\
            " --outFileNamePrefix "+runDir
    #pdb.set_trace()
    subprocess.call(cmd, shell=True)
    
    #3. mapping to reference (2pass)
    genomeDir_2pass=data_GATK_addr+"/Chr15_2pass/"
    cmd = "mkdir -p "+genomeDir_2pass
    subprocess.call(cmd, shell=True)
    cmd =   "STAR --runThreadN 20 --runMode genomeGenerate "+\
            "--genomeDir "+genomeDir_2pass+" "+\
            "--genomeFastaFiles "+ Genome_addr+" "+\
            "--sjdbGTFfile "+gtf_addr+" "\
            "--sjdbFileChrStartEnd "+ runDir+"/SJ.out.tab "+\
            "--sjdbOverhang 75"
    #pdb.set_trace()
    subprocess.call(cmd, shell=True)
    
    #4.
    runDir_2pass=data_GATK_addr+"/2pass/"
    cmd = "mkdir -p "+runDir_2pass
    subprocess.call(cmd, shell=True)
    #cmd = "cd "+runDir_2pass
    #subprocess.call(cmd, shell=True)
    cmd =   "STAR --runThreadN 20 --genomeDir "+genomeDir_2pass+" "\
            "--readFilesIn "+Read_addr+\
            " --outFileNamePrefix "+runDir_2pass
    #pdb.set_trace()
    subprocess.call(cmd, shell=True)
    
    #5. use Picard to:
    #pdb.set_trace() 
    #sam_input = "/data1/shunfu1/SNPCalling/data/rsem/Chr15.genome.sorted_n.sam"
    sam_input = runDir_2pass+"/Aligned.out.sam "

    cmd ="java -jar "+Picard_addr+" AddOrReplaceReadGroups "+\
         "I=%s "%sam_input+\
         "O="+runDir_2pass+"/rg_added_sorted.bam "+\
         "SO=coordinate "+\
         "RGID=id "+\
         "RGLB=library "+\
         "RGPL=platform "+\
         "RGPU=machine "+\
         "RGSM=sample"
    subprocess.call(cmd, shell=True)
    
    #runDir_2pass+"/dedupped.bam " --> !!bam file used by our caller
    cmd ="java -jar "+Picard_addr+" MarkDuplicates "+\
         "I="+runDir_2pass+"/rg_added_sorted.bam "+\
         "O="+runDir_2pass+"/dedupped.bam "+\
         "CREATE_INDEX=true "+\
         "VALIDATION_STRINGENCY=SILENT "+\
         "M="+runDir_2pass+"/output.metrics"
    subprocess.call(cmd, shell=True)
    
    #6. split n cigar reads
    cmd =   "java -jar "+Picard_addr+" CreateSequenceDictionary "+\
            "R="+Genome_addr+" "+\
            "O="+Genome_dict_addr
    subprocess.call(cmd, shell=True)
    
    cmd = "samtools faidx "+Genome_addr
    subprocess.call(cmd, shell=True)
    
    cmd = "java -jar "+GATK_addr+" "+\
          "-T SplitNCigarReads "+\
          "-R "+Genome_addr+" "+\
          "-I "+runDir_2pass+"/dedupped.bam "+\
          "-o "+runDir_2pass+"/split.bam "+\
          "-rf ReassignOneMappingQuality "+\
          "-RMQF 255 "+\
          "-RMQT 60 "+\
          "-U ALLOW_N_CIGAR_READS"
    
    subprocess.call(cmd, shell=True)
    #"""
    
    #7. Variant Caller 
    #pdb.set_trace()
    
    data_GATK_out_addr=data_GATK_addr+"/GATK_out/"
    cmd = "mkdir -p "+data_GATK_out_addr
    subprocess.call(cmd, shell=True) 

    runDir_2pass=data_GATK_addr+"/2pass/"
    Res_addr=data_GATK_out_addr+"/raw_variants.vcf"
    
    cmd =   "java -jar "+GATK_addr+" "+\
            "-T HaplotypeCaller "+\
            "-R "+Genome_addr+" "+\
            "-I "+runDir_2pass+"/split.bam "+\
            "-dontUseSoftClippedBases "+\
            "-stand_emit_conf 20 "+\
            "-stand_call_conf 20 "+\
            "-o "+Res_addr
    subprocess.call(cmd, shell=True) 
    
    return [Res_addr]

'''
modified based on do_gatk_best_practice,
more flexible on input parameter & instructures

usage:
--do_gatk_best_practice_modi
--DefRefPath drp #replace the one of Address.py
--GATKdir gd #default data_GATK only a folder name, put related alignments & res here, under drp
--Genome fa #e.g. address of ref genome, no mutation
[--SkipPrepareGenomeRefIdx] #skip 1pass STAR --runMode genomeGenerate if done previously
--GTF gtfAddr #may contain more annotations (e.g. hg19.gtf while chr15.fa used)
--Rd readAddress #SE used currently (to extend to --Rd1 --Rd2)
[--toolFolder toolFolder] #GATK tools' location

example:
python batch_run_case2_rsem.py --do_gatk_best_practice_modi --DefRefPath /data1/shunfu1/SNPCalling/data/ --GATKdir data_GATK_chr15_reads100k --Genome /data1/shunfu1/SNPCalling/data/genome/chr15.fa --GTF /data1/shunfu1/SNPCalling/data/gene_annotation/gencode.v19.annotation.gtf --Rd /data1/shunfu1/SNPCalling/data/sim_reads_100k.fq

output:
  <drp>/<gd>/<genome>
            /1pass/
            /<genome>_2pass/
            /2pass/ (alignment files)
            /GATK_out/ (vcf file)
'''
def do_gatk_best_practice_modi(args):

    #pdb.set_trace()

    Ref_dir=Default_Proj_Path #above data folder

    if '--toolFolder' in args:
        toolFolder = args[args.index('--toolFolder')+1]
        GATK_addr=toolFolder+"/GATK/GenomeAnalysisTK.jar"
        Picard_addr=toolFolder+"/picard-tools-1.138/picard.jar"
    else:
        GATK_addr=Ref_dir+"/tools/GATK/GenomeAnalysisTK.jar"
        Picard_addr=Ref_dir+"/tools/picard-tools-1.138/picard.jar"
        #Picard_addr=Ref_dir+"/tools/picard-tools-2.9.0/picard.jar"

    ########## Parameters

    if '--DefRefPath' in args:
        drp=args[args.index('--DefRefPath')+1]
    else:
        drp=Default_Ref_Path

    if '--GATKdir' in args:
        gd = args[args.index('--GATKdir')+1]
    else:
        gd = 'data_GATK'

    data_GATK_addr=drp+'/%s/'%gd #e.g. "/data_GATK_useGTF/"#"/data_GATK_rsem/"    
    cmd = "mkdir -p "+data_GATK_addr
    subprocess.call(cmd, shell=True)   

    if '--Genome' in args:
        fa_addr = args[args.index('--Genome')+1]
        Genome_addr = fa_addr
        GenomeName = Genome_addr.split('/')[-1][:-3] #skip .fa
        Genome_dict_addr=Genome_addr[:-3]+'.dict' #data_GATK_addr+"/Chr15.dict" 
    else:
        fa_addr=drp+"/Chr15.fa"
        Genome_addr = data_GATK_addr+"/Chr15.fa"
        cmd = "cp "+fa_addr+" "+ Genome_addr
        subprocess.call(cmd, shell=True)
        GenomeName = Genome_addr.split('/')[-1][:-3]
        Genome_dict_addr=Genome_addr[:-3]+'.dict' #data_GATK_addr+"/Chr15.dict" 

    if '--SkipPrepareGenomeRefIdx' in args:
        skip = 1
    else:
        skip = 0

    if '--GTF' in args:
        gtf_addr = args[args.index('--GTF')+1]
    else:
        gtf_addr = drp+"/hg19_chr15-UCSC.gtf"
    
    if '--Rd' in args:
        Read_addr = args[args.index('--Rd')+1]
    else:
        Read_addr = drp+"/Tar_read_l100.fastq"

    ########## GATK best practice    
    
    #"""
    #1. mapping to reference
    genomeDir=data_GATK_addr+"/%s/"%GenomeName
    cmd = "mkdir -p "+genomeDir
    subprocess.call(cmd, shell=True)
    cmd =   "STAR --runThreadN 20 "+\
            "--runMode genomeGenerate "+\
            "--genomeDir "+genomeDir+" "+\
            "--genomeFastaFiles "+Genome_addr+" "+\
            "--sjdbGTFfile "+gtf_addr
    #pdb.set_trace()
    if skip==0: subprocess.call(cmd, shell=True)
    
    #2.
    runDir=data_GATK_addr+"/1pass/"
    cmd = "mkdir -p "+runDir
    subprocess.call(cmd, shell=True)
    #cmd = "cd "+runDir
    #subprocess.call(cmd, shell=True)
    cmd =   "STAR --runThreadN 20 --genomeDir "+genomeDir+\
            " --readFilesIn "+Read_addr+\
            " --outFileNamePrefix "+runDir
    #pdb.set_trace()
    subprocess.call(cmd, shell=True)
    
    #3. mapping to reference (2pass)
    genomeDir_2pass=data_GATK_addr+"/%s_2pass/"%GenomeName
    cmd = "mkdir -p "+genomeDir_2pass
    subprocess.call(cmd, shell=True)
    cmd =   "STAR --runThreadN 20 --runMode genomeGenerate "+\
            "--genomeDir "+genomeDir_2pass+" "+\
            "--genomeFastaFiles "+ Genome_addr+" "+\
            "--sjdbGTFfile "+gtf_addr+" "\
            "--sjdbFileChrStartEnd "+ runDir+"/SJ.out.tab "+\
            "--sjdbOverhang 75"
    #pdb.set_trace()
    subprocess.call(cmd, shell=True)
    
    #4.
    runDir_2pass=data_GATK_addr+"/2pass/"
    cmd = "mkdir -p "+runDir_2pass
    subprocess.call(cmd, shell=True)
    #cmd = "cd "+runDir_2pass
    #subprocess.call(cmd, shell=True)
    cmd =   "STAR --runThreadN 20 --genomeDir "+genomeDir_2pass+" "\
            "--readFilesIn "+Read_addr+\
            " --outFileNamePrefix "+runDir_2pass
    #pdb.set_trace()
    subprocess.call(cmd, shell=True)
    
    #5. use Picard to:
    #pdb.set_trace() 
    #sam_input = "/data1/shunfu1/SNPCalling/data/rsem/Chr15.genome.sorted_n.sam" #cross check
    sam_input = runDir_2pass+"/Aligned.out.sam "

    cmd ="java -jar "+Picard_addr+" AddOrReplaceReadGroups "+\
         "I=%s "%sam_input+\
         "O="+runDir_2pass+"/rg_added_sorted.bam "+\
         "SO=coordinate "+\
         "RGID=id "+\
         "RGLB=library "+\
         "RGPL=platform "+\
         "RGPU=machine "+\
         "RGSM=sample"
    subprocess.call(cmd, shell=True)
    
    #runDir_2pass+"/dedupped.bam " --> !!bam file used by our caller
    cmd ="java -jar "+Picard_addr+" MarkDuplicates "+\
         "I="+runDir_2pass+"/rg_added_sorted.bam "+\
         "O="+runDir_2pass+"/dedupped.bam "+\
         "CREATE_INDEX=true "+\
         "VALIDATION_STRINGENCY=SILENT "+\
         "M="+runDir_2pass+"/output.metrics"
    subprocess.call(cmd, shell=True)
    
    #6. split n cigar reads
    if os.path.exists(Genome_dict_addr) == False:
        cmd =   "java -jar "+Picard_addr+" CreateSequenceDictionary "+\
                "R="+Genome_addr+" "+\
                "O="+Genome_dict_addr
        subprocess.call(cmd, shell=True)
    
    cmd = "samtools faidx "+Genome_addr
    subprocess.call(cmd, shell=True)

    #pdb.set_trace()
    
    cmd = "java -jar "+GATK_addr+" "+\
          "-T SplitNCigarReads "+\
          "-R "+Genome_addr+" "+\
          "-I "+runDir_2pass+"/dedupped.bam "+\
          "-o "+runDir_2pass+"/split.bam "+\
          "-rf ReassignOneMappingQuality "+\
          "-RMQF 255 "+\
          "-RMQT 60 "+\
          "-U ALLOW_N_CIGAR_READS"
    
    subprocess.call(cmd, shell=True)
    #"""
    
    #7. Variant Caller 
    #pdb.set_trace()
    
    data_GATK_out_addr=data_GATK_addr+"/GATK_out/"
    cmd = "mkdir -p "+data_GATK_out_addr
    subprocess.call(cmd, shell=True) 

    Res_addr=data_GATK_out_addr+"/raw_variants.vcf"
    
    cmd =   "java -jar "+GATK_addr+" "+\
            "-T HaplotypeCaller "+\
            "-R "+Genome_addr+" "+\
            "-I "+runDir_2pass+"/split.bam "+\
            "-dontUseSoftClippedBases "+\
            "-stand_emit_conf 20 "+\
            "-stand_call_conf 20 "+\
            "-o "+Res_addr
    subprocess.call(cmd, shell=True) 
    
    return [Res_addr]

if __name__ == "__main__":
    
    #pdb.set_trace()

    args = sys.argv

    if '--do_gatk_best_practice_modi' in args:
        do_gatk_best_practice_modi(args)
    else:    
        main()
