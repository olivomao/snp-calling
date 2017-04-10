"""
- case 1 & case 6 pipeline combined
    -- case 6: mainly processing on read alignment
    -- case 1: mainly processing on snp caller & filtering
    
- main() is non-parallel version

- 1/24: use filter_snp_lam_half_filt2 (altCount expression level focuses on snp base)
- 1/31: use 'dedupped.sam' instead of 'split.sam' to avoid 'H' pattern in Cigar String

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
from util import *

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

#coverage_address = merge_coverage(coverage_address_m, coverage_address_p)
def merge_coverage(coverage_address_m, coverage_address_p, flag=True):

    coverage_address = coverage_address_m[:-6]+'.txt' #e.g. coverage_m.txt --> coverage.txt

    if flag==False:
        return coverage_address

    try:
        from itertools import izip as zip
    except ImportError: # will be 3.x series
        pass

    with open(coverage_address_m, 'r') as cm, open(coverage_address_p, 'r') as cp, open(coverage_address, 'w') as c:

        for mline, pline in zip(cm, cp):

            if mline=='' or pline=='':
                continue

            mtoks = mline.split()
            ptoks = pline.split()
            cov_sum = float(mtoks[4])+float(ptoks[4])

            newline = '\t'.join(mtoks[0:4])+'\t%f\n'%cov_sum

            c.write(newline)

    return coverage_address


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
        coverage_address_m = Default_Ref_Path + COV_fn_m
        COV_fn_filtered_m = COV_fn_m[:-4] + '_qt'+repr(qt[0])+'.txt'
        coverage_address_filtered_m = Default_Ref_Path + COV_fn_filtered_m
        #use filtered coverage.txt to make snps generated in high exp-level region
        FilterCoverage(coverage_address_m, coverage_address_filtered_m, line_cover_threshold) 
        
        Stat_p.set_qt(qt) 
        Stat_p.get_acc_cov_hd(COV_fn_p)
        Stat_p.set_acc_cov_qt()
        line_cover_threshold = Stat_p.acc_cover_qt[0]
        coverage_address_p = Default_Ref_Path + COV_fn_p
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
               flag=True,
               flag_m=True,
               flag_p=True):
    if flag==True:   
        GenTarget(ref_address, coverage_address_filtered_m, Num_SNP, tar_address_m, SNP_address_m, flag_m ) # Generate 2 random target2 (for paternal and maternal) sequence from ref_address by insering Num_SNP SNPs in the exonic positions
        GenTarget(ref_address, coverage_address_filtered_p, Num_SNP, tar_address_p, SNP_address_p, flag_p )
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

#when reads are generated from m and p, it's possible there're generated reads of same id
#to make all ids distinct so that the alignments from these reads are not considered as multiple alignment
#
#modified readFQ_address will replace original file
#logFile is readFQ_address.log (what read ids are modified)
def enforce_unique_readID(readFQ_address, flag=True):

    if flag==False:
        return

    src = readFQ_address
    dst = src + '.tmp'
    log = src + '.log'

    names = {}

    nLines=sum([1 for l in open(src,'r')]); T=nLines/100; p=0; q=0;
    with open(src, 'r') as sF, \
         open(dst, 'w') as dF, \
         open(log, 'w') as lF:

        for line in sF:
            p += 1
            if p>=T: p=0; q+=1; sys.stdout.write('\r'); sys.stdout.write('%d %% processed (enforce_unique_readID)'%q); sys.stdout.flush()

            if line[0]=='@':
                name = line.split()[0]
                if name in names:
                    names[name]+=1
                    name = name + '_%d'%names[name]
                    lF.write('duplicated readID modified as:'+name+'\n')
                else:
                    names[name]=0
                    name = name
                dF.write(name+'\n')
            else:
                dF.write(line)

    #pdb.set_trace()

    cmd = 'mv %s %s'%(dst, src)
    print('')
    run_cmd(cmd)

    return

def do_gen_read_alignment(ref_address,
                          readFQ_address,
                          case=6,
                          flag=True):
    if case==1:
        print('do_gen_read_alignment case1')
        aligner ='tophat'
        tophat_dir = '/tophat_out/'
        if flag==True:
            sam_address = Align(ref_address, readFQ_address, aligner, Default_Ref_Path+tophat_dir)
        else:
            sam_address = Default_Ref_Path + tophat_dir + '/accepted_hits.sam' #for test purpose
    
    elif case==6:
        print('do_gen_read_alignment case6')
        
        if flag==True:
            get_bam_from_gatk_best_practice()
            
        flder = Default_Ref_Path +'/data_GATK/2pass/'
        #name = 'split.bam' #Aligned.out.sam
        name = 'dedupped.bam'        
        #prepare sam file    
        if 'bam' in name:
            bam_address = flder + name
            sam_address = flder + name[:-4] + '.sam'
            if os.path.exists(sam_address)==False: #need to convert bam to sam
                SamtoolsProgram = SamPath + '/samtools'
                subprocess.call(SamtoolsProgram +  ' view -@ 20 -h ' + bam_address + ' > ' + sam_address, shell=True )
        else:
            sam_address = flder + name
            
    else:
        print('do_gen_read_alignment case?')
    return [sam_address]

def do_rsem(ref_address,
            readFQ_address,
            BED_address,
            gtf_address,
            EXON_fn,
            L,
            flag = True):
    if flag==True:
        RSEM_result_address = RSEM(ref_address, readFQ_address, BED_address, gtf_address)
        #RSEM_result_address = Default_Ref_Path + '/rsem/Chr15.isoforms.results'
        
        exon_address = BED2Exon(BED_address, EXON_fn)
        #exon_address = Default_Ref_Path + EXON_fn

        RSEM2Coverage(RSEM_result_address, exon_address,  L)
        count_rsem_address = Default_Ref_Path + '/count_rsem.txt'
        
    else: #do_rsem==False
        count_rsem_address = Default_Ref_Path + '/count_rsem.txt'
    return [count_rsem_address]
    
def do_gen_count(sam_address,
                 ref_address,                 
                 #coverage_address,#use this if use_rsem==False
                 count_rsem_address,
                 use_rsem=True,
                 flag=True): 
                     
    sam_fn = sam_address.split('/')[-1]
    sam_dir = sam_address[0:len(sam_address)-len(sam_fn)]
    count_fn = '/count_'+sam_fn[:-4]+'.txt'
        
    if flag==True:        
        if use_rsem==False:
            print("need to combine coverage_address_m.txt and coverage_address_p.txt")
            pdb.set_trace()
            #generate_count_file(sam_address, ref_address, coverage_address, count_fn) #en_debug_0817
        else:
            generate_count_file(sam_address, ref_address, count_rsem_address, count_fn)
    count_abs_address = sam_dir + count_fn
    
    return [count_abs_address]

def do_caller(sam_address,
              Threshold_num_reads,
              count_abs_address,
              flag=True):
    
    sam_fn = sam_address.split('/')[-1]
    ft = '.txt' #file type
    caller_op_addr = '/caller_output_'+sam_fn[:-4]+'_T'+repr(Threshold_num_reads)+ft
    caller_op_exception_addr = '/caller_output_exception_'+sam_fn[:-4]+'_T'+repr(Threshold_num_reads)+ft
    caller_op_snp_found_addr = '/caller_output_snp_found_'+sam_fn[:-4]+'_T'+repr(Threshold_num_reads)+ft    
    
    if flag==True:
        final_caller(count_abs_address,
                     caller_op_addr,
                     caller_op_exception_addr,
                     caller_op_snp_found_addr,
                     Threshold_num_reads)
                     
    return [caller_op_addr,
            caller_op_exception_addr,
            caller_op_snp_found_addr]

def do_snp_res_stat(sam_address,
                    caller_op_snp_found_addr,
                    Threshold_num_reads,
                    SNP_address_m,
                    SNP_address_p,
                    flag=True):
    
    caller_op_snp_found_fn = caller_op_snp_found_addr.split('/')[-1]
    snp_res_address = Default_Ref_Path + caller_op_snp_found_fn
    snp_res_stat_addr = Default_Ref_Path + caller_op_snp_found_fn[:-4]+'_[snp_res_stat].txt'
    
    snp_res_stat = []
    if flag==True:        
        snp_res_stat = do_snp_res_statistics(snp_res_address, SNP_address_m, SNP_address_p, snp_res_stat_addr)
    return [snp_res_stat_addr, snp_res_stat]

def do_filt_snp(sam_address,
                caller_op_snp_found_addr,
                count_abs_address,
                flag=True):
    
    caller_op_snp_found_addr = caller_op_snp_found_addr.split('/')[-1]
    
    if flag == True:
        #prepare input/output file names
        snp_res_address = Default_Ref_Path + caller_op_snp_found_addr
        #lam half filtering
        filt_snp_res_address2 = Default_Ref_Path + caller_op_snp_found_addr[:-4] + '_filt.txt'
        
        sam_fn = sam_address.split('/')[-1]
        sam_dir = sam_address[0:len(sam_address)-len(sam_fn)]
        count_altInfo_address = sam_dir + '/count_'+sam_fn[:-4]+'_altInfo.txt'
        
        #altCount modified
        filter_snp_lam_half_filt2(snp_res_address,
                                  filt_snp_res_address2,
                                  count_altInfo_address,
                                  count_abs_address)  
    else:
        filt_snp_res_address2 = Default_Ref_Path + caller_op_snp_found_addr[:-4] + '_filt.txt'
        
    return [filt_snp_res_address2]

def main(case = 6,
         para_comp = False,
         num_p = 1,
         Threshold_num_reads=1,
         samExp=0):

    #pdb.set_trace()
    
    print('Default data folder is %s'%Default_Ref_Path)
    
    Stat_m = sim_stat('/sim_stat_dmp_m.txt')
    Stat_p = sim_stat('/sim_stat_dmp_p.txt')
    
    ref_address = Default_Ref_Path + '/Chr15.fa'
    BED_address = Default_Ref_Path + '/hg19_chr15-UCSC.bed'
    
    EXP_fn_m = '/exp_m.txt' #for true exp
    EXP_fn_p = '/exp_p.txt'
    
    COV_fn_m = '/coverage_m.txt'
    COV_fn_p = '/coverage_p.txt'
    
    Num_SNP = 1000 #1000 #20 #for true snp
    
    tar_address_m = Default_Ref_Path + '/Tar_m.txt'
    tar_address_p = Default_Ref_Path + '/Tar_p.txt'

    SNP_address_m = Default_Ref_Path + '/SNP_m.txt' 
    SNP_address_p = Default_Ref_Path + '/SNP_p.txt'
    
    N=10000000 #10000000 #100000 #Number of Reads
    L=100       #Read Length
    error_rate = 0
    
    [exp_address_m, exp_address_p, BED_sorted_address] = do_gen_exp(BED_address,
                                                                    EXP_fn_m,
                                                                    EXP_fn_p,
                                                                    flag=False)

    if samExp==1:
        cmd = 'cp %s %s'%(exp_address_m, exp_address_p)
        run_cmd(cmd)
    
    [coverage_address_m, coverage_address_p] = do_gen_cov( exp_address_m,
                                                           exp_address_p,
                                                           BED_sorted_address,
                                                           COV_fn_m,
                                                           COV_fn_p,
                                                           Stat_m,
                                                           Stat_p,
                                                           flag=False)

    coverage_address = merge_coverage(coverage_address_m, coverage_address_p, flag=False) #note: filtered_m and filtered_p is only used for snp gen purpose; to get lambda info throughout genome, we use original cov info.
                                                           
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
               flag=False,
               flag_m=True,
               flag_p=False)#generate only snps from m

    [readFQ_address] = do_gen_reads( BED_address,
                                     tar_address_m,
                                     tar_address_p,
                                     exp_address_m,
                                     exp_address_p,
                                     N,
                                     L,
                                     error_rate,
                                     flag=False)

    #pdb.set_trace()

    enforce_unique_readID(readFQ_address, flag=False)
    
    [sam_address] = do_gen_read_alignment(ref_address,
                               readFQ_address,
                               case=6,
                               flag=False)
                          
    
    #pdb.set_trace()
    
    gtf_address = Default_Ref_Path + '/hg19_chr15-UCSC.gtf'
    EXON_fn='/exon.txt'     
    [count_rsem_address] = do_rsem( ref_address,
                                    readFQ_address,
                                    BED_address,
                                    gtf_address,
                                    EXON_fn,
                                    L,
                                    flag = False)
    
    #pdb.set_trace()
    if para_comp==False:
        sam_address = '/data1/shunfu1/SNPCalling/data/rsem/Chr15.genome.sorted_n.sam'
        print('test using rsem sam: %s'%sam_address)
        pdb.set_trace()
        [count_abs_address] = do_gen_count(  sam_address,
                                             ref_address,                 
                                             #coverage_address,#use this if use_rsem==False
                                             count_rsem_address,
                                             use_rsem=True,
                                             flag=True)
        
        pdb.set_trace()
        [caller_op_addr, caller_op_exception_addr, caller_op_snp_found_addr] = do_caller( sam_address,
                                                                                          Threshold_num_reads,
                                                                                          count_abs_address,
                                                                                          flag=True)
        [snp_res_stat_addr, snp_res_stat] = do_snp_res_stat(sam_address,
                                                            caller_op_snp_found_addr,
                                                            Threshold_num_reads,
                                                            SNP_address_m,
                                                            SNP_address_p,
                                                            flag=False)
                                                            
        #pdb.set_trace()
        [filt_snp_res_address2] = do_filt_snp(  sam_address,
                                                caller_op_snp_found_addr,
                                                count_abs_address,
                                                flag=False)
                                                
        [snp_res_stat_addr2, snp_res_stat2] = do_snp_res_stat(  sam_address,
                                                                filt_snp_res_address2,
                                                                Threshold_num_reads,
                                                                SNP_address_m,
                                                                SNP_address_p,
                                                                flag=False)
    else: #para

        #ref_address = '/data1/shunfu1/SNPCalling/data/Chr15.fa'
        #covAddress = '/data1/shunfu1/SNPCalling/data/count_rsem.txt'
        #sam_address = '/data1/shunfu1/SNPCalling/data/data_GATK/2pass/dedupped.sam'
        #output_address = '/data1/shunfu1/SNPCalling/data/' #not used yet

        #pdb.set_trace()

        testArg = '-r %s -c %s -s %s -O %s -T 1 -p 20'%(ref_address, count_rsem_address, sam_address, 'output_address_N/A')
        cmd = 'python batch_run_parallel_modi.py %s'%testArg
        #run_cmd(cmd)
                                                            
    #pdb.set_trace()
    
    return

def get_bam_from_gatk_best_practice():
    
    Ref_dir=Default_Proj_Path
    #"/home/sreeramkannan/singleCell/SNP_Calling_Summer15/"
    GATK_addr=Ref_dir+"/tools/GATK/GenomeAnalysisTK.jar"
    #"/home/sreeramkannan/singleCell/SNP_Calling_Summer15/tools/GATK/GenomeAnalysisTK.jar"
    Picard_addr=Ref_dir+"/tools/picard-tools-1.138/picard.jar"
    #"/home/sreeramkannan/singleCell/SNP_Calling_Summer15/tools/picard-tools-1.138/picard.jar"
    
    data_GATK_addr=Default_Ref_Path+"/data_GATK/"    
    cmd = "mkdir -p "+data_GATK_addr
    subprocess.call(cmd, shell=True)
    
    """
    # for step 7
    data_GATK_out_addr=data_GATK_addr+"/GATK_out/"
    cmd = "mkdir -p "+data_GATK_out_addr
    subprocess.call(cmd, shell=True) 
    """
    
    fa_addr=Default_Ref_Path+"/Chr15.fa"
    Genome_addr = data_GATK_addr+"/Chr15.fa"
    cmd = "cp "+fa_addr+" "+ Genome_addr
    subprocess.call(cmd, shell=True)
    
    Read_addr=Default_Ref_Path+"/Tar_read_l100.fastq"
    Genome_dict_addr=data_GATK_addr+"/Chr15.dict"
    
    #1. mapping to reference
    genomeDir=data_GATK_addr+"/Chr15/"
    cmd = "mkdir -p "+genomeDir
    subprocess.call(cmd, shell=True)
    cmd =   "STAR --runMode genomeGenerate "+\
            "--genomeDir "+genomeDir+" "+\
            "--genomeFastaFiles "+Genome_addr+\
            " --runThreadN 20"
    subprocess.call(cmd, shell=True)
    
    #2.
    runDir=data_GATK_addr+"/1pass/"
    cmd = "mkdir -p "+runDir
    subprocess.call(cmd, shell=True)
    #cmd = "cd "+runDir
    #subprocess.call(cmd, shell=True)
    cmd =   "STAR --genomeDir "+genomeDir+\
            " --readFilesIn "+Read_addr+\
            " --outFileNamePrefix "+runDir+\
            " --runThreadN 20"
    subprocess.call(cmd, shell=True)
    
    #3. mapping to reference (2pass)
    genomeDir_2pass=data_GATK_addr+"/Chr15_2pass/"
    cmd = "mkdir -p "+genomeDir_2pass
    subprocess.call(cmd, shell=True)
    cmd =   "STAR --runMode genomeGenerate "+\
            "--genomeDir "+genomeDir_2pass+" "\
            "--genomeFastaFiles "+ Genome_addr+" "\
            "--sjdbFileChrStartEnd "+ runDir+"/SJ.out.tab "+\
            "--sjdbOverhang 75"+\
            " --runThreadN 20"
    subprocess.call(cmd, shell=True)
    
    #4.
    runDir_2pass=data_GATK_addr+"/2pass/"
    cmd = "mkdir -p "+runDir_2pass
    subprocess.call(cmd, shell=True)
    #cmd = "cd "+runDir_2pass
    #subprocess.call(cmd, shell=True)
    cmd =   "STAR --genomeDir "+genomeDir_2pass+" "\
            "--readFilesIn "+Read_addr+\
            " --outFileNamePrefix "+runDir_2pass+\
            " --runThreadN 20"
    subprocess.call(cmd, shell=True)
    
    #5. use Picard to: 
    cmd ="java -jar "+Picard_addr+" AddOrReplaceReadGroups "+\
         "I="+runDir_2pass+"/Aligned.out.sam "+\
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
    
    """
    #7. Variant Caller

    Res_addr=data_GATK_out_addr+"/raw_variants.vcf"
    
    #generate .dict file
    #java -jar $Picard_addr CreateSequenceDictionary \
    #R=$Genome_addr \
    #O=$Genome_dict_addr

    #generate .fai file
    #samtools faidx $Genome_addr
    #
    java -jar $GATK_addr \
    -T   HaplotypeCaller \
    -R   $Genome_addr \
    -I   $runDir_2pass/split.bam \
    -dontUseSoftClippedBases \
    -stand_emit_conf 20 \
    -stand_call_conf 20 \
    
    """
    
    return
        
if __name__ == "__main__":
    
    #pdb.set_trace()
    
    main()