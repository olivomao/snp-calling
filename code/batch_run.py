# -*- coding: utf-8 -*-
"""
Created on Thu Aug  6 14:48:40 2015

@author: olivo
"""
from Address import *
from Synthesis import *
from ReadProcess import *
from count_read_lambda import *
from final_caller19thjuly_m import *

from debug_MACRO import *
from sim_stat import *
from find_pos import *
import pdb
import subprocess

# clear folder
pdb.set_trace()
print('Default data folder is %s'%Default_Ref_Path)
c = int(input('clear default data folter? 1-yes 0-no:'))
if c==1:
    subprocess.call('rm -r '+ Default_Ref_Path + '/*', shell = True)
    if isServerSide == True:
        pdb.set_trace()
        print('Please make sure the following folder/data exist:' + Default_Proj_Path+'/data_genome/')
    subprocess.call('cp '+Default_Proj_Path+'/data_genome/* '+Default_Ref_Path, shell = True)

# configuration

Stat = sim_stat('/sim_stat_dmp.txt')

## ---------- synthesis data
ref_address = Default_Ref_Path + '/Chr15.fa'
BED_address = Default_Ref_Path + '/hg19_chr15-UCSC.bed'

EXP_fn = '/exp.txt' #for true exp
COV_fn = '/coverage.txt'

Num_SNP = 1000 #20 #for true snp
tar_address_m = Default_Ref_Path + '/Tar_m.txt'
tar_address_p = Default_Ref_Path + '/Tar_p.txt'
SNP_address_m = Default_Ref_Path + '/SNP_m.txt' 
SNP_address_p = Default_Ref_Path + '/SNP_p.txt'

N=10000000    #Number of Reads
L=100       #Read Length
error_rate = 0

#[BED_sorted_address, exp_address] = BED2ExpressionLevel(BED_address, EXP_fn) # generate a random expression level file
exp_address = Default_Ref_Path + EXP_fn #for test purpose
BED_sorted_address = Default_Ref_Path + 'hg19_chr15-UCSC-sorted.bed'
#coverage_address = ExpressionLevel2Coverage(BED_sorted_address, exp_address, COV_fn, Stat) # Find the exonic parts with coverage not less than a ceryain value
coverage_address = Default_Ref_Path + COV_fn #for test purpose
#pdb.set_trace()

"""
if en_debug_0814 == 1:
    
    
    #to generate data for checking by test_SNP_covered_by_reads()
    
    #pdb.set_trace()
    
    qt = [90, 80, 70, 60, 50]
    Stat.set_qt(qt) 
    Stat.get_acc_cov_hd(COV_fn)
    Stat.set_acc_cov_qt()
    Stat.dump_stats_Synthesis()
    acc_cover_qt = Stat.acc_cover_qt
    
    #pdb.set_trace()
    
    for i in range(len(qt)):
        
        #pdb.set_trace()
        rel_case_dir = '/case_qt' + repr(qt[i]) + '_N' + repr(N) + '/'
        case_dir = Default_Ref_Path + rel_case_dir
        subprocess.call('mkdir -p '+ case_dir, shell = True)
        
        line_cover_threshold = acc_cover_qt[i]
        COV_fn_filtered = '/coverage_qt'+repr(qt[i])+'.txt'
        coverage_address_filtered = Default_Ref_Path + COV_fn_filtered

        #use filtered coverage.txt to make snps generated in high exp-level region
        FilterCoverage(coverage_address, coverage_address_filtered, line_cover_threshold) 
    
        GenTarget(ref_address, coverage_address_filtered, Num_SNP, tar_address_m, SNP_address_m ) # Generate 2 random target2 (for paternal and maternal) sequence from ref_address by insering Num_SNP SNPs in the exonic positions
        GenTarget(ref_address, coverage_address_filtered, Num_SNP, tar_address_p, SNP_address_p )
        
        [readBED_address_m, readFA_address_m, readFQ_address_m] = ReadGeneration(tar_address_m, BED_address, exp_address,  N, L, error_rate)
        #readFQ_address_m = Default_Ref_Path + '/Tar_m_read_l'+ repr(L) +'.fastq' #for test purpose
        [readBED_address_p, readFA_address_p, readFQ_address_p] = ReadGeneration(tar_address_p, BED_address, exp_address,  N, L, error_rate)
        #readFQ_address_p = Default_Ref_Path + '/Tar_p_read_l'+ repr(L) +'.fastq' #for test purpose
        
        readFQ_address = Default_Ref_Path + '/Tar_read_l'+ repr(L) +'.fastq'  #for test purpose
        subprocess.call('cat '  + readFQ_address_m + ' ' + readFQ_address_p + ' > ' + readFQ_address, shell=True )
        
        #move files
        subprocess.call('cp ' + BED_address + ' ' + case_dir, shell=True)
        subprocess.call('mv ' + coverage_address_filtered + ' ' + case_dir, shell=True)
        subprocess.call('mv ' + tar_address_m + ' ' + case_dir, shell=True)
        subprocess.call('mv ' + SNP_address_m + ' ' + case_dir, shell=True)        
        subprocess.call('mv ' + tar_address_p + ' ' + case_dir, shell=True)
        subprocess.call('mv ' + SNP_address_p + ' ' + case_dir, shell=True)
        subprocess.call('mv ' + readFQ_address_m + ' ' + case_dir, shell=True)
        subprocess.call('mv ' + readFA_address_m + ' ' + case_dir, shell=True)
        subprocess.call('mv ' + readBED_address_m + ' ' + case_dir, shell=True)
        subprocess.call('mv ' + readFQ_address_p + ' ' + case_dir, shell=True)
        subprocess.call('mv ' + readFA_address_p + ' ' + case_dir, shell=True)
        subprocess.call('mv ' + readBED_address_p + ' ' + case_dir, shell=True)
        subprocess.call('mv ' + readFQ_address + ' ' + case_dir, shell=True)                
       
    pdb.set_trace()
    

else:
"""

qt = [50]
#Stat.set_qt(qt) 
#Stat.get_acc_cov_hd(COV_fn)
#Stat.set_acc_cov_qt()
#Stat.dump_stats_Synthesis()
#acc_cover_qt = Stat.acc_cover_qt
    
#line_cover_threshold = acc_cover_qt[0]
COV_fn_filtered = '/coverage_qt'+repr(qt[0])+'.txt'
coverage_address_filtered = Default_Ref_Path + COV_fn_filtered

#use filtered coverage.txt to make snps generated in high exp-level region
#FilterCoverage(coverage_address, coverage_address_filtered, line_cover_threshold) 
    
#GenTarget(ref_address, coverage_address_filtered, Num_SNP, tar_address_m, SNP_address_m ) # Generate 2 random target2 (for paternal and maternal) sequence from ref_address by insering Num_SNP SNPs in the exonic positions
#GenTarget(ref_address, coverage_address_filtered, Num_SNP, tar_address_p, SNP_address_p )
        
#[readBED_address_m, readFA_address_m, readFQ_address_m] = ReadGeneration(tar_address_m, BED_address, exp_address,  N, L, error_rate)
readFQ_address_m = Default_Ref_Path + '/Tar_m_read_l'+ repr(L) +'.fastq' #for test purpose
#[readBED_address_p, readFA_address_p, readFQ_address_p] = ReadGeneration(tar_address_p, BED_address, exp_address,  N, L, error_rate)
readFQ_address_p = Default_Ref_Path + '/Tar_p_read_l'+ repr(L) +'.fastq' #for test purpose
        
readFQ_address = Default_Ref_Path + '/Tar_read_l'+ repr(L) +'.fastq'  #for test purpose
#subprocess.call('cat '  + readFQ_address_m + ' ' + readFQ_address_p + ' > ' + readFQ_address, shell=True )

#pdb.set_trace()
## ---------- read alignment & abundance estimation
aligner ='tophat'
gtf_address = Default_Ref_Path + '/hg19_chr15-UCSC.gtf'
tophat_dir = '/tophat_out'
EXON_fn='/exon.txt'

#en_debug_0811 sam_address = Align(ref_address, readFQ_address_m + ' ' + readFQ_address_p, aligner, Default_Ref_Path+tophat_dir)
#sam_address = Align(ref_address, readFQ_address, aligner, Default_Ref_Path+tophat_dir)
sam_address = Default_Ref_Path + tophat_dir + '/accepted_hits.sam' #for test purpose

#RSEM_result_address = RSEM(ref_address, readFQ_address, BED_address, gtf_address)
#RSEM_result_address = Default_Ref_Path + '/rsem/Chr15.isoforms.results'

#exon_address = BED2Exon(BED_address, EXON_fn)
#exon_address = Default_Ref_Path + EXON_fn

#RSEM2Coverage(RSEM_result_address, exon_address,  L)
#count_rsem_address = Default_Ref_Path + '/count_rsem.txt' #for test purpose

## ---------- generate count file
count_fn = '/count.txt'
#generate_count_file(sam_address, ref_address, coverage_address, count_fn) #en_debug_0817
count_rel_address = tophat_dir + count_fn #for test purpose

## ---------- caller
pdb.set_trace()
final_caller(count_rel_address, '/caller_output.txt')
