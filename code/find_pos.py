# -*- coding: utf-8 -*-
"""
Created on Tue Aug 11 04:24:57 2015

@author: olivo
"""

import math
import csv
import operator
import pickle
from Address import *
import pdb
from debug_MACRO import *
from progress.bar import Bar
import subprocess
from sim_stat import *

def process_row_cov(pos_idx, row, OUT_file):
    #from Address import *
# row example:
# ['11849',   '66621617',        '66621737',        '26267475',       '738443411.775']
    e_stt = int(row[1])
    e_stp = int(row[2])
    
    res_covered = 0
    
    if pos_idx >= e_stt and pos_idx <= e_stp:
        
        res_covered = 1        
        
        OUT_file.write('find pos %d in this block:\n'%pos_idx)
        OUT_file.write('\t'.join(row) + '\n')
        OUT_file.write('\n')
            
        print('\nfind pos %d in this block:'%pos_idx)
        print('\t'.join(row) + '')
        print('')
    
    return res_covered

def process_row(pos_idx, row, OUT_file):
    #from Address import *
# row example:
# ['chr15', '20192908', '20193370', 'ENST00000558565.2', '0', '-', '20192908', '20193370', '0', '2', '313,46,', '0,416,']
    e_stt = int(row[1])
    e_stp = int(row[2])
    rid = row[3] #.split(',')[0]
    blocks_stt = row[11].split(',')
    blocks_sz = row[10].split(',')    
    n_blocks = int(row[9]) #len(blocks_stt)
    
    res_covered = 0
    
    for i in range(n_blocks):
        block_stt = e_stt + int(blocks_stt[i])
        block_stp = block_stt + int(blocks_sz[i]) - 1
        
        if pos_idx >= block_stt and pos_idx <= block_stp:
            
            res_covered += 1            
            
            OUT_file.write('find pos %d in this block:\n'%pos_idx)
            OUT_file.write('\t'.join(row) + '\n')
            OUT_file.write('rid=%s  e_stt~e_stp:%d~%d\n'%(rid, e_stt, e_stp))
            OUT_file.write('\tblock%d %d~%d\n'%(i+1, block_stt, block_stp))
            OUT_file.write('\n')
            
            print('\nfind pos %d in this block:'%pos_idx)
            print('\t'.join(row) + '')
            print('rid=%s  e_stt~e_stp:%d~%d'%(rid, e_stt, e_stp))
            print('\tblock%d %d~%d'%(i+1, block_stt, block_stp))
            print('')
    
    return res_covered
    
"""
objective: to find if a pos is covered by any block indicated in COV_fn
           if so, print the relevant line in COV_fn, detailed block stt & stp
"""
def find_pos_cov(Ref_Path,
                 rel_case_dir,                 
                 COV_fn,
                 pos_idx,
                 OUT_file):
    #from Address import *

    COV_address = Ref_Path + rel_case_dir + COV_fn
    
    #num_lines = sum(1 for line in open(COV_address))
    #bar = Bar('Processing', max=num_lines)
    
    with open(COV_address) as COV_file:
        
        reader = csv.reader(COV_file, delimiter='\t')
        
        OUT_file.write('==========%s\n'%COV_fn)
        print('==========%s'%COV_fn)
        
        num_occur = 0
        
        for row in reader:
            #bar.next()
            #pdb.set_trace()
            tmp_occur = process_row_cov(pos_idx, row, OUT_file)
            if tmp_occur > 0:
                num_occur = num_occur + 1
        #bar.finish()
                
        OUT_file.write('num_occur = %d\n'%num_occur)
        print('num_occur = %d'%num_occur)
        
    return num_occur

"""
objective: to find if a pos is covered by any block indicated in BED_fn
           if so, print the relevant line in BED_fn, detailed block stt & stp
"""
def find_pos(Ref_Path,
             rel_case_dir,                 
             BED_fn,
             pos_idx,
             OUT_file):
    #from Address import *

    #num_lines = sum(1 for line in open(BED_fn))
    #bar = Bar('Processing', max=num_lines)
    
    with open(Ref_Path + rel_case_dir + BED_fn) as BED_file:
        
        reader = csv.reader(BED_file, delimiter='\t')
        
        OUT_file.write('==========%s\n'%BED_fn)
        print('==========%s'%BED_fn)
        
        num_occur = 0
        
        for row in reader:
            #bar.next()
            #pdb.set_trace()
            tmp_occur = process_row(pos_idx, row, OUT_file) 
            if tmp_occur > 0:
                num_occur = num_occur + 1
        #bar.finish()
        
        OUT_file.write('num_occur = %d\n'%num_occur)
        print('num_occur = %d'%num_occur)
    
    return num_occur
    
def find_pos_0(Ref_Path,
               rel_case_dir,
               SNP_fn,
               rel_res_dir,
               COV_fn,
               BED_fn,
               Stat):
    #from Address import *
    
    #pdb.set_trace()
    
    res_dir = Ref_Path+rel_case_dir+rel_res_dir#extract COV_fn without '.txt'
    
    subprocess.call( 'mkdir -p ' + res_dir,  shell=True )
    print('Folder loc %s'%res_dir)
    c = int(input('clear data in this folter? 1-yes 0-no:'))
    if c==1:
        subprocess.call( 'rm -r '+ res_dir + '/*.*',  shell=True )
                   
    with open(Ref_Path+rel_case_dir+SNP_fn) as SNP_file:
        
        #COV_fn = ['coverage.txt']        
        
        #BED_fn = ['hg19_chr15-UCSC.bed',
        #          'Tar_m_read_l100.bed',
        #          'Tar_p_read_l100.bed']
        
        #pdb.set_trace()
        
        reader = csv.reader(SNP_file, delimiter='\t')
        
        for row in reader:
            print(row)
            pos_idx = int(row[0])
            OUT_fn = res_dir + '/find_pos_' + str(pos_idx) + '.txt'
            OUT_file = open(OUT_fn, 'w+')
            
            #pdb.set_trace()
            for i in range(len(COV_fn)):
                find_pos_cov(Ref_Path,
                     rel_case_dir,                 
                     COV_fn[i],
                     pos_idx,
                     OUT_file)
                
            for i in range(len(BED_fn)):            
                tmp_res = find_pos(Ref_Path,
                                   rel_case_dir,                 
                                   BED_fn[i],
                                   pos_idx,
                                   OUT_file)
                #pdb.set_trace()
                #if tmp_res>0 and 'read' in BED_fn[i] and 'SNP_m' in SNP_fn:
                #    Stat.num_snp_m_covered_by_reads += 1
                #    break
                #if tmp_res>0 and 'read' in BED_fn[i] and 'SNP_p' in SNP_fn:
                #    Stat.num_snp_p_covered_by_reads += 1
                #    break
            
            OUT_file.close()
            #pdb.set_trace()
            
    return
    
#assume relevant files are already there
def test_SNP_covered_by_reads():
    pdb.set_trace()
    
    #Ref_Path = '/home/olivo/Desktop/SNP-Calling-Summer15/data_0814_stat_check_snp_covered_by_reads_N10000/'
    Ref_Path = '/home/sreeramkannan/singleCell/SNP_Calling_Summer15/data_0814/'
    
    vqt = [90]
    N = 10000

    for qt in vqt:            
        #rel_case_dir = '/case_qt' + repr(qt) + '_N' + repr(N) + '/'
        rel_case_dir = ''
        #COV_fn = '/coverage_qt'+repr(qt)+'.txt'
        COV_fn = '/coverage.txt'
        
        SNP_m_fn= 'snp_pos_for_analysis_cross_check_case7.txt'
        #SNP_p_fn= 'SNP_p.txt'
        COV_fn = [COV_fn]
        BED_fn = ['hg19_chr15-UCSC.bed',
                  'Tar_m_read_l100.bed',
                  'Tar_p_read_l100.bed']    
        
        Stat = sim_stat()
        Stat.set_dmp_addr(Ref_Path + rel_case_dir + '/sim_stat_dmp.txt')    
                  
        Stat.num_snp_m_covered_by_reads = 0
        Stat.num_snp_p_covered_by_reads = 0
        
        find_pos_0(Ref_Path,
                   rel_case_dir,
                   SNP_m_fn,
                   '/check_snp_pos_cross_check_case7/m/',
                   COV_fn,
                   BED_fn,
                   Stat)
                   
        #find_pos_0(Ref_Path,
        #           rel_case_dir,
        #           SNP_p_fn,
        #           '/check_snp_pos/p/',
        #           COV_fn,
        #           BED_fn,
        #           Stat)
                   
        #pdb.set_trace()
        
        Stat.dmp_scaler(Stat.num_snp_m_covered_by_reads,
                     'num of m SNPs covered by reads (N='+repr(N)+', qt='+repr(qt)+')')
                     
        #Stat.dmp_scaler(Stat.num_snp_p_covered_by_reads,
        #             'num of p SNPs covered by reads (N='+repr(N)+', qt='+repr(qt)+')')
        
        #pdb.set_trace()
    
    return
    
def test_SNP_coverage_statistics():
    
    pdb.set_trace()
    
    Ref_Path = '/home/sreeramkannan/singleCell/SNP_Calling_Summer15/data_0814/'
    
    vqt = [90]
    N = 100000

    for qt in vqt:            
        #rel_case_dir = '/case_qt' + repr(qt) + '_N' + repr(N) + '/'
        rel_case_dir = ''
        #COV_fn = '/coverage_qt'+repr(qt)+'.txt'
        COV_fn = '/coverage_qt90.txt'
        
        SNP_m_fn= 'snp_pos_for_analysis.txt'        
        COV_fn = [COV_fn]
        BED_fn = ['hg19_chr15-UCSC.bed',
                  'Tar_m_read_l100.bed',
                  'Tar_p_read_l100.bed']    
        
        Stat = sim_stat()
        Stat.set_dmp_addr(Ref_Path + rel_case_dir + '/sim_stat_dmp.txt')    
                  
        Stat.num_snp_m_covered_by_reads = 0
        
        find_pos_0(Ref_Path,
                   rel_case_dir,
                   SNP_m_fn,
                   '/check_snp_pos/',
                   COV_fn,
                   BED_fn,
                   Stat)
                   
        #pdb.set_trace()
        
        Stat.dmp_scaler(Stat.num_snp_m_covered_by_reads,
                     'num of m SNPs covered by reads (N='+repr(N)+', qt='+repr(qt)+')')
        
        #pdb.set_trace()
    
    return

def test_SNP_covered_by_reads_0928():
    
    pdb.set_trace()
    
    #Ref_Path = '/home/sreeramkannan/singleCell/SNP_Calling_Summer15/data_0814/'
    #Ref_Path = '/home/sreeramkannan/singleCell/SNP_Calling_Summer15/data_SNP1k_Reads10M/'
    Ref_Path = '/home/sreeramkannan/singleCell/SNP_Calling_Summer15/data_SNP20_Reads100K_diffMPExp/'
    
    rel_case_dir = ''
    COV_fn_m = '/coverage_m.txt'
    COV_fn_p = '/coverage_p.txt'
    
    SNP_fn= 'snp_pos_for_analysis.txt'
    
    COV_fn = [COV_fn_m,
              COV_fn_p]
    BED_fn = ['hg19_chr15-UCSC.bed',
              'Tar_m_read_l100.bed',
              'Tar_p_read_l100.bed']    
    
    Stat = sim_stat()
    Stat.set_dmp_addr(Ref_Path + rel_case_dir + '/sim_stat_dmp.txt')    
              
    Stat.num_snp_m_covered_by_reads = 0
    Stat.num_snp_p_covered_by_reads = 0
    
    find_pos_0(Ref_Path,
               rel_case_dir,
               SNP_fn,
               '/check_snp_pos_1124/',
               COV_fn,
               BED_fn,
               Stat)
    
    return

if __name__ == "__main__":
    
    print('test_SNP_covered_by_reads')
    #test_SNP_covered_by_reads()
    test_SNP_covered_by_reads_0928() #mis-detection analysis after alt mapping modi for count_read_lambda
    
    #print('test_SNP_coverage_statistics')
    #test_SNP_coverage_statistics()
    pdb.set_trace()
