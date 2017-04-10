# -*- coding: utf-8 -*-
"""
Created on Wed Sep 23 22:09:52 2015

@author: olivo
"""

import pdb
import csv
import re
import sys
from progress.bar import Bar

def get_snp_dic(snp_addr):
    
    print('debug at get_snp_dic')
    pdb.set_trace()
    
    snp_dic = {}
    
    with open(snp_addr) as snp_file:        
        reader = csv.reader(snp_file, delimiter='\t')
        for row in reader:
            rPos = int(row[0])
            rB = row[1].split()[0].upper()
            nB = row[3].split()[0].upper()
            if rPos in snp_dic.keys():
                snp_dic[rPos].append([rB, nB])
            else:
                snp_dic[rPos] = [[rB, nB]]
                
    return snp_dic
    
def get_ref_dna(ref_address):
    
    print('debug at get_ref_dna')    
    pdb.set_trace()
    
    ref_file=open(ref_address,'rU')
    REF=''
    #counter=0    
    for line in ref_file:
        #if counter==0:
        #    header=line
        #    counter=1
        Segment = re.search(r'([acgtnACGTN]+)',line)
        if Segment!=None: 
            Seg = Segment.group(1)
        if Seg == line[:len(line)-1]:
            REF=REF+Seg
            
    ref_file.close()
    
    return REF

def process_bed_line_info(line1): #bed_line_info = 

    #print('debug at process_bed_line_info')
    #pdb.set_trace()
    
    itms = line1.split('\t')
    rid = itms[3]
    direction = itms[5]
    stt_pos = int(itms[6])
    
    num_block = int(itms[9])
    blk_len_s = itms[10].split(',')
    blk_len_s = [int(blk_len_s[i]) for i in range(num_block)]
    blk_stt_pos_s = itms[11].split(',')
    blk_stt_pos_s = [int(blk_stt_pos_s[i]) for i in range(num_block)]    
    
    return [rid, direction, stt_pos, num_block, blk_len_s, blk_stt_pos_s]

def process_fq_line_info(line21, line22, line23, line24): #fq_line_info =

    #print('debug at process_fq_line_info')
    #pdb.set_trace()
    
    rid = line21[1:len(line21)-1]
    generated_read = line22[0:100]
    
    return [rid, generated_read]

def comp_base(b):
    b = b.upper()
    
    res = ''
    if b == 'A':
        res = 'T'
    elif b == 'T':
        res = 'A'
    elif b == 'C':
        res = 'G'
    elif b == 'G':
        res = 'C'
    else:
        print('exception at comp_base')
        pdb.set_trace()
        
    return res

def check_reads(snp_addr, read_bed_addr, read_fq_addr, ref_dna_addr):
    
    print('debug at check_reads')
    pdb.set_trace()
    
    read_err_dic = {}
    
    snp_dic = get_snp_dic(snp_addr)
    
    ref_dna = get_ref_dna(ref_dna_addr)
    
    num_lines = sum(1 for line in open(read_bed_addr))
    bar = Bar('Progress of check_reads', max=num_lines)

    pdb.set_trace()
    
    read_bed_file = open(read_bed_addr, 'r')
    read_fq_file = open(read_fq_addr, 'r')
    
    while True:
        #pdb.set_trace()
        
        line1 = read_bed_file.readline()
        #sys.stdout.write(line1)
        
        line21 = read_fq_file.readline()
        #sys.stdout.write(line21)
        line22 = read_fq_file.readline()
        #sys.stdout.write(line22)
        line23 = read_fq_file.readline()
        #sys.stdout.write(line23)
        line24 = read_fq_file.readline()
        #sys.stdout.write(line24)
        
        if not line1 and not line21:
            break        
        
        bed_line_info = process_bed_line_info(line1)
        
        fq_line_info = process_fq_line_info(line21, line22, line23, line24)
        
        bar.next()
        #pdb.set_trace()
        
        if bed_line_info[0]==fq_line_info[0]:
            
            #if bed_line_info[0] == 'ENST00000526079.1_e_2_chr15_23447189':
            #    pdb.set_trace()
                
            read_base_counter = 0
            direction = bed_line_info[1]
            genome_stt_pos = bed_line_info[2]
            num_block = bed_line_info[3]
            blk_len_s = bed_line_info[4]
            blk_stt_pos_s = bed_line_info[5]
            generated_read = fq_line_info[1]
            
            for ith_blk in range(num_block):
                for blk_pos in range(blk_len_s[ith_blk]):
                    genome_pos = genome_stt_pos + blk_stt_pos_s[ith_blk] + blk_pos
                    #if genome_pos == 28765945:
                    #    pdb.set_trace()
                    ref_base = ref_dna[genome_pos].upper()
                    read_base = generated_read[read_base_counter].upper()
                    
                    #if bed_line_info[0] == 'ENST00000526079.1_e_2_chr15_23447189' and read_base_counter==45:
                    #    pdb.set_trace()
                        
                    if direction == '-':
                        read_base = comp_base(generated_read[99-read_base_counter].upper())
                    if ref_base != read_base:
                        if genome_pos not in snp_dic.keys():
                            print('debug: check new read err')                        
                            pdb.set_trace()                        
                            rerr = read_err_dic.setdefault(genome_pos, [])# count = counts.setdefault(genome_pos, [])
                            #rerr.append([bed_line_info, fq_line_info, ith_blk, blk_pos])
                            rerr.append('err') #([bed_line_info, fq_line_info, ith_blk, blk_pos])
                        else:
                            print('debug: check new read snp')                        
                            pdb.set_trace()  
                            rerr = read_err_dic.setdefault(genome_pos, [])# count = counts.setdefault(genome_pos, [])
                            rerr.append('snp') #([bed_line_info, fq_line_info, ith_blk, blk_pos, 'snp'])
                    #if bed_line_info[0] == 'ENST00000526079.1_e_2_chr15_23447189':
                    #    print('%d: ref_base=%s ref_dna=%s'%(read_base_counter, ref_base, read_base))
                    read_base_counter = read_base_counter + 1
        else: #exception
            print('read_bed_file and read_fq_file not compatible (diff read id)')
            pdb.set_trace()
        
        #if not line1 and not line21:
        #    break
    
    read_bed_file.close()
    read_fq_file.close()

    bar.finish()
    pdb.set_trace()    
    
    return read_err_dic

if __name__ == "__main__":
    
    working_dir = '/home/olivo/Downloads/0907_count_alt_mapping/bkp/data_0814/'
    
    snp_fn = '/SNP_p.txt'
    snp_addr = working_dir + snp_fn
    
    read_bed_fn = '/Tar_p_read_l100.bed'
    read_bed_addr = working_dir + read_bed_fn
    
    read_fq_fn = '/Tar_p_read_l100.fastq'
    read_fq_addr = working_dir + read_fq_fn
    
    ref_dna_fn = '/Chr15.fa'
    ref_dna_addr = working_dir + ref_dna_fn
    
    read_err_dic = check_reads(snp_addr, read_bed_addr, read_fq_addr, ref_dna_addr)
    
    print('ready to exit')
    pdb.set_trace()
    pdb.set_trace()