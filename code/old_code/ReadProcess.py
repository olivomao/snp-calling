"""
- 10/24: RSEM2Coverage is modified:
         an exon genome pos may occur >= 2 times, merge the same genome pos by using dic of ExonPosExp
         (key, val)=(genome/exon pos, accumulated expression level/base_count)
"""
from Address import *
#import Bowtie
#import Rules

import os
import subprocess
import re
import time
import sys
from operator import itemgetter, attrgetter
import random
from numpy import array
import math
#from matplotlib import pyplot as plt
import numpy as np
import heapq

import datetime
import pdb
from debug_MACRO import *

 #--------------------------------------------------------------------------------------------------------

def Align(ref_address, readFQ_address, aligner,out_address, num_p=1):
    # aligner can be bowtie or bowtie2 or tophat
    #pdb.set_trace()
    b=re.search('([\S]+)/([\S]+)\.([\S]+)', ref_address)
    Path = b.group(1)
    
    # index the reference/transcriptome
    Index_address = Path + '/index/' + b.group(2) 
     
    if (aligner == 'bowtie2') or (aligner =='tophat'): 
        B_Path = Bowtie2Path
        if not os.path.isfile(Index_address + '.1.bt2'):
            #subprocess.call( 'mkdir ' + Path + '/index', shell=True)
            subprocess.call( 'mkdir -p ' + Path + '/index', shell=True)
            Index_address = Bowtie.index(B_Path, ref_address, aligner)
        else:
            print ' index already exists!'
    elif aligner == 'bowtie':
        B_Path = BowtiePath
        if not os.path.isfile(Index_address + '.1.ebwt'):
            subprocess.call( 'mkdir ' + Path + '/index', shell=True)
            Index_address = Bowtie.index(B_Path, ref_address, aligner)
        else:
            print ' index already exists!'
    else:
        print 'invalid aligner!'
         
    # align using Bowtie to Transcriptome and generate .sam format
    if (aligner == 'bowtie') or (aligner == 'bowtie2'):
        b=re.search('([\S]+)\.([\S]+)', readFQ_address)
        sam_address = b.group(1)+ '.sam'
        Bowtie.Bowtie(B_Path, Index_address, readFQ_address, sam_address, ' -a -v3 --sam ', aligner)
    elif aligner == 'tophat':
        
        if num_p > 1:
            print('Tophat to be run with num of threads=%d'%num_p)
            args = ' -p ' + repr(num_p) + ' -N 3 --read-edit-dist 4 '
        else:
            args = ' -N 3 --read-edit-dist 4 '
        Bowtie.TopHat(TopHatPath, Index_address, readFQ_address, out_address, args)
        #pdb.set_trace() #debug
        bam_address = out_address  + '/accepted_hits.bam'
        sam_address = out_address  + '/accepted_hits.sam'
        SamtoolsProgram = SamPath + '/samtools'
        subprocess.call(SamtoolsProgram +  ' view -h ' + bam_address + ' > ' + sam_address, shell=True )
        
    return sam_address

#--------------------------------------------------------------------------------------------------------

def BED2Exon(BED_address, EXON_fn='/exon.txt'):
    # given .bed file, generates exon.txt file
    # for each exon, provides a list of all transcript IDs the exon belongs to, as well as 
    # the transcript length before and after the start position of the exon
    
    #pdb.set_trace()
    
    b=re.search('([\S]+)/([\S]+)\.([\S]+)', BED_address)
    exon_address =  b.group(1) + EXON_fn #'/exon.txt'
    
    Dict = {}
    EXONS = [ ]
    
    counter =0
    bed_file = open(BED_address, 'rU')
    for line in bed_file:
        x = line.split()
        tr_ID = x[3]
        tr_start = int( x[1] )
        number_exon = int( x[9] )
        EXON_len_strng = x[10].split(',') 
        EXON_len = [ int(EXON_len_strng[i]) for i in range(number_exon) ]
        
        #pdb.set_trace()
        
        EXON_start_strng = x[11].split(',')
        EXON_start = [ int(EXON_start_strng[i]) for i in range(number_exon) ]
        
        #pdb.set_trace()
        
        for i in range(number_exon):
            exon_start = tr_start + int(EXON_start[i])
            exon_end = exon_start + int(EXON_len[i])
            if (exon_start,exon_end) in Dict.keys():
                EXONS [ Dict[ (exon_start,exon_end) ] ].append( ( tr_ID, sum(EXON_len[:i]), sum(EXON_len[i+1:]) ) )
            else:
                Dict[ (exon_start,exon_end) ] = counter
                EXONS.append( [ ( tr_ID, sum(EXON_len[:i]), sum(EXON_len[i+1:]) ) ] )
                counter += 1
        
        
    EXONS_sortet = sorted( Dict.keys() , key=lambda x: x[0])
    
    exon_file = open(exon_address,'w+')
    for exon in EXONS_sortet:
        exon_file.write( repr(exon[0]) + '\t' + repr(exon[1]) + '\t') 
        for tr in EXONS[Dict[exon]]:
            exon_file.write( tr[0] + ',' + repr(tr[1]) + ',' + repr(tr[2])  + '\t')
        #pdb.set_trace()
        exon_file.write('\n')
    exon_file.close()
    
    return exon_address  

def BED2Exon1(BED_address, exon_address):
    # given .bed file, generates exon.txt file
    # for each exon, provides a list of all transcript IDs the exon belongs to, as well as 
    # the transcript length before and after the start position of the exon
    
    #pdb.set_trace()
    
    #b=re.search('([\S]+)/([\S]+)\.([\S]+)', BED_address)
    #exon_address =  b.group(1) + EXON_fn #'/exon.txt'
    
    Dict = {}
    EXONS = [ ]
    
    counter =0
    bed_file = open(BED_address, 'rU')
    for line in bed_file:
        if line[0]=='#': continue
        
        x = line.split()
        tr_ID = x[3]
        tr_start = int( x[1] )
        number_exon = int( x[9] )
        EXON_len_strng = x[10].split(',') 
        EXON_len = [ int(EXON_len_strng[i]) for i in range(number_exon) ]
        
        #pdb.set_trace()
        
        EXON_start_strng = x[11].split(',')
        EXON_start = [ int(EXON_start_strng[i]) for i in range(number_exon) ]
        
        #pdb.set_trace()
        
        for i in range(number_exon):
            exon_start = tr_start + int(EXON_start[i])
            exon_end = exon_start + int(EXON_len[i])
            if (exon_start,exon_end) in Dict.keys():
                EXONS [ Dict[ (exon_start,exon_end) ] ].append( ( tr_ID, sum(EXON_len[:i]), sum(EXON_len[i+1:]) ) )
            else:
                Dict[ (exon_start,exon_end) ] = counter
                EXONS.append( [ ( tr_ID, sum(EXON_len[:i]), sum(EXON_len[i+1:]) ) ] )
                counter += 1
        
        
    EXONS_sortet = sorted( Dict.keys() , key=lambda x: x[0])
    
    exon_file = open(exon_address,'w+')
    for exon in EXONS_sortet:
        exon_file.write( repr(exon[0]) + '\t' + repr(exon[1]) + '\t') 
        for tr in EXONS[Dict[exon]]:
            exon_file.write( tr[0] + ',' + repr(tr[1]) + ',' + repr(tr[2])  + '\t')
        #pdb.set_trace()
        exon_file.write('\n')
    exon_file.close()
    
    return exon_address  
          
    #--------------------------------------------------------------------------------------------------------   

def RSEM(ref_address, readFQ_address, BED_address, gtf_address):
    pdb.set_trace()
    # estimate abundance level using RSEM (http://deweylab.biostat.wisc.edu/rsem/README.html)
    b=re.search('([\S]+)/([\S]+)\.([\S]+)', ref_address)
    path = b.group(1) + '/rsem'
    #subprocess.call( 'mkdir ' + path, shell=True )
    subprocess.call( 'mkdir -p ' + path, shell=True )

    name = path + '/' + b.group(2)
    rsem_index_command = rsemPath + '/rsem-prepare-reference'
    rsem_estim_command = rsemPath + '/rsem-calculate-expression'
   
    #pdb.set_trace()
    if en_debug == 0:    
        #subprocess.call( rsem_index_command + ' --gtf ' + gtf_address + ' ' + ref_address + ' '  + name, shell=True )
        subprocess.call( rsem_index_command + ' --star -p 20 --star-path %s --gtf '%STAR_path + gtf_address + ' ' + ref_address + ' '  + name, shell=True )
    else:
        #subprocess.call( rsem_index_command + ' --bowtie --gtf ' + gtf_address + ' ' + ref_address + ' '  + name, shell=True )
        subprocess.call( rsem_index_command + ' --star -p 20 --star-path %s --gtf '%STAR_path + gtf_address + ' ' + ref_address + ' '  + name, shell=True )
    #subprocess.call( rsem_estim_command + ' ' + readFQ_address+ ' ' + name + ' ' + name, shell=True )
    subprocess.call( rsem_estim_command + ' --star -p 20 --star-path %s '%STAR_path + readFQ_address+ ' ' + name + ' ' + name, shell=True )
    RSEM_result_address = name + '.isoforms.results'
    return RSEM_result_address 
    

#PE reads
#external
#python ReadProcess.py --RSEM2 --ref fa_file --r1 r1_file --r2 r2_file --gtf gtf_file
def RSEM2(args): #(ref_address, readFQ_address, BED_address, gtf_address):

    pdb.set_trace()

    ref_address = args[args.index('--ref')+1]
    r1 = args[args.index('--r1')+1]
    r2 = args[args.index('--r2')+1]
    gtf_address = args[args.index('--gtf')+1]
    # estimate abundance level using RSEM (http://deweylab.biostat.wisc.edu/rsem/README.html)
    b=re.search('([\S]+)/([\S]+)\.([\S]+)', ref_address)
    path = b.group(1) + '/rsem'
    #subprocess.call( 'mkdir ' + path, shell=True )
    subprocess.call( 'mkdir -p ' + path, shell=True )

    name = path + '/' + b.group(2)
    rsem_index_command = rsemPath + '/rsem-prepare-reference'
    rsem_estim_command = rsemPath + '/rsem-calculate-expression'
   
    print('skip rsem-prepare-reference if done'); pdb.set_trace()
    if en_debug == 0:    
        #subprocess.call( rsem_index_command + ' --gtf ' + gtf_address + ' ' + ref_address + ' '  + name, shell=True )
        pass; #subprocess.call( rsem_index_command + ' --star -p 20 --star-path %s --gtf '%STAR_path + gtf_address + ' ' + ref_address + ' '  + name, shell=True )
    else:
        #subprocess.call( rsem_index_command + ' --bowtie --gtf ' + gtf_address + ' ' + ref_address + ' '  + name, shell=True )
        pass; #subprocess.call( rsem_index_command + ' --star -p 20 --star-path %s --gtf '%STAR_path + gtf_address + ' ' + ref_address + ' '  + name, shell=True )
    #subprocess.call( rsem_estim_command + ' ' + readFQ_address+ ' ' + name + ' ' + name, shell=True )
    #subprocess.call( rsem_estim_command + ' --star -p 20 --star-path %s '%STAR_path + readFQ_address+ ' ' + name + ' ' + name, shell=True )
    cmd = '%s '%rsem_estim_command + \
          '--star -p 20 --star-path %s '%STAR_path + \
          '--paired-end %s %s '%(r1,r2) + \
          '%s '%name + \
          '%s/exp '%path #sample_name
    subprocess.call(cmd, shell=True)

    RSEM_result_address = name + '.isoforms.results'

    return RSEM_result_address 

#subsample top n reads from read file into dst
#usage
#python ReadProcess.py --SubsampleRead -r read_file -d dst -n n
def SubsampleRead(args): #read_file, dst, n):

    pdb.set_trace()

    read_file = args[args.index('-r')+1]
    dst = args[args.index('-d')+1]
    n = int(args[args.index('-n')+1])

    with open(dst, 'w') as f:

        i = 0
        with open(read_file, 'r') as r:

            for line in r:

                if line[0]=='#': continue 

                if line[0]=='@':
                    #pdb.set_trace()
                    i += 1
                    if i==n:break

                    f.write(line)
                    
                else:
                    f.write(line)

    return

     
def RSEM2Coverage(RSEM_result_address, exon_address,  L): 
    # Similar to ExpectedLevel2Coverage
    # computes the expected coverage depth based on estimated RPKM obtained from RSEM
    #pdb.set_trace()
    EXONS ={}
    exon_file = open(exon_address,'rU')
    for line in exon_file:
        x=line.split()
        exon_start = int(x[0])
        exon_end = int (x[1])
        exon = (exon_start, exon_end)
        z = []
        for i in range(2,len(x)):
            y=x[i].split(',')
            z.append( [y[0], int(y[1]), int(y[2])] )
        EXONS[exon] = z
    exon_file.close()
    EXONS_sortet = sorted( EXONS.keys() , key=lambda x: x[0])

    #pdb.set_trace()
    RSEM_file = open(RSEM_result_address, 'rU')
    RSEM_file.readline() # the first line is just header
    
    Avrg_Coverage = {}
    Tr = {}
    for line in RSEM_file:
        y = line.split()
        transcript_ID = y[0]
        transcript_expected_count = float( y[4] )
        transcript_len = int (y[2])
        Tr[transcript_ID] = [transcript_expected_count, transcript_len]

    RSEM_file.close() 
    
    #pdb.set_trace()
    b=re.search('([\S]+)/([\S]+)\.([\S]+)', exon_address)
    count_rsem_address = b.group(1) + '/count_rsem.txt'
    count_rsem_file = open(count_rsem_address,'w+')

    max_base_count = 0
    ExonPosExp = {}
    for exon in EXONS_sortet:
        for global_pos in range(exon[0],exon[1]):
            base_count = 0
            for transcript in EXONS[exon]:
                local_pos = (global_pos - exon[0]) + transcript[1]
                base_count += Transcriptome_contribution_to_exon(local_pos, Tr[transcript[0]], L)
            
            ExonPos = ExonPosExp.setdefault(global_pos, [])
            if ExonPos != []:
                #if base_count > 0:
                #    pdb.set_trace()
                ExonPos[0] = ExonPos[0] + base_count
            else:
                #pdb.set_trace()
                ExonPos.append(base_count)
            #count_rsem_file.write( repr(global_pos) + '\t' + repr(base_count) + '\n')
            if base_count > max_base_count:
                max_base_count = base_count

    #pdb.set_trace()
    ExonPosSorted = sorted( ExonPosExp.keys() ) # , key=lambda x: x[0])
    for gp in ExonPosSorted:
        e_l = float(ExonPosExp.get(gp,0)[0])
        count_rsem_file.write( repr(gp) + '\t' + repr(e_l) + '\n')
        
    count_rsem_file.close()
    print('max L0: %f'%max_base_count)
    
    return count_rsem_address

#use external address
def RSEM2Coverage1(RSEM_result_address, exon_address,  count_rsem_address, L): 
    # Similar to ExpectedLevel2Coverage
    # computes the expected coverage depth based on estimated RPKM obtained from RSEM
    #pdb.set_trace()
    EXONS ={}
    exon_file = open(exon_address,'rU')
    for line in exon_file:
        x=line.split()
        exon_start = int(x[0])
        exon_end = int (x[1])
        exon = (exon_start, exon_end)
        z = []
        for i in range(2,len(x)):
            y=x[i].split(',')
            z.append( [y[0], int(y[1]), int(y[2])] )
        EXONS[exon] = z
    exon_file.close()
    EXONS_sortet = sorted( EXONS.keys() , key=lambda x: x[0])

    #pdb.set_trace()
    RSEM_file = open(RSEM_result_address, 'rU')
    RSEM_file.readline() # the first line is just header
    
    Avrg_Coverage = {}
    Tr = {}
    for line in RSEM_file:
        y = line.split()
        transcript_ID = y[0]
        transcript_expected_count = float( y[4] )
        transcript_len = int (y[2])
        Tr[transcript_ID] = [transcript_expected_count, transcript_len]

    RSEM_file.close() 
    
    #pdb.set_trace()
    #b=re.search('([\S]+)/([\S]+)\.([\S]+)', exon_address)
    #count_rsem_address = b.group(1) + '/count_rsem.txt'
    count_rsem_file = open(count_rsem_address,'w+')

    max_base_count = 0
    ExonPosExp = {}
    for exon in EXONS_sortet:
        for global_pos in range(exon[0],exon[1]):
            base_count = 0
            for transcript in EXONS[exon]:
                local_pos = (global_pos - exon[0]) + transcript[1]
                base_count += Transcriptome_contribution_to_exon(local_pos, Tr[transcript[0]], L)
            
            ExonPos = ExonPosExp.setdefault(global_pos, [])
            if ExonPos != []:
                #if base_count > 0:
                #    pdb.set_trace()
                ExonPos[0] = ExonPos[0] + base_count
            else:
                #pdb.set_trace()
                ExonPos.append(base_count)
            #count_rsem_file.write( repr(global_pos) + '\t' + repr(base_count) + '\n')
            if base_count > max_base_count:
                max_base_count = base_count

    #pdb.set_trace()
    ExonPosSorted = sorted( ExonPosExp.keys() ) # , key=lambda x: x[0])
    for gp in ExonPosSorted:
        e_l = float(ExonPosExp.get(gp,0)[0])
        count_rsem_file.write( repr(gp) + '\t' + repr(e_l) + '\n')
        
    count_rsem_file.close()
    print('max L0: %f'%max_base_count)
    
    return count_rsem_address

    #--------------------------------------------------------------------------------------------------------
    
def Transcriptome_contribution_to_exon(local_pos, Tr, L):
    # for side effects
    # Tr = [transcript_expected_count, transcript_len]
    # local_pos = position w.r.t. transcrip_start
    transcript_expected_count = Tr[0]
    transcript_len = Tr[1]
    
    transcript_avrg_depth = transcript_expected_count * L / float(transcript_len)
    if transcript_len < L:
        transcript_peak = 0
    elif transcript_len == L:
        transcript_peak = transcript_avrg_depth
    else:
        transcript_peak = min(2, float(transcript_len)/(transcript_len - L) ) * transcript_avrg_depth 
    cov_depth = min(local_pos, transcript_len - local_pos, L) * transcript_peak/float(L) 
    
    return cov_depth    

 #--------------------------------------------------------------------------------------------------------
 
def main():
    #pdb.set_trace() #debug
    Path = Default_Ref_Path
    ref_address = Path + '/Chr15.fa'
    aligner ='tophat'
    L = 100
    
    
    readFQ_address_m = Path + '/Tar_m_read_l100.fastq'
    readFQ_address_p = Path + '/Tar_p_read_l100.fastq'
    readFQ_address = Path + '/Tar_read_l100.fastq'
    subprocess.call(cat +  ' '  + readFQ_address_m + ' ' + readFQ_address_p + ' > ' + readFQ_address, shell=True )
    
    
#    sam_address = Align(ref_address, readFQ_address_m, aligner, Path+'/tophat_out_m')
#    sam_address = Align(ref_address, readFQ_address_p, aligner, Path+'/tophat_out_p')
    sam_address = Align(ref_address, readFQ_address_m + ' ' + readFQ_address_p, aligner, Path+'/tophat_out')
    

    BED_address = Path + '/Chr15.bed'
    gtf_address = Path + '/hg19_chr15-UCSC.gtf'

    RSEM_result_address = RSEM(ref_address, readFQ_address, BED_address, gtf_address)
    exon_address = BED2Exon(BED_address)
    RSEM2Coverage(RSEM_result_address, exon_address,  L)
        
    
    
    
#    count_RSEM_address = short_Path +  '/count_rsem.txt'
#    cov_pdf_address = short_Path + '/cov_pdf.txt'
#    if os.path.isfile(count_RSEM_address) and flag:
#        print 'Abundance estimation is merged with the read count before!\n'
#    else:
#        exon_address = BED2Exon(BED_address)
#        output = RSEM2Coverage(RSEM_result_address, exon_address, count_address, L)
#        count_RSEM_address = output[0]
#        cov_pdf_address = output[1]
        #coverage_RSEM_address = RSEM2Coverage(RSEM_result_address, BED_address, L)
        #count_RSEM_address = InsEstExpInCOUNT(coverage_RSEM_address, count_address2, L)
#        print 'Abundance estimation is merging with the read count now!\n'
#        flag =0
    

'''
usage:

# rsem quantification
python ReadProcess.py --RSEM2 --ref fa_file --r1 r1_file --r2 r2_file --gtf gtf_file

# subsample a read file
python ReadProcess.py --SubsampleRead -r read_file -d dst -n n


'''
if __name__ == '__main__':

    args = sys.argv 

    if '--RSEM2' in args:
        RSEM2(args)
    elif '--SubsampleRead' in args:
        SubsampleRead(args)
    else:
        main()
    