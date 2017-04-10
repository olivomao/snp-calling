'''
used to analyze sam file
'''

import sys, pdb
import numpy as np
import re

'''
--stat
-s sam_file: need to be sorted according to name
-q mapq_dmp_file : to check histogram of mapq
-z zw_dmp_file : RSEM pads a new tag ZW:f:value, where value is a single precision floating number representing the posterior probability of that alignment being the true mapping of a read
--remove0zw : for RSEM, throw away alignment with 0 zw or no zw values

Get stat from sam file
-- # of uniq aligned reads
-- # of mult aligned reads
-- dump mapq fields (-q) to external mapq_dmp_file for further histogram plot (by hist_view)
-- for RSEM: dump ZW:f:val (-z) to external zw_dmp_file for further histogram plot (by hist_view)
-- for RSEM: skip count alignments w/ ZW:f:0 or no ZW:f:val field using --remove0zw option
'''
def stat(args):

    sam_file = args[args.index('-s')+1]

    nLines=sum([1 for l in open(sam_file,'r')]); T=nLines/100; p=0; q=0;
    print('%d lines in %s'%(nLines, sam_file))

    #MAPQ
    if '-q' in args:
        mapq_dmp_file = args[args.index('-q')+1]
        mapq_dmp_file = open(mapq_dmp_file, 'w')
    else:
        mapq_dmp_file = None

    #ZW val
    if '-z' in args:
        zw_dmp_file = args[args.index('-z')+1]
        zw_dmp_file = open(zw_dmp_file, 'w')
        cnt_no_zw = 0
        cnt_zw = 0
    else:
        zw_dmp_file = None
        cnt_no_zw = 0
        cnt_zw = 0

    #--remove0zw
    if '--remove0zw' in args:
        remove0zw = True
        cnt_remove0zw = 0
    else:
        remove0zw = False
        cnt_remove0zw = 0

    num_uniq_aligned_reads = 0
    num_mult_aligned_reads = 0

    with open(sam_file, 'r') as sF:

        read_group = []

        for line in sF:
            p += 1
            if p>=T: p=0; q+=1; sys.stdout.write('\r'); sys.stdout.write('%d %% processed (analyze_sam.stat)'%q); sys.stdout.flush()

            if line[0]=='@': continue # skip head info

            #pdb.set_trace()

            #MAPQ
            if mapq_dmp_file is not None:
                mapq = int(line.split()[4])
                mapq_dmp_file.write('%d\n'%mapq)

            #ZW val
            if zw_dmp_file is not None or remove0zw==True: #check zw val
                b=re.search('ZW:f:([\S]+)', line)
                if b is not None:
                    v = b.group(1)
                    #print(v)
                    #if float(v)==0.0: pdb.set_trace()
                    if zw_dmp_file is not None: zw_dmp_file.write('%e\n'%float(v))
                    cnt_zw += 1
                else:
                    v = -1
                    cnt_no_zw += 1
                    #pdb.set_trace()

            if remove0zw==True and float(v)<=0.0:
                #print('%d %s'%(cnt_remove0zw, line))
                cnt_remove0zw += 1
                continue

            #single/mult aligned reads
            if read_group == []:
                read_group.append(line)
            else:
                prev_read_id = read_group[-1].split()[0]
                curr_read_id = line.split()[0]
                if prev_read_id == curr_read_id:
                    read_group.append(line)
                else:
                    if len(read_group)==1:
                        num_uniq_aligned_reads += 1
                    else:
                        num_mult_aligned_reads += 1
                    read_group = [line]

        if len(read_group)==1:
            num_uniq_aligned_reads += 1
        else:
            num_mult_aligned_reads += 1

    num_tot_aligned_reads = num_uniq_aligned_reads + num_mult_aligned_reads
    print('\nnum_uniq_aligned_reads: %d out of %d (%f)'%(num_uniq_aligned_reads, num_tot_aligned_reads, float(num_uniq_aligned_reads)/ num_tot_aligned_reads))
    print('num_mult_aligned_reads: %d out of %d (%f)'%(num_mult_aligned_reads, num_tot_aligned_reads, float(num_mult_aligned_reads)/ num_tot_aligned_reads))

    if mapq_dmp_file is not None:
        mapq_dmp_file.close()

    if zw_dmp_file is not None:
        zw_dmp_file.close()
        print('%d zw (%f) and %d no zw (%f)'%(cnt_zw, float(cnt_zw)/(cnt_zw+cnt_no_zw), cnt_no_zw, float(cnt_no_zw)/(cnt_zw+cnt_no_zw)))

    if remove0zw:
        print('%d lines contain remove0zw'%(cnt_remove0zw))

    return

def cov(st, type='int'):

    if type=='int':
        return int(st)
    elif type=='float':
        return float(st)

'''
--hist_view
-f file
-t type

plot histogram from file (each val per line, type specified by -t) at local pc
'''
def hist_view(args):
    import matplotlib.pyplot as plt

    file = args[args.index('-f')+1]
    type = args[args.index('-t')+1]

    fig, ax = plt.subplots()

    mapq_list = []

    with open(file, 'r') as f:
        for line in f:
            mapq_list.append(cov(line.strip(), type))

    ax.hist([mapq_list])
    ax.set_yscale('log')
    plt.show()

    #pdb.set_trace()

def main():

    args = sys.argv

    '''
    args = ['--stat',
            '-s', 'tmp/analyze_sam/test.sam',
            '-q', 'tmp/analyze_sam/mapq_dmp.txt',
            '-z', 'tmp/analyze_sam/zw_dmp.txt']
    '''

    '''
    args = ['--hist_view',
            '-f', 'tmp/analyze_sam/rsem_sam_large_mapq.txt',#'tmp/analyze_sam/dedupped_large_mapq.txt',#
            '-t', 'int']
    '''

    '''
    args = ['--hist_view',
            '-f', 'tmp/analyze_sam/rsem_sam_large_zw.txt', #'tmp/analyze_sam/rsem_sam_mapq.txt',
            '-t', 'float']
    '''

    if '--stat' in args:
        stat(args)
    elif '--hist_view' in args:
        hist_view(args)

    #pdb.set_trace()

    return

if __name__ == "__main__":
    main()