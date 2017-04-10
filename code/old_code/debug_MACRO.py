import csv
from progress.bar import Bar

#---------- MACROs
en_debug = 1
#en_debug_0811 = 1
#
# Synthesis -> ReadGeneration: modify the names of read files (.bed and .fasta)
# use merged m & p reads to be used by tophat for read alignment
# final caller exception
#
#en_debug_0812 = 1
# 
# bed & exp consistency
# coverage.txt issue (debug and modification)
#
#en_debug_0814
#
# sim_stat.py
# modify sim_res_statistics.py
# filter coverage
# ReadGeneration use BED (not BED sorted) files
# generate data files for different qt & N values to check how SNPs are covered by reads (test_SNP_covered_by_reads, find_pos_0)

en_check_multi_repeats = 0 #to understand why reads are grouped into nonseg/seg and noninverse/inverse

#en_debug_0817


#---------- functions

def FilterCoverage(coverage_address, coverage_address_filtered, line_cover_threshold):
    
    coverage_file_filtered = open(coverage_address_filtered, 'w+')
    
    idx0 = 0
    acc_len0 = 0
    
    num_lines = sum(1 for line in open(coverage_address))
    bar = Bar('filter coverage', max=num_lines) 
    
    with open(coverage_address) as coverage_file:
        
            reader = csv.reader(coverage_file, delimiter='\t')
            for row in reader:
                bar.next()
                
                e_stt = int(row[1])
                e_stp = int(row[2])
                e_len = e_stp - e_stt
                e_cov = float(row[4])
                
                if e_cov > line_cover_threshold:
    
                    idx0 = idx0 + 1                
                    acc_len0 = acc_len0 + e_len
                    
                    row_new = [repr(idx0), repr(e_stt), repr(e_stp), repr(acc_len0), repr(e_cov)]
                    coverage_file_filtered.write('\t'.join(row_new)+'\n')
                
                """
                check if there's no threshold
                idx0 = idx0 + 1                
                acc_len0 = acc_len0 + e_len
                
                coverage_file_filtered.write('\t'.join(row)+'\n')
                coverage_file_filtered.write('\trecalculated idx=%d, acc_len=%d\n\n'%(idx0, acc_len0))
                """
                
            bar.finish()
    
    coverage_file_filtered.close()
    
    return