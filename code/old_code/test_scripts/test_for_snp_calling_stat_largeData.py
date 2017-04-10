'''
call case6 (ours, in parallel, multi thresholds) and case2 (GATK)
'''

import sys
from util import run_cmd
import pdb
from batch_run_case1plus6 import merge_coverage
import os

#condition = 'snp20_reads100k_11_varyT'
#condition = 'snp20_reads100k_4_T1_useRsemSam_filtZW0.1'
#condition = 'snp20_reads100k_7_T1_useRsemSam_filtZW0.1'
#condition = 'snp20_reads100k_10_useRsemSam_filtZW0.1andSameZW'
condition = ''
#condition = 'snp20_reads100k_10_varyT'
#condition = 'snp1k_reads10m_varyT_rsemCount2'
#condition = 'snp1k_reads10m_varyT_idealCov'
#condition = 'snp1k_reads10m_T1_useRsemSam_filtZW0.1'

#condition = 'snp1k_reads10m_sameExp_varyT_rsemCount'
#condition = 'snp1k_reads10m_sameExp_varyT_rsemCount'
#condition = 'snp1k_reads10m_0_rsemCount_useRsemSam_filtZW0.1andSameZW'

#Thre = [1,5]
#Thre = [1,3,5,10,25]
Thre = [1]#, 100, 200, 300]
iter = -1

ref_address = '/data1/shunfu1/SNPCalling/data/Chr15.fa'
covAddress = '/data1/shunfu1/SNPCalling/data/count_rsem.txt'
#covAddress = '/data1/shunfu1/SNPCalling/data/coverage.txt' #covAddress = merge_coverage('/data1/shunfu1/SNPCalling/data/coverage_m.txt', '/data1/shunfu1/SNPCalling/data/coverage_p.txt', flag=True)
sam_address = '/data1/shunfu1/SNPCalling/data/data_GATK/2pass/dedupped.sam'
#sam_address = '/data1/shunfu1/SNPCalling/data/rsem/Chr15.genome.sorted_n.sam'
#pdb.set_trace()
if 'rsem' in sam_address and os.path.exists(sam_address)==False:
    sam_fn = sam_address.split('/')[-1][:-4]
    sam_dir = sam_address[0:len(sam_address)-len(sam_fn)-4]
    
    cmd = 'rsem-tbam2gbam %s/Chr15 %s/Chr15.transcript.bam %s/Chr15.genome.bam -p 20'%(sam_dir, sam_dir, sam_dir)
    run_cmd(cmd)

    cmd = 'samtools sort -n -o %s/Chr15.genome.sorted_n.sam -@ 20 %s/Chr15.genome.bam'%(sam_dir, sam_dir)
    run_cmd(cmd)

output_address = '/data1/shunfu1/SNPCalling/data/' #not used yet

#generate data
#cmd = 'python batch_run_case1plus6.py'
#run_cmd(cmd)

#snp calling
for thre in Thre:

    iter += 1
    if iter==0:
        dupRun = 0
    else:
        dupRun = 1


    testArg = '-r %s -c %s -s %s -O %s -T %d -p 20 --dupRun %d'%(ref_address, covAddress, sam_address, output_address, thre, dupRun)
    cmd = 'python batch_run_parallel_modi.py %s'%testArg
    run_cmd(cmd)

    #pdb.set_trace()
    cmd = 'python test_snp_analysis.py --condition %s '%condition+\
          ' --copyRes '+\
          ' --thre %d '%thre+\
          ' --para --compare 0'
    #run_cmd(cmd)

#GATK & comparison
#pdb.set_trace()
cmd = 'python batch_run_case2.py'
#run_cmd(cmd)

#pdb.set_trace()
cmd = 'python test_snp_analysis.py --condition %s '%condition+\
      ' --copyRes '+\
      ' --thre 1 '+\
      ' --para --compare 1 '+\
      ' --sam %s'%sam_address
#run_cmd(cmd)
