from util import run_cmd

c1f = '/data1/shunfu1/SNPCalling/data/data_GATK/2pass//split_sorted_sam_dedupped/count_y/count_dedupped.txt'
#c1f = '/data1/shunfu1/SNPCalling/data/data_GATK/2pass//count_dedupped.txt'

c2f = '/data1/shunfu1/SNPCalling/data/data_GATK/2pass//split_sorted_sam_dedupped/count_y_altInfo/count_dedupped_altInfo.txt'
#c2f = '/data1/shunfu1/SNPCalling/data/data_GATK/2pass//count_dedupped_altInfo.txt'

cmd = 'python snp_analysis.py --loadSnpInfo  -L1 O -F1 /data1/shunfu1/SNPCalling/data/caller_output_snp_found_Chr15.genome.sorted_n_T1_para_filt.txt '+\
      ' -L2 m -F2 /data1/shunfu1/SNPCalling/data/SNP_m.txt '+\
      ' -L3 p -F3 /data1/shunfu1/SNPCalling/data/SNP_p.txt '+\
      ' -L4 G -F4 /data1/shunfu1/SNPCalling/data/data_GATK/GATK_out/raw_variants.vcf.txt '+\
      ' -C1 %s '%c1f+\
      ' -C2 %s '%c2f+\
      ' --snpLog /data1/shunfu1/SNPCalling/data/snpLog_[Chr15.genome.sorted_n][filtSameZW_groupByDedupCountAlt]_T1_c1.txt '+\
      ' --snpSum /data1/shunfu1/SNPCalling/data/snpSum_[Chr15.genome.sorted_n][filtSameZW_groupByDedupCountAlt]_T1_c1.txt '+\
      ' --snpSum2 /data1/shunfu1/SNPCalling/data/snpSum2_[Chr15.genome.sorted_n][filtSameZW_groupByDedupCountAlt]_T1_c1.txt'
run_cmd(cmd)