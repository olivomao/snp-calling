#test filter_snp_lam_half_filt3
from snp_res_statistics import filter_snp_lam_half_filt3
from util import run_cmd, run_cmds

#snp_res_address = '/data1/shunfu1/SNPCalling/data/caller_output_snp_found_Chr15.genome.sorted_n_T1_para_filt.txt'
#snp_res_address = '/data1/shunfu1/SNPCalling/data/caller_output_snp_found_Chr15.genome.sorted_n_T1_para.txt'
snp_res_address = '/data1/shunfu1/SNPCalling/data/caller_output_snp_found_dedupped_T1_para.txt'
#count_altInfo_address = '/data1/shunfu1/SNPCalling/data/rsem/split_sorted_sam_Chr15.genome.sorted_n/count_y_altInfo/count_Chr15.genome.sorted_n_altInfo.txt'
count_altInfo_address = '/data1/shunfu1/SNPCalling/data/data_GATK/2pass/split_sorted_sam_dedupped/count_y_altInfo/count_dedupped_altInfo.txt'
#count_abs_address = '/data1/shunfu1/SNPCalling/data/rsem/split_sorted_sam_Chr15.genome.sorted_n/count_y/count_Chr15.genome.sorted_n.txt'
count_abs_address = '/data1/shunfu1/SNPCalling/data/data_GATK/2pass/split_sorted_sam_dedupped/count_y/count_dedupped.txt'
filtSameAb = True
rocThre = [0.0, 0.1, 0.5, 1.0, 5.0, 10.0]
#rocThre = [0.0]
nJobs = 20

filter_snp_lam_half_filt3(snp_res_address,
                          count_altInfo_address,
                          count_abs_address,
                          filtSameAb=filtSameAb,
                          rocThre=rocThre)



#for group uniq/non uniq snp purpose

c1f = '/data1/shunfu1/SNPCalling/data/data_GATK/2pass//split_sorted_sam_dedupped/count_y/count_dedupped.txt'
#c1f = '/data1/shunfu1/SNPCalling/data/data_GATK/2pass//count_dedupped.txt'

c2f = '/data1/shunfu1/SNPCalling/data/data_GATK/2pass//split_sorted_sam_dedupped/count_y_altInfo/count_dedupped_altInfo.txt'
#c2f = '/data1/shunfu1/SNPCalling/data/data_GATK/2pass//count_dedupped_altInfo.txt'

cmds = []
for rocT in rocThre:

    filt_snp_file = snp_res_address[:-4]+'_rocT_%.2f.txt'%rocT

    cmd = 'python snp_analysis.py --loadSnpInfo  -L1 O -F1 %s '%filt_snp_file+\
          ' -L2 m -F2 /data1/shunfu1/SNPCalling/data/SNP_m.txt '+\
          ' -L3 p -F3 /data1/shunfu1/SNPCalling/data/SNP_p.txt '+\
          ' -L4 G -F4 /data1/shunfu1/SNPCalling/data/data_GATK/GATK_out/raw_variants.vcf.txt '+\
          ' -C1 %s '%c1f+\
          ' -C2 %s '%c2f+\
          ' --snpLog /data1/shunfu1/SNPCalling/data/snpLog_[dedupped][filt3_filtSameAb%s_rocT%.2f]_T1_c1.txt '%(filtSameAb, rocT)+\
          ' --snpSum /data1/shunfu1/SNPCalling/data/snpSum_[dedupped][filt3_filtSameAb%s_rocT%.2f]_T1_c1.txt '%(filtSameAb, rocT)+\
          ' --snpSum2 /data1/shunfu1/SNPCalling/data/snpSum2_[dedupped][filt3_filtSameAb%s_rocT%.2f]_T1_c1.txt '%(filtSameAb, rocT)
    cmds.append(cmd)
run_cmds(cmds,nJobs)