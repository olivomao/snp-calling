"""
run variant calling algo (note4) using .sam from STAR aligner
in order to compare our caller with GATK's Haplotype Caller
(e.g. miss-detection, false positive)
"""

from count_read_lambda import *
from final_caller19thjuly_m import *
import pdb
from Address import *

ref_address = "/home/sreeramkannan/singleCell/SNP_Calling_Summer15/data_0814/Chr15.fa"

if 1:
    #sam_address = "/home/sreeramkannan/singleCell/SNP_Calling_Summer15/data_0814/data_GATK/2pass/Aligned.out.sam"
    #sam_address = "/home/sreeramkannan/singleCell/SNP_Calling_Summer15/data_0814/tophat_out/accepted_hits.sam"
    #sam_address = "/home/sreeramkannan/singleCell/SNP_Calling_Summer15/data_0814/data_GATK/1pass/Aligned.out.sam"
    #sam_name = "Aligned.out.sam"
    #sam_name = "accepted_hits.sam"
    sam_address = "/home/sreeramkannan/singleCell/SNP_Calling_Summer15/data_0814/data_GATK/2pass/Aligned.out.sam"
else:
    flder = "/home/sreeramkannan/singleCell/SNP_Calling_Summer15/data_0814/data_GATK/2pass/"
    name = "rg_added_sorted"
    #name = "dedupped"
    #name = "split"
    bam_address = flder+name+".bam"
    sam_address = flder+name+".sam"
    SamtoolsProgram = SamPath + '/samtools'
    subprocess.call(SamtoolsProgram +  ' view -h ' + bam_address + ' > ' + sam_address, shell=True )
    #subprocess.call(SamtoolsProgram +  ' sort -O sam -o ' + sam_address + ' -T temp ' + bam_address, shell=True )
    pdb.set_trace()

coverage_address = "/home/sreeramkannan/singleCell/SNP_Calling_Summer15/data_0814/coverage.txt"

#count_fn = '/count_cross_check_accepted_hits_modCase7.txt'
#count_fn = '/count_cross_check_star1pass.txt'
#count_fn = '/count_cross_check_star1pass_modCase7.txt'
#count_fn = '/count_cross_check_rg_added_sorted.txt'
#count_fn = '/count_cross_check_dedupped.txt'
#count_fn = '/count_cross_check_split_trimS3.txt'
#count_fn = '/count_cross_check_split.txt'
#count_fn = '/count_altMap_cross_check_star1pass.txt'
#count_fn = '/count_altMap_cross_check_split_0915.txt'
#count_fn = '/count_altMap_cross_check_rg_added_sorted.txt'
count_fn = '/count_altMap_cross_check_star2pass.txt'

pdb.set_trace()
generate_count_file(sam_address, ref_address, coverage_address, count_fn) #en_debug_0817
#count_rel_address = "/tophat_out/" + count_fn
#count_rel_address = "/data_GATK/1pass/" + count_fn
count_rel_address = "/data_GATK/2pass/" + count_fn

#caller_output_res_fn = '/caller_output_altMap_cross_check_star1pass.txt'
#caller_output_exception_fn = '/caller_output_exception_altMap_cross_check_star1pass.txt'
#caller_output_snp_found_fn = '/caller_output_snp_found_altMap_cross_check_star1pass.txt'
#caller_output_res_fn = '/caller_output_altMap_cross_check_split_0915.txt'
#caller_output_exception_fn = '/caller_output_exception_altMap_cross_check_split_0915.txt'
#caller_output_snp_found_fn = '/caller_output_snp_found_altMap_cross_check_split_0915.txt'
#caller_output_res_fn = '/caller_output_altMap_cross_check_dedupped.txt'
#caller_output_exception_fn = '/caller_output_exception_altMap_cross_check_dedupped.txt'
#caller_output_snp_found_fn = '/caller_output_snp_found_altMap_cross_check_dedupped.txt'
#caller_output_res_fn = '/caller_output_altMap_cross_check_rg_added_sorted.txt'
#caller_output_exception_fn = '/caller_output_exception_altMap_cross_check_rg_added_sorted.txt'
#caller_output_snp_found_fn = '/caller_output_snp_found_altMap_cross_check_rg_added_sorted.txt'
caller_output_res_fn = '/caller_output_altMap_cross_check_star2pass.txt'
caller_output_exception_fn = '/caller_output_exception_altMap_cross_check_star2pass.txt'
caller_output_snp_found_fn = '/caller_output_snp_found_altMap_cross_check_star2pass.txt'
#pdb.set_trace()

## ---------- caller
pdb.set_trace()
final_caller(count_rel_address, 
             caller_output_res_fn, 
             caller_output_exception_fn,
             caller_output_snp_found_fn)
