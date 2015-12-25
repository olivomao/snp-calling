import pdb
import subprocess
from final_caller19thjuly_m import *

"""
Created on Fri Oct  2 17:37:58 2015

@author: olivo
"""

print('split count file for parallel processing')
pdb.set_trace()

Default_Ref_Path = '/home/sreeramkannan/singleCell/SNP_Calling_Summer15/data_SNP1k_Reads10M/'
tophat_dir = '/tophat_out/'
count_split_dir = Default_Ref_Path+tophat_dir+'/count_y03/'
subprocess.call('mkdir -p '+count_split_dir, shell = True)

count_abs_address = Default_Ref_Path+tophat_dir+'/count_y/count_y03.txt'

num_p = 20

do_split = False
if do_split == True:
    num_lines = sum(1 for line in open(count_abs_address))

    num_lines_per_file = int(num_lines / num_p) + 1

    subprocess.call('split -l '+repr(num_lines_per_file)+' -d '+count_abs_address, shell = True) #split -l 167880 -d $count_fn
    subprocess.call('mv x* '+count_split_dir, shell = True) #mv x* $count_split_dir

cmd = 'parallel python final_caller19thjuly_m.py '+\
      tophat_dir+'/count_y03/ '+\
      'x{} '+\
      'caller_op_y03_{}.txt '+\
      'caller_op_exception_y03_{}.txt '+\
      'caller_op_snp_y03_{}.txt '+\
      ':::'
for i in range(num_p):
    cmd = cmd + ' %02d'%i
    
subprocess.call(cmd, shell = True)

merge_caller_res(Default_Ref_Path,
                 'caller_op_y03_', 'caller_op_exception_y03_', 'caller_op_snp_y03_', num_p,
                 'caller_op_y03_para.txt', 'caller_op_exception_y03_para.txt', 'caller_op_snp_y03_para.txt')
