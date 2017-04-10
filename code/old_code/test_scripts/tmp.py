'''
copy snpLog, snpSum, snpSum2 from data to code/tmp/ to local pc for analysis
'''
from util import *
import os

src_rootFld = '/data1/shunfu1/SNPCalling/'
dst_rootFld = '/home/shunfu1/SNP_Calling_Summer15/snp-calling/code/tmp/'
prefix = 'snp20_reads100k_'

pdb.set_trace()

for i in range(12):

    src_fld = '%s%s%d/'%(src_rootFld, prefix, i)
    dst_fld = '%s%s%d/'%(dst_rootFld, prefix, i)

    f1 = src_fld + 'snpLog.txt'
    f2 = src_fld + 'snpSum.txt'
    f3 = src_fld + 'snpSum2.txt'

    if os.path.exists(f1)==True and os.path.exists(f2)==True and os.path.exists(f3)==True:
        cmd = 'mkdir -p %s/'%(dst_fld)
        run_cmd(cmd)

        cmd = 'cp %s %s/snpLog.txt'%(f1, dst_fld)
        run_cmd(cmd)

        cmd = 'cp %s %s/snpSum.txt'%(f2, dst_fld)
        run_cmd(cmd)

        cmd = 'cp %s %s/snpSum2.txt'%(f3, dst_fld)
        run_cmd(cmd)



