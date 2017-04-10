#merge count_y{}.txt into count_y.txt
#debug why using para and not using para brings different snp res

from util import *

dir = '/data1/shunfu1/SNPCalling/data/data_GATK/2pass/split_sorted_sam_dedupped/count_y/'
dst = '/data1/shunfu1/SNPCalling/data/data_GATK/2pass/split_sorted_sam_dedupped/count_y/count_y.txt'

pdb.set_trace()

for i in range(20):
	src = dir + 'count_y%02d.txt'%i
	cmd = 'cat %s >> %s'%(src, dst)
	#pdb.set_trace()
	run_cmd(cmd)

pdb.set_trace()