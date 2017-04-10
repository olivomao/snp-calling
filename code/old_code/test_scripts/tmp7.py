#debug rsem sam
#test count_read_lambda

from count_read_lambda import *

sam = 'tmp/rsemSam_debug/sample_rsemSam.sam'
ref = 'tmp/data/Chr15.fa'
cov_addr = 'tmp/data/count_rsem.txt'

generate_count_file(sam, ref, cov_addr)
