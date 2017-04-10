import sys, pdb

from tool_address import FA2FQ_address
from old_code.batch_run_case1plus6 import enforce_unique_readID
from old_code.util import run_cmd
from old_code_multiShannon.sim_reads import sim_reads_g_b

def main():

    pdb.set_trace()
    
    return

'''
example:

python sim_data_generator.py --read_generation -g /data1/shunfu1/SNPCalling/snp20_reads100k_10/Tar_m.txt 
                             -b /data1/shunfu1/SNPCalling/snp20_reads100k_10/hg19_chr15-UCSC-sorted.bed
                             -O /data1/shunfu1/SNPCalling/snp20_reads100k_10_abSNPcode/reads_N100k_L100_Err0 
                             -n 100000 -l 100 -r 0 
                             --suffix _m --toFq
'''
def read_generation(args):

    sim_reads_g_b(args)

    #add suffix
    out_dir = args[args.index('-O')+1]

    suf = args[args.index('--suffix')+1]      
    
    src_dst_pairs = []

    reads = '%s/reads.fa'%(out_dir); readsNew = '%s/reads%s.fa'%(out_dir, suf)
    src_dst_pairs.append([reads, readsNew])
    
    if '--exp_path' not in args:
        exp = '%s/intermediate/exp.txt'%(out_dir); expNew = '%s/intermediate/exp%s.txt'%(out_dir, suf)
        src_dst_pairs.append([exp, expNew])

    reads_bed = '%s/intermediate/reads.bed'%(out_dir); reads_bedNew = '%s/intermediate/reads%s.bed'%(out_dir, suf)
    src_dst_pairs.append([reads_bed, reads_bedNew])

    for i in range(len(src_dst_pairs)):
        run_cmd('mv %s %s'%(src_dst_pairs[i][0], src_dst_pairs[i][1]))

    #fa to fq
    readsFq = '%s/reads%s.fq'%(out_dir, suf)
    run_cmd( 'perl ' + FA2FQ_address + ' ' + readsNew + ' > '  + readsFq )

    return

def merge_reads(args):

    reads_m = args[args.index('-m')+1]
    reads_p = args[args.index('-p')+1]
    merged_reads = args[args.index('-o')+1]

    run_cmd('cat '  + reads_m + ' ' + reads_p + ' > ' + merged_reads)

    if '--uniqID' in args:
        enforce_unique_readID(merged_reads, flag=True)
        
    return

'''
usage:

# use RNASeqReadSimulator to generate reads per target allele
#
# we re-use the simulator written in multi shannon
# for snp calling purpose, we need to add suffix (e.g. _m, _p) and convert fa to fq
#
# note:
#
# .fa, \.sorted.bed
#                               --> out_dir/intermediate/exp[_suffix].txt (if --exp_path unspecified), reads[_suffix].bed 
#                               --> out_dir/reads[_suffix].[fa or fq] (SE)
# 
# .sorted.bed needs 1st line to be "# trNum (trNum) avgTrLen (avgTrLen)" for '-c readCoverage'
#

python sim_data_generator.py --read_generation -g genomeFile -b trBedFile -O out_dir
                             (-n numReads or -c readCoverage) -l readLength -r errRate
                             --suffix suf --toFq
                             [--exp_path exp_path]

# pool reads sampled from maternal and paternal alleles together
#

python sim_data_generator.py --merge_reads -m reads_m.fq -p reads_p.fq -o merged_reads.fq [--uniqID]

'''
if __name__ == "__main__":

    args = sys.argv

    if '--read_generation' in args:
        read_generation(args)
    elif '--merge_reads' in args:
        merge_reads(args)
    else:
        main()