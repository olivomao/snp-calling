'''
generate simulated reads using rsem

'''

from util import *
from Address import *
from batch_run_case1plus6 import enforce_unique_readID
import sys, pdb, random

def main():

    #get_exp_theta0('/data1/shunfu1/SNPCalling/data/rsem/exp.stat/exp.theta')
    
    #modi_tpm('/data1/shunfu1/SNPCalling/data/rsem/exp.isoforms.results', 1)

    '''
    args = '--gtf /data1/shunfu1/SNPCalling/data/gene_annotation/gencode.v19.annotation.gtf '+\
           '-r /data1/shunfu1/SNPCalling/data/ref_snp/chr15_[from_snp_m].fa '+\
           '-o /data1/shunfu1/SNPCalling/data/rsem_m/chr15_m'
    #prepare_rsem_ref(args.split())

    args = '--gen_sim_reads_rsem -R /data1/shunfu1/SNPCalling/data/rsem_m/chr15_m '+\
           '-E /data1/shunfu1/SNPCalling/data/rsem/exp '+\
           '-O /data1/shunfu1/SNPCalling/data/rsem_m_reads10m/ '+\
           '--modiTPM 0 '+\
           '-N 10000000'
    gen_sim_reads_rsem(args.split())
    '''

    '''
    args = '--gtf /data1/shunfu1/SNPCalling/data/gene_annotation/gencode.v19.annotation.gtf '+\
           '-r /data1/shunfu1/SNPCalling/data/ref_snp/chr15_[from_snp_p].fa '+\
           '-o /data1/shunfu1/SNPCalling/data/rsem_p/chr15_p'
    #prepare_rsem_ref(args.split())

    args = '--gen_sim_reads_rsem -R /data1/shunfu1/SNPCalling/data/rsem_p/chr15_p '+\
           '-E /data1/shunfu1/SNPCalling/data/rsem/exp '+\
           '-O /data1/shunfu1/SNPCalling/data/rsem_p_reads10m/ '+\
           '--modiTPM 1 '+\
           '-N 10000000'
    gen_sim_reads_rsem(args.split())
    '''

    #merge reads
    r1 = '/data1/shunfu1/SNPCalling/data/rsem_m_reads10m/sim_reads_1.fq'
    r2 = '/data1/shunfu1/SNPCalling/data/rsem_m_reads10m/sim_reads_2.fq'
    r_a  = '/data1/shunfu1/SNPCalling/data/rsem_m_reads10m/sim_reads_pooled.fq'
    args = '--merge_pe_reads --r1 %s --r2 %s -o %s'%(r1, r2, r_a)
    merge_pe_reads(args.split())

    r1 = '/data1/shunfu1/SNPCalling/data/rsem_p_reads10m/sim_reads_1.fq'
    r2 = '/data1/shunfu1/SNPCalling/data/rsem_p_reads10m/sim_reads_2.fq'
    r_b  = '/data1/shunfu1/SNPCalling/data/rsem_p_reads10m/sim_reads_pooled.fq'
    args = '--merge_pe_reads --r1 %s --r2 %s -o %s'%(r1, r2, r_b)
    merge_pe_reads(args.split())

    #pool merged reads
    r_c = '/data1/shunfu1/SNPCalling/data/sim_reads_10m.fq'
    args = '--cat_fq -1 %s -2 %s -o %s'%(r_a, r_b, r_c)
    cat_fq(args.split())

    return

'''
a wrapper of rsem-prepare-reference; note: mutated genome to use

usage:
--prepare_rsem_ref --gtf gtf_file
                   -r    ref_fa
                   -o    outPrefix
'''

def prepare_rsem_ref(args):

    #pdb.set_trace()

    gtf_file = args[args.index('--gtf')+1]
    ref_address = args[args.index('-r')+1]
    outPrefix = args[args.index('-o')+1]
    outFolder = '/'+'/'.join(outPrefix.split('/')[:-1])
    run_cmd('mkdir -p %s'%outFolder)

    rsem_index_command = rsemPath + '/rsem-prepare-reference'

    cmd = '%s '%rsem_index_command+\
          '--star -p 20 --star-path %s '%STAR_path +\
          '--gtf %s '%gtf_file +\
          '%s '%ref_address +\
          '%s '%outPrefix
    run_cmd(cmd)

    return


'''
usage:
--gen_sim_reads_rsem -R inputRsemRef
                     -E inputRsemExp
                     -O outputFolder
                     --modiTPM 0/1/etc (0 - OFF, 1 - rand, etc)
                     -N numReads

note:
- inputRsemRef: contains ref generation prefix using rsem
- inputRsemExp: contains exp estimation prefix using real reads
- outputFolder: store generated sim reads
--modi TPM
  0 - no modification
  1 - each tr gets rand TPM within original min~max range
'''
def gen_sim_reads_rsem(args):

    #prev from refShannon
    #rsem-simulate-reads \
    #/data1/shunfu1/ref_shannon_modi/data/sgRefShannon/rsem/align_star/snyder/ref \
    #exp.stat/exp.model \
    #exp.isoforms.results \
    #theta0(first value of the third line of the file 'exp.stat/exp.theta' ==> 0.0429317239932184) \
    #N (the 4th number of the first line of the file 'sample_name.stat/sample_name.cnt' ==> 150M) \
    #output_name(path/to/sim_reads ==> /data1/shunfu1/ref_shannon_modi/data/sgRefShannon/rsem/sim_star/snyder/sim_reads) \
    #--seed 0

    #pdb.set_trace()

    inputRsemRef = args[args.index('-R')+1]

    inputRsemExp = args[args.index('-E')+1]
    sample_name = [i for i in inputRsemExp.split('/') if i != ''][-1]
    exp_model = '%s.stat/%s.model'%(inputRsemExp, sample_name)
    exp_res = '%s.isoforms.results'%(inputRsemExp)
    exp_theta = '%s.stat/%s.theta'%(inputRsemExp, sample_name)
    exp_theta0 = get_exp_theta0(exp_theta) ##float

    outputFolder = args[args.index('-O')+1]
    run_cmd('mkdir -p %s'%outputFolder)
    outputReadsPre = outputFolder + '/sim_reads'

    modiTPM = int(args[args.index('--modiTPM')+1])
    exp_res = modi_tpm(exp_res, modiTPM) ##)

    N = int(args[args.index('-N')+1])


    cmd = 'rsem-simulate-reads '+\
          '%s '%inputRsemRef+\
          '%s '%exp_model+\
          '%s '%exp_res+\
          '%f '%exp_theta0+\
          '%d '%N+\
          '%s '%outputReadsPre+\
          '--seed 0'
    #pdb.set_trace()
    run_cmd(cmd)

    return

#first value of the third line of the file 'exp.stat/exp.theta'
def get_exp_theta0(exp_theta_file):

    with open(exp_theta_file, 'r') as f:
        #pdb.set_trace()

        line = f.readline()
        line = f.readline()
        line = f.readline()
        exp_theta0 = float(line.split()[0])
        
        #pdb.set_trace()

    return exp_theta0

#if modiTPM==1, modi TPM of exp_res file, store into a new file and return its address
def modi_tpm(exp_res, modiTPM):

    if modiTPM==0:
        return exp_res
    elif modiTPM>1:
        print('modi_tmp exception undefined modiTPM: %d'%modiTPM)
        pdb.set_trace()

    exp_res_new = exp_res + '.modiTPM'
    tpm_vals = []
    #pdb.set_trace()
    with open(exp_res, 'r') as f:
        f.readline()
        for line in f:
            tpm = float(line.split()[5])
            tpm_vals.append(tpm)
        print('max tpm: %f'%max(tpm_vals))

    #pdb.set_trace()
    with open(exp_res, 'r') as f1, open(exp_res_new, 'w') as f2:

        line = f1.readline()
        f2.write(line)

        for line in f1:
            tokens = line.split()
            rand_idx = random.randint(0, len(tpm_vals)-1) #a and b inclusive
            tokens[5] = '%.2f'%tpm_vals[rand_idx]
            tokens[3] = '.'; tokens[4] = '.'; tokens[6] = '.'; tokens[7] = '.'
            del tpm_vals[rand_idx]
            f2.write('\t'.join(tokens)+'\n')
    return exp_res_new

'''
merge reads: PE --> SE

currently we treat PE reads as SE reads
esp convert readA ('/1'-->'a') and readB ('/2'-->'b')

usage:
--merge_pe_reads --r1 r1 --r2 r2 -o path/to/outFile
'''
def merge_pe_reads(args):

    r1 = args[args.index('--r1')+1]
    r2 = args[args.index('--r2')+1]
    outFile = args[args.index('-o')+1]

    with open(r1,'r') as f1, open(r2,'r') as f2, open(outFile, 'w') as f3:

        for line in f1:
            if line[0]=='@':
                #pdb.set_trace()
                line = line.strip()[:-2]+'a'
                f3.write(line+'\n')
            else:
                f3.write(line)
        print('%s added'%r1)

        for line in f2:
            if line[0]=='@':
                #pdb.set_trace()
                line = line.strip()[:-2]+'b'
                f3.write(line+'\n')
            else:
                f3.write(line)
        print('%s added'%r2)

    #ensure uniq id
    enforce_unique_readID(outFile, flag=True)

    return
'''
cat several fq reads files (e.g. each from merge_pe_reads) into one with uniq read ids

usage:
--cat_fq -1 r1 -2 r2 ... -o r
'''
def cat_fq(args):

    cmd = 'cat '
    idx = 1
    while '-%d'%idx in args:
        cmd += args[args.index('-%d'%idx)+1] + ' '
        idx += 1

    outFile = args[args.index('-o')+1]
    cmd += '> '+outFile
    run_cmd(cmd)

    #ensure uniq id
    enforce_unique_readID(outFile, flag=True)
    print('')

    return

if __name__ == "__main__":

    args = sys.argv

    #pdb.set_trace()

    if '--merge_pe_reads' in args:
        merge_pe_reads(args)
    elif '--cat_fq' in args:
        cat_fq(args)
    else:
        main()