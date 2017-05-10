import pdb
from old_code.util import run_cmd

def run_varscan():

    varscan_dir = '/home/shunfu1/software/VarScan/'
    #ref_file = '/data1/shunfu1/SNPCalling//data_large_0_idealCov/Chr15.fa' 
    ref_file = '/data1/shunfu1/SNPCalling/SimSNPs_MultiRun_5/reads_N100K_L100_Err0.00/rsem/Chr15.transcripts.fa'
    #bam_file = '/data1/shunfu1/SNPCalling/SimSNPs_MultiRun_5/reads_N100K_L100_Err0.00/gatk/2pass/dedupped.bam'
    #bam_file = '/data1/shunfu1/SNPCalling/SimSNPs_MultiRun_5/reads_N100K_L100_Err0.00/rsem/Chr15.genome.sorted.bam'
    bam_file = '/data1/shunfu1/SNPCalling/SimSNPs_MultiRun_5/reads_N100K_L100_Err0.00/rsem/Chr15.transcript.sorted.bam'
    res_file_name = 'run_5_n100k_err000_rsem_transcriptome_res.txt'

    cmd = 'samtools mpileup -f %s %s | java -jar %s/VarScan.v2.3.9.jar pileup2snp > %s/%s'%\
          (ref_file, bam_file, varscan_dir, varscan_dir, res_file_name)
    pdb.set_trace()
    run_cmd(cmd)

    return

def analyze():

    #cmp_res_file = '/home/shunfu1/software/VarScan/run_5_n100k_err000_dedupped_res.txt'
    cmp_res_file = '/home/shunfu1/software/VarScan/run_5_n100k_err000_rsem_res.txt'

    count = '/data1/shunfu1/SNPCalling/SimSNPs_MultiRun_5/SNP_p.txt'
    pos_idx = 0

    cmp_res = set()
    with open(cmp_res_file, 'r') as f:
        f.readline()
        for line in f:
            if line[0]=='#': continue
            tokens = line.split()
            cmp_res.add(int(tokens[1])-1)

    pos_count = set()
    with open(count, 'r') as f:
        for line in f:
            if line[0]=='#': continue
            tokens = line.split()
            pos_count.add(int(tokens[pos_idx]))

    #pdb.set_trace()

    res = cmp_res.intersection(pos_count)
    print(len(res))

    pdb.set_trace()

    return

if __name__ == "__main__":

    run_varscan()

    #analyze()


