import pdb
from old_code.util import run_cmd

def main():

    ref_file = '/data1/shunfu1/SNPCalling/data_large_0_idealCov/Chr15.fa'
    bed_file = '/data1/shunfu1/SNPCalling/data_large_0_idealCov/hg19_chr15-UCSC.bed'
    bam_file = '/data1/shunfu1/SNPCalling/SimSNPs_MultiRun_5/reads_N100K_L100_Err0.00/rsem/Chr15.genome.sorted.bam'

    out_dir = '/data1/shunfu1/SNPCalling/SimSNPs_MultiRun_5/reads_N100K_L100_Err0.00/samtools/'; run_cmd('mkdir -p %s'%out_dir)
    mpileup_file = '%s/res.mpileup'%out_dir
    bcf_file = '%s/res.bcf'%out_dir


    #cmd = 'samtools mpileup -f %s -l %s --output-MQ %s > %s'%\
    #      (ref_file, bed_file, bam_file, mpileup_file)
    #cmd = 'samtools mpileup -f %s %s > %s'%\
    #      (ref_file, bam_file, mpileup_file)
    #pdb.set_trace()
    #run_cmd(cmd)

    #cmd = 'bcftools call --threads 20 -o %s %s'%\
    #      (bcf_file, mpileup_file)
    cmd = 'bcftools mpileup -Ou -f %s %s | bcftools call -Ou -mv | bcftools filter -s LowQual > %s'%\
          (ref_file, bam_file, bcf_file)
    pdb.set_trace()
    run_cmd(cmd)

    return

if __name__ == "__main__":

    main()