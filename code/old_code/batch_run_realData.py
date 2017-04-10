from util import *
from snp_oper import *
import sys, pdb, copy

def gen_rsem_quant_sam():

    '''
    1. rsem quantification
    '''
    dir = '/data1/shunfu1/SNPCalling/data/'

    ref_addr = '%s/hg19.fa'%dir 
    r1 = '%s/reads_large/reads_1.fastq'%dir
    r2 = '%s/reads_large/reads_2.fastq'%dir
    gtf = '%s/gene_annotation/gencode.v19.annotation.gtf'%dir

    cmd = 'python ReadProcess.py --RSEM2 --ref %s '%ref_addr+\
          '--r1 %s '%r1+\
          '--r2 %s '%r2+\
          '--gtf %s'%gtf
    run_cmd(cmd)

    '''
    2. prepare rsem sam
    '''
    rsem_dir = '%s/rsem/'%dir
    sample_name = 'exp'

    tbam = '%s/%s.transcript.bam'%(rsem_dir, sample_name)
    gbam = '%s/%s.genome.bam'%(rsem_dir, sample_name)
    osam = '%s/%s.genome.sorted_n.sam'%(rsem_dir, sample_name)

    cmd = 'rsem-tbam2gbam %s/hg19 '%rsem_dir +\
           '%s %s -p 20'%(tbam, gbam) 
    #run_cmd(cmd)

    cmd = 'samtools sort -n -@ 20 -o %s %s'%(osam, gbam)
    #run_cmd(cmd)

    return

'''
a wrapper to call sim_reads.py to convert gtf (transcriptome) to bed
only exons in gtf are considered, 

for example, convert
'/data1/shunfu1/SNPCalling/data/gene_annotation/gencode.v19.annotation.gtf' into 
chromosome based bed files (e.g. gencode.v19.annotation.target.chr15.sorted.bed) at same folder.

note:
- we find gencode chr15 has more transcripts (e.g. 7k+) than ucsc chr15 bed
- we don't deal with thickStart/ thickEnd in bed (i.e. the two entries in bed file contains dummy information)

'''
def gtf2bed():

    #tr_gtf = '/data1/shunfu1/SNPCalling/genome/hg19_chr15-UCSC.gtf'
    tr_gtf = '/data1/shunfu1/SNPCalling/data/gene_annotation/gencode.v19.annotation.gtf'
    #out_dir = '/data1/shunfu1/SNPCalling/genome/'
    out_dir = '/data1/shunfu1/SNPCalling/data/gene_annotation/'

    cmd = 'python /home/shunfu1/ref_shannon_modi/code_sgRefShannon/sim_reads.py --gtf2bed '+\
          '-i %s '%tr_gtf+\
          '-O %s '%out_dir

    run_cmd(cmd)

    return

def main():

    return

'''
usage:

#filtered vcf (e.g. chr15, snv only) --> target snps (e.g. SNP_m -- 1st allele & SNP_p -- 2nd allele)
python batch_run_realData.py --vcf2snp -i vcf_file -o path/to/<snp prefix>

#generate target genome (m or p) w/ extracted snp info
python batch_run_realData.py --gen_genome_of_snps -r ref_genome -c target_chr -s snp_file -t path/to/target_genome
'''

if __name__ == "__main__":

    args = sys.argv

    if '--vcf2snp' in args:
        vcf2snp(args)
    elif '--gen_genome_of_snps' in args:
        gen_genome_of_snps(args)
    else:
        main()

    
    