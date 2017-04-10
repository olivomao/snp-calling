import sys, pdb
from old_code.batch_run_case2_rsem import do_gatk_best_practice_modi
from old_code.snp_res_statistics import isVCF, extractVCF

'''
based on old_code-->batch_run_case2_rsem-->do_gatk_best_practice_modi,
more flexible on input parameter & instructures

usage:
--do_gatk_best_practice_modi
--DefRefPath drp #replace the one of Address.py
--GATKdir gd #default data_GATK only a folder name, put related alignments & res here, under drp
--Genome fa #e.g. address of ref genome, no mutation
[--SkipPrepareGenomeRefIdx] #skip 1pass STAR --runMode genomeGenerate if done previously
--GTF gtfAddr #may contain more annotations (e.g. hg19.gtf while chr15.fa used)
--Rd readAddress #SE used currently (to extend to --Rd1 --Rd2)
[--toolFolder toolFolder] #GATK tools' location; e.g. import os cwd = os.getcwd(); e.g. GATK's location: toolFolder+"/GATK/GenomeAnalysisTK.jar"

example:
python snp_call_gatk.py --do_gatk_best_practice_modi
                        --DefRefPath /data1/shunfu1/SNPCalling/data/
                        --GATKdir data_GATK_chr15_reads100k
                        --Genome /data1/shunfu1/SNPCalling/data/genome/chr15.fa
                        --GTF /data1/shunfu1/SNPCalling/data/gene_annotation/gencode.v19.annotation.gtf
                        --Rd /data1/shunfu1/SNPCalling/data/sim_reads_100k.fq

output:
<drp>/<gd>/<genome>
          /1pass/
          /<genome>_2pass/
          /2pass/ (alignment files)
          /GATK_out/ (vcf file and extracted snps)
'''

def main():

    args = sys.argv

    called_variants = do_gatk_best_practice_modi(args)[0]
    #called_variants='%s/%s/GATK_out/raw_variants.vcf'%(args[args.index('--DefRefPath')+1], args[args.index('--GATKdir')+1])

    #pdb.set_trace()
    if isVCF(called_variants):
        called_snps = extractVCF(called_variants) # x.vcf (contain snps and other variants) --> x.vcf.txt (contain snps)

    return

if __name__ == "__main__":

    main()