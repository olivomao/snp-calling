import os, pdb, sys
from old_code.util import run_cmd

'''
usage:

python test_gatk_batch.py
--srcDir srcDir
--dstDir dstDir
--numReads numReads
--ReadLen  L
--errRate errRate
--readsLabel readsLabel

output:

dstDir/readsLabel/gatk/<genome>
                      /1pass/
                      /<genome>_2pass/
                      /2pass/ (alignment files)
                      /GATK_out/ (vcf file and extracted snps)
'''
def main():

    args = sys.argv

    if '--srcDir' in args:
        srcDir = args[args.index('--srcDir')+1]
    else:
        srcDir = '/data1/shunfu1/SNPCalling/data_large_0_idealCov/'

    if '--dstDir' in args:
        dstDir = args[args.index('--dstDir')+1]
    else:
        dstDir = '/data1/shunfu1/SNPCalling/data_large_0_abSNP_batch/'

    num_p = 20

    if '--numReads' in args:
        numReads = int(args[args.index('--numReads')+1])
    else:
        numReads = 10000000 #100000 #10000000

    if '--ReadLen' in args:
        L = int(args[args.index('--ReadLen')+1])
    else:
        L = 100
    
    if '--errRate' in args:
        errRate = float(args[args.index('--errRate')+1])
    else:
        errRate = 0.01

    if '--readsLabel' in args:
        readsLabel = args[args.index('--readsLabel')+1]
    else:
        readsLabel = 'reads_N%s_L%s_Err%.2f'%(numReads, L, errRate)

    chrom = 'Chr15'    

    true_SNP = ['%s/SNP_m.txt'%dstDir, 
                '%s/SNP_p.txt'%dstDir]

    #info from abSNP
    callerThre=1
    rocT = [0.0, 1.0] # filt poisson

    samFn = '%s.genome.sorted_n'%chrom

    rsem_out_dir = '%s/%s/rsem/'%(dstDir, readsLabel)

    cntFile = '%s/split_sorted_sam_%s/count_y/count_%s.txt'%(rsem_out_dir, samFn, samFn)
    cntFileAltInfo = '%s/split_sorted_sam_%s/count_y_altInfo/count_%s_altInfo.txt'%(rsem_out_dir, samFn, samFn)

    #auto configuration
    DefRefPath = '%s/%s/'%(dstDir, readsLabel)
    GATKdir = 'gatk'    

    refGenome = '%s/%s.fa'%(srcDir, chrom)
    gtfFile = '%s/hg19_chr15-UCSC.gtf'%srcDir   

    readAddress = '%s/%s/merged_reads.fq'%(dstDir, readsLabel)

    toolFolder = os.getcwd()+'/tools/'    

    #snp call
    cmd = 'python snp_call_gatk.py '+\
          '--do_gatk_best_practice_modi '+\
          '--DefRefPath %s '%DefRefPath+\
          '--GATKdir %s '%GATKdir+\
          '--Genome %s '%refGenome+\
          '--GTF %s '%gtfFile+\
          '--Rd %s '%readAddress+\
          '--toolFolder %s '%toolFolder
    #pdb.set_trace()
    run_cmd(cmd)

    #evaluator
    for i in range(len(rocT)):
        snpResSubDir = '%s/%s/snp_res_gatk/vs_abSNP_callerT%d_filt3rocT%.2f/'%(dstDir, readsLabel, callerThre, rocT[i])
        run_cmd('mkdir -p %s'%snpResSubDir)

        snp_res_abSNP = '%s/%s/snp_res/T%d/caller_output_snp_found_%s_T%d_para_rocT_%.2f.txt'%(dstDir, readsLabel, callerThre ,samFn, callerThre, rocT[i])
        snp_res_GATK = '%s/%s/gatk/GATK_out/raw_variants.vcf.txt'%(dstDir, readsLabel)

        snpLog = '%s/snpLog.txt'%(snpResSubDir)
        snpSum = '%s/snpSum.txt'%(snpResSubDir)
        snpSum2 = '%s/snpSum2.txt'%(snpResSubDir)

        cmd = 'python evaluator.py --loadSnpInfo '+\
                                  '-L1 m -F1 %s '%true_SNP[0]+\
                                  '-L2 p -F2 %s '%true_SNP[1]+\
                                  '-L3 O -F3 %s '%snp_res_abSNP+\
                                  '-L4 G -F4 %s '%snp_res_GATK+\
                                  '-C1 %s '%cntFile+\
                                  '-C2 %s '%cntFileAltInfo+\
                                  '--snpLog %s '%snpLog+\
                                  '--snpSum %s '%snpSum+\
                                  '--snpSum2 %s '%snpSum2
        #pdb.set_trace()
        run_cmd(cmd)

    
    return

if __name__ == "__main__":

    main()