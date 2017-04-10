import pdb
from old_code.util import run_cmd

#use prev data (e.g. fa, snps, exp, target alleles), to generate new reads and do snp calls by abSNP code

def main():

    #snp1k_reads10m_abSNP_batch
    '''
    srcDir = '/data1/shunfu1/SNPCalling/data_large_0_idealCov/'
    dstDir = '/data1/shunfu1/SNPCalling/data_large_0_sbSNP_batch/'
    numReads = 10000000
    L = 100
    errRate = 0.01
    readsLabel = 'reads_N%s_L%s_Err%.2f'%(numReads, L, errRate)
    '''
    
    #snp20_reads100k_10_abSNP_batch
    #'''
    srcDir = '/data1/shunfu1/SNPCalling/snp20_reads100k_10/'
    dstDir = '/data1/shunfu1/SNPCalling/snp20_reads100k_10/'

    numReads = 100000
    L = 100
    errRate = 0.00
    readsLabel = '' #'reads_N%s_L%s_Err%.2f'%(numReads, L, errRate)
    #'''
    
    chrom = 'Chr15'

    refGenome = '%s/%s.fa'%(srcDir, chrom)

    genomeFile = ['%s/Tar_m.txt'%srcDir,
                  '%s/Tar_p.txt'%srcDir]

    trBedFile = '%s/hg19_chr15-UCSC-sorted.bed'%srcDir
    gtfFile = '%s/hg19_chr15-UCSC.gtf'%srcDir     

    read_out_dir = '%s/%s/'%(dstDir, readsLabel)
    rsem_out_dir = '%s/rsem/'%(dstDir)
    snp_out_dir = '%s/snp_res_abSNP_[filtSameAb1]/'%(dstDir)

    exon_address = '%s/exon.txt'%dstDir #intermediate output
    rsemCov = '%s/count_rsem.txt'%dstDir #intermediate output   

    rsemPrefixName = chrom
    rsemGenomePrefix = '%s/%s'%(rsem_out_dir, rsemPrefixName)
    tbam = '%s.transcript.bam'%rsemGenomePrefix
    gbam = '%s.genome.bam'%rsemGenomePrefix
    gsam = '%s.genome.sorted_n.sam'%rsemGenomePrefix    

    suf = ['_m', '_p']
    exp_path = ['%s/exp_m.txt'%srcDir, 
                '%s/exp_p.txt'%srcDir]

    true_SNP = ['%s/SNP_m.txt'%srcDir, 
                '%s/SNP_p.txt'%srcDir]

    samFn = 'Chr15.genome.sorted_n'
    cntFile = '%s/rsem/split_sorted_sam_%s/count_y/count_%s.txt'%(dstDir, samFn, samFn)
    cntFileAltInfo = '%s/rsem/split_sorted_sam_%s/count_y_altInfo/count_%s_altInfo.txt'%(dstDir, samFn, samFn)

    snp_res_0 = '%s/caller_output_snp_found_Chr15.genome.sorted_n_T1_para.txt'%snp_out_dir
    snp_res_1 = '%s/caller_output_snp_found_Chr15.genome.sorted_n_T1_para_rocT_1.00.txt'%snp_out_dir

    snpLog = '%s/snpLog.txt'%snp_out_dir
    snpSum = '%s/snpSum.txt'%snp_out_dir
    snpSum2 = '%s/snpSum2.txt'%snp_out_dir

    pdb.set_trace()

    #read generation
    for i in range(2):

        cmd = 'python sim_data_generator.py --read_generation '+\
              '-g %s -b %s -O %s '%(genomeFile[i], trBedFile, read_out_dir)+\
              '-n %d -l %d -r %f '%(numReads, L, errRate)+\
              '--suffix %s --toFq '%suf[i]+\
              '--exp_path %s'%(exp_path[i])
        #pdb.set_trace()
        #run_cmd(cmd)

    #merge reads
    cmd = 'python sim_data_generator.py --merge_reads -m %s/reads_m.fq -p %s/reads_p.fq -o %s/merged_reads.fq --uniqID'% \
          (read_out_dir, read_out_dir, read_out_dir)
    #pdb.set_trace()
    #run_cmd(cmd)

    #quantification

    cmd = 'python quantification.py --RSEM1 --ref %s -r %s/merged_reads.fq -g %s -O %s -p %s'% \
          (refGenome, read_out_dir, gtfFile, rsem_out_dir, rsemPrefixName)
    #pdb.set_trace()
    #run_cmd(cmd)

    rsemIsoformsResults = '%s/%s.isoforms.results'%(rsem_out_dir, rsemPrefixName)

    cmd = 'python quantification.py --rsemCoverage -i %s -b %s -e %s -c %s -L %d'% \
          (rsemIsoformsResults, trBedFile, exon_address, rsemCov, L)
    #pdb.set_trace()
    #run_cmd(cmd)

    cmd = 'python quantification.py --tbam2gbam -p %s -t %s -g %s'% \
          (rsemGenomePrefix, tbam, gbam)
    #pdb.set_trace()
    #run_cmd(cmd)

    #snp call
    #pdb.set_trace()
    cmd = 'python snp_call_para.py  -r %s '%refGenome+\
                                   '-c %s '%rsemCov+\
                                   '-s %s '%gsam+\
                                   '-O %s '%snp_out_dir+\
                                   '-T 1 '+\
                                   '-p 20 '+\
                                   '--dupRun 0 '+\
                                   '--code_folder old_code '+\
                                   '--filt3_sameAb 1 '+\
                                   '--filt3_rocT 0,1'
    #run_cmd(cmd)

    #evaluator
    pdb.set_trace()

    cmd = 'python evaluator.py --loadSnpInfo '+\
                                  '-L1 m -F1 %s '%true_SNP[0]+\
                                  '-L2 p -F2 %s '%true_SNP[1]+\
                                  '-L3 O -F3 %s '%snp_res_0+\
                                  '-L4 G -F4 %s '%snp_res_1+\
                                  '-C1 %s '%cntFile+\
                                  '-C2 %s '%cntFileAltInfo+\
                                  '--snpLog %s '%snpLog+\
                                  '--snpSum %s '%snpSum+\
                                  '--snpSum2 %s '%snpSum2
    run_cmd(cmd)

    return

if __name__ == "__main__":

    main()