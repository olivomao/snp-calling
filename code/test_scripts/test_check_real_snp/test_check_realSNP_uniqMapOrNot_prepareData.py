import pdb
from old_code.util import run_cmd
from tool_address import FA2FQ_address

'''
 prepare data necessary to run abSNP or GATK pipeline

 especially we want to prepare:
 - snps m/p: real, filtered so that all snps are guaranteed to be covered by some read
 - target genome m/p: snps applied onto reference genome
 - reads m/p fa/fq: guaranteed to cover all snps, m/p merged
'''

#### Configuration ####

RNASeqReadSimulatorPath='tools/RNASeqReadSimulator/' 

old_code_Dir = 'old_code/'

srcDir = '/data1/shunfu1/SNPCalling/data_real_HeteSNPsCodingRegion/'

chrom = 'chr15'

#input
trBed = '%s/gencode.v19.annotation.target.chr15.sorted.bed'%(srcDir)
refGenome = '%s/%s.fa'%(srcDir, chrom)

snp_files = ['%s/heteSNPs/chr15_snp_m_codingRegion_noHomo.txt'%srcDir,
             '%s/heteSNPs/chr15_snp_p_codingRegion_noHomo.txt'%srcDir]

#read bed files (snps mostly covered)
outBeds = ['%s/heteSNPs/read_m_SNPs_covered.bed'%srcDir,
           '%s/heteSNPs/read_p_SNPs_covered.bed'%srcDir]

outReadFa = ['%s/heteSNPs/read_m_SNPs_covered.fa'%srcDir,
             '%s/heteSNPs/read_p_SNPs_covered.fa'%srcDir]

outReadFq = ['%s/heteSNPs/read_m_SNPs_covered.fq'%srcDir,
             '%s/heteSNPs/read_p_SNPs_covered.fq'%srcDir]
outReadFq_merged = '%s/heteSNPs/read_merged_SNPs_covered.fq'%srcDir

readLen = 100
errRate = 0.0

snpReadCovFile = ['%s/heteSNPs/read_m_SNPs_covered.snp_read_cov'%srcDir,
                  '%s/heteSNPs/read_p_SNPs_covered.snp_read_cov'%srcDir]

filt_snp_files = ['%s/heteSNPs/chr15_snp_m_codingRegion_noHomo_allCoveredByReads.txt'%srcDir,
                  '%s/heteSNPs/chr15_snp_p_codingRegion_noHomo_allCoveredByReads.txt'%srcDir]

target_genome_filt_snp = ['%s/heteSNPs/chr15_m_HeteFiltSNPs.fa'%srcDir,
                          '%s/heteSNPs/chr15_p_HeteFiltSNPs.fa'%srcDir]

#### Simulation ####

for i in range(2):
    cmd = 'python sim_data_generator.py --readBed_generation_at_sel_snps '+\
                                       '-s %s '%snp_files[i]+\
                                       '-b %s '%trBed+\
                                       '-o %s '%outBeds[i]+\
                                       '-L %d '%readLen+\
                                       '-c %s '%chrom #'-o2 log '
    #pdb.set_trace()
    run_cmd(cmd)

for i in range(2):
    cmd = 'python sim_data_generator.py --snp_read_cov1 '+\
                                       '-s %s '%snp_files[i]+\
                                       '-m %s '%outBeds[0]+\
                                       '-p %s '%outBeds[1]+\
                                       '-o %s'%snpReadCovFile[i]
    #pdb.set_trace()
    run_cmd(cmd)

for i in range(2):
    cmd = 'python sim_data_generator.py --sel_snps_covered '+\
                                       '-i %s '%snp_files[i]+\
                                       '-o %s '%filt_snp_files[i]+\
                                       '-a %s '%snpReadCovFile[i]
    #pdb.set_trace()
    run_cmd(cmd)

for i in range(2):
    cmd = 'python %s/batch_run_realData.py --gen_genome_of_snps '%old_code_Dir+\
                                       '-r %s '%refGenome+\
                                       '-c %s '%chrom+\
                                       '-s %s '%filt_snp_files[i]+\
                                       '-t %s'%target_genome_filt_snp[i]
    #pdb.set_trace()
    run_cmd(cmd)

for i in range(2):
    cmd = 'python %s/getseqfrombed.py -r %.2f -l %d %s %s > %s'% \
          (RNASeqReadSimulatorPath, errRate, readLen, outBeds[i], target_genome_filt_snp[i], outReadFa[i])
    #pdb.set_trace()
    run_cmd(cmd)

#fa to fq
for i in range(2):
    cmd = 'perl ' + FA2FQ_address + ' ' + outReadFa[i] + ' > '  + outReadFq[i]
    #pdb.set_trace()
    run_cmd(cmd)

#merge reads
cmd = 'python sim_data_generator.py --merge_reads -m %s -p %s -o %s --uniqID'%(outReadFq[0], outReadFq[1], outReadFq_merged)
#pdb.set_trace()
run_cmd(cmd)