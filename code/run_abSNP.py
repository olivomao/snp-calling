# a wrapper for external users using abSNP from Github

import sys, pdb
from old_code.util import run_cmd

'''
python run_abSNP.py --snp_simulator
                    --bedSorted <bedSorted>
                    --outDir <outDir>
                    --qt <quantile>
                    --refGenome <refGenomeFile>
                    --NumSNP <numSNP>

Command descriptions:
--bedSorted <bedSorted>
  absolute path of a sorted bed file (e.g. a reference transcriptome) to specify where SNPs shall be generated
  to sort a bed file (e.g by its chrom & genome start pos), you can use e.g. "sort -k 1,1 -k 2,2n in.BED > out.BED".
--outDir <outDir>
  absolute path of the directory to store result files, including:
  exp_m.txt, exp_p.txt, 
  cov_m.txt, cov_p.txt, 
  cov_m_<qt>.txt, cov_p_<qt>.txt
  SNP_m.txt, SNP_p.txt
  Tar_m.txt, Tar_p.txt
  m, p stand for maternal & paternal alleles
--qt <quantile>
  an integer value in [0,99]; e.g. --qt 90 means SNPs are generated in top (100-90)% highly expressed genes
--refGenome <refGenomeFile>
  absolute path of the reference genome containing one chromosome
--NumSNP <numSNP>
  an integer value in [1, inf) to specify the number of SNPs generated per allele

Description:
Given a reference transcriptome (in bed format), assign rand expression levels for each gene independently for each maternal and paternal allele.
Then random snps restricted to these regions are generated, in genes that are top (100-quantile)% highly expressed for each allele.
The genomes of target alleles are also generated.

Output files:
-- exp_m.txt, exp_p.txt
   expression levels of RNA transcript for maternal  & paternal alleles 
- SNP_m.txt, SNP_p.txt
   simulated snps, in format of (0-based locus, reference base, '-->', target/mutated base)
- Tar_m.txt, Tar_p.txt
   genomes of target alleles that contain the simulated SNPs
[removed files]:
- cov_m.txt, cov_p.txt
   expression coverage per locus for maternal & paternal alleles
- cov_m_<qt>.txt, cov_p_<qt>.txt;
   expression coverage per locus (only the loci where corresponding transcripts are top (100-qt)% highly expressed) for maternal & paternal alleles 
'''
def snp_simulator(args):
    from sim_data_generator import gen_sim_exp_snp_tar

    bedSorted = args[args.index('--bedSorted')+1]
    outDir = args[args.index('--outDir')+1]; run_cmd('mkdir -p %s'%outDir)
    qt = int(args[args.index('--qt')+1])
    refGenome = args[args.index('--refGenome')+1]
    numSNP = int(args[args.index('--NumSNP')+1])

    args2 = '--gen_sim_exp_snp_tar --bedSorted %s -O %s --qt %d --refGenome %s --NumSNP %d'%(bedSorted, outDir, qt, refGenome, numSNP)

    gen_sim_exp_snp_tar(args2.split())

    ### clear
    run_cmd('rm %s/coverage*'%outDir)

    return

'''
usage
python run_abSNP.py --read_simulator
                    --tarGenome_m <tarGenome_m>
                    --tarGenome_p <tarGenome_p>
                    --exp_m <exp_m>
                    --exp_p <exp_p>
                    --bedSorted <bedSorted>
                    --outDir <outDir>
                    --numReads <numReads>
                    --readLength <readLength>
                    --errRate <errRate>

Description:
Sample reads from the transcriptome of target maternal and paternal alleles. 
The amount of reads sampled from per RNA transcript is specified by expression files. 
The number of reads per allele, the read length and error rate are also required.

Input
--tarGenome_m <tarGenome_m>
    path of the target maternal genome, containing one chromosome, in fasta format. Generated from SNP simulator.
--tarGenome_p <tarGenome_p>
    path of the target paternal genome, containing one chromosome, in fasta format. Generated from SNP simulator.
--exp_m <exp_m>
    path of the expression level file of RNA transcript for maternal allele. Generated from SNP simulator.
--exp_p <exp_p>
    path of the expression level file of RNA transcript for paternal allele. Generated from SNP simulator.
--bedSorted <bedSorted>
    path of the reference transcriptome in sorted bed format. To be consistent with the one used in SNP simulator.
--outDir <outDir>
    path of the directory to store the generated simulated reads in fastq format (e.g. merged_reads.fq).
--numReads <numReads>
    integer value. Specify the number of reads to be generated per target allele.
--readLength <readLength>
    integer value. Specify the read length (e.g. 100).
--errRate <errRate>
    float value. Specify the read error rate (e.g. 0.01)

Output

merged_reads.fq
the generated simulated reads in fastq format, including reads sampled from both alleles.

'''
def read_simulator(args):

    genomeFile = []
    genomeFile.append(args[args.index('--tarGenome_m')+1])
    genomeFile.append(args[args.index('--tarGenome_p')+1])
    exp_path = []
    exp_path.append(args[args.index('--exp_m')+1])
    exp_path.append(args[args.index('--exp_p')+1])
    trBedFile = args[args.index('--bedSorted')+1]
    outDir = args[args.index('--outDir')+1]; run_cmd('mkdir -p %s'%outDir)
    numReads = int(args[args.index('--numReads')+1])
    L = int(args[args.index('--readLength')+1])
    errRate = float(args[args.index('--errRate')+1])
    suf = ['_m', '_p']

    #read generation
    for i in range(2):

        cmd = 'python sim_data_generator.py --read_generation '+\
              '-g %s -b %s -O %s '%(genomeFile[i], trBedFile, outDir)+\
              '-n %d -l %d -r %f '%(numReads, L, errRate)+\
              '--suffix %s --toFq '%suf[i]+\
              '--exp_path %s'%(exp_path[i])
        #pdb.set_trace()
        run_cmd(cmd)

    #merge reads
    cmd = 'python sim_data_generator.py --merge_reads -m %s/reads_m.fq -p %s/reads_p.fq -o %s/merged_reads.fq --uniqID'% \
          (outDir, outDir, outDir)
    #pdb.set_trace()
    run_cmd(cmd)

    #clear
    run_cmd('rm -r %s/intermediate'%outDir)
    run_cmd('rm %s/reads*.*'%outDir)
    run_cmd('rm %s/*.log'%outDir)

    return

'''
usage:

description:
python run_abSNP.py --call_snps
                    --refGenome <refGenome>
                    --readFile <readFile>
                    --gtfFile <gtfFile>
                    --bedSorted <bedSorted>
                    --p <numProcess>
                    --readLength <readLength>
                    --alpha <alpha>
                    --outDir <outDir>

Input

--refGenome <refGenome>
    path of the reference genome containing one chromosome
--readFile <readFile>
    path of the sampled read file in fastq format
--gtfFile <gtfFile>
    path of the reference transcriptome (set of RNA transcripts, which are considered as target regions to detect SNPs) in gtf format
--bedSorted <bedSorted>
    path of the reference transcriptome in sorted bed format. To be consistent with the gtfFile.
--p <numProcess>
    an integer value specifyng how many CPU processes to use. Default 1. Recommend to use 10 to 20.
--readLength <readLength>
    integer value. Specify the read length (e.g. 100).
--alpha <alpha>
    a float value between 0 and 1. 0 corresponds to minimum false positive and 1 corresponds to max sensitivity.
--outDir <outDir> 
    path of the directory to store the file of called SNP candidates.

Output

snp_candidates.txt
    SNPs called by abSNP. format, with tab separated columns. col-0: 0-based genome locus; col-1: reference base; col-2: '\-->' (a dummy symbol); col-3: target base

'''
def call_snps(args):

    refGenome = args[args.index('--refGenome')+1]
    readFile = args[args.index('--readFile')+1]
    gtfFile = args[args.index('--gtfFile')+1]
    trBedFile = args[args.index('--bedSorted')+1]
    num_p = int(args[args.index('--p')+1])
    L = int(args[args.index('--readLength')+1])
    rocT = float(args[args.index('--alpha')+1])

    outDir = args[args.index('--outDir')+1];run_cmd('mkdir -p %s'%outDir)

    #### auto config    

    rsem_out_dir = '%s/rsem/'%outDir; run_cmd('mkdir -p %s'%rsem_out_dir)
    rsemPrefixName = 'chrom'
    rsemGenomePrefix = '%s/%s'%(rsem_out_dir, rsemPrefixName)
    rsemIsoformsResults = '%s/%s.isoforms.results'%(rsem_out_dir, rsemPrefixName)
    exon_address = '%s/exon.txt'%rsem_out_dir #intermediate output
    rsemCov = '%s/rsemCoverage.txt'%rsem_out_dir #intermediate output   
    tbam = '%s.transcript.bam'%rsemGenomePrefix
    gbam = '%s.genome.bam'%rsemGenomePrefix
    gsam = '%s.genome.sorted_n.sam'%rsemGenomePrefix 

    snp_out_dir = '%s/snp_res/'%outDir; run_cmd('mkdir -p %s'%snp_out_dir)
    samFn = '%s.genome.sorted_n'%rsemPrefixName

    callerThre=1 
    filt3_sameAb=1
    clear_intFiles=1
    code_folder = 'old_code/' 

    #quantification
    cmd = 'python quantification.py --RSEM1 --ref %s -r %s -g %s -O %s -p %s --num_p %d'% \
          (refGenome, readFile, gtfFile, rsem_out_dir, rsemPrefixName, num_p)
    run_cmd(cmd)#;pdb.set_trace()
    
    cmd = 'python quantification.py --rsemCoverage -i %s -b %s -e %s -c %s -L %d'% \
          (rsemIsoformsResults, trBedFile, exon_address, rsemCov, L)
    run_cmd(cmd)#;pdb.set_trace()

    cmd = 'python quantification.py --tbam2gbam -p %s -t %s -g %s'% \
          (rsemGenomePrefix, tbam, gbam)
    #pdb.set_trace()
    run_cmd(cmd)#;pdb.set_trace()

    #snp call
    cmd = 'python snp_call_para.py  -r %s '%refGenome+\
                                   '-c %s '%rsemCov+\
                                   '-s %s '%gsam+\
                                   '-O %s '%snp_out_dir+\
                                   '-T %d '%callerThre+\
                                   '-p %d '%num_p+\
                                   '-q %d '%num_p+\
                                   '--dupRun 0 '+\
                                   '--code_folder %s '%code_folder+\
                                   '--filt3_sameAb %d '%filt3_sameAb+\
                                   '--filt3_rocT %.2f '%rocT+\
                                   '--clear_intFiles %d'%clear_intFiles
    #pdb.set_trace()
    run_cmd(cmd)#;pdb.set_trace()

    #clear
    snp_res_filt_file = '%s/caller_output_snp_found_%s_T%d_para_rocT_%.2f.txt'%(snp_out_dir, samFn, callerThre, rocT)
    snp_res_filt_file_new = '%s/snp_candidates.txt'%(outDir)
    run_cmd('mv %s %s'%(snp_res_filt_file, snp_res_filt_file_new))
    run_cmd('rm -r %s'%rsem_out_dir)
    run_cmd('rm -r %s'%snp_out_dir)

    return

'''
python run_abSNP.py --check
                    --snp_res <snpRes>
                    --snp_m <snpM>
                    --snp_p <snpP>

description:
We provde a simple and quick check of SNP calling performance in terms of mis-detection and false positive.

input:
--snp_res <snpRes> path of SNP candidates called by abSNP
--snp_m <snpM\> path of ground truth SNPs from target maternal individual
--snp_p <snpP\> path of ground truth SNPs from target paternal individual

output:
Statistics of mis-detection and false positives will be printed on screen.

'''
def brief_check(args):

    from old_code.snp_res_statistics import get_snp_res_statistics, group_snps

    snp_res_address = args[args.index('--snp_res')+1]
    snp_m_address = args[args.index('--snp_m')+1]
    snp_p_address = args[args.index('--snp_p')+1]
    
    snps_list_sorted = get_snp_res_statistics(snp_res_address,
                                              snp_m_address,
                                              snp_p_address)
                                              
    [m_snps_cd, m_snps_md, p_snps_cd, p_snps_md, r_snps_fp] = group_snps(snps_list_sorted)
    
    print('# of mis-detection (m & p):%d'%(len(m_snps_md)+len(p_snps_md)))
    print('...# of mis-detection (m):%d'%len(m_snps_md))
    print('...# of mis-detection (p):%d'%len(p_snps_md))
    
    print('# of false-positive:%d'%len(r_snps_fp))
    
    print('# of correct-detection (m & p):%d'%(len(m_snps_cd)+len(p_snps_cd)))
    print('...# of correct-detection (m):%d'%len(m_snps_cd))
    print('...# of correct-detection (p):%d'%len(p_snps_cd))
    
    #pdb.set_trace()

    return

def main(args):

    if '--snp_simulator' in args:
        snp_simulator(args)
    elif '--read_simulator' in args:
        read_simulator(args)
    elif '--call_snps' in args:
        call_snps(args)
    elif '--check' in args:
        brief_check(args)
    else:
        print('unknown mode for run_abSNP. check usage from:')
        run_abSNP_usage()

    return

def run_abSNP_usage():

    usage_str = 'README or https://github.com/shunfumao/abSNP'

    print(usage_str)

    return

if __name__ == "__main__":
    
    args = sys.argv

    main(args)