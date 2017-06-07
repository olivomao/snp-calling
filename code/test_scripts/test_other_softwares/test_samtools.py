import pdb, sys, os
from old_code.util import run_cmd, run_cmds, from_fasta
from intervaltree import Interval, IntervalTree
from old_code.snp_res_statistics import extractVCF

'''
usage:
python test_samtools.py --convert_samtools --snpIn snpInFile(vcf gz format)

description:
vcf to our snp res (snpInFile + del '.gz' + '.txt')
'''
def convert_samtools(args):

    #pdb.set_trace()

    resFile_a = args[args.index('--snpIn')+1]

    cmd = 'gunzip %s'%resFile_a
    run_cmd(cmd)

    resFile_b = resFile_a[:-3]

    res=extractVCF(resFile_b)

    return

'''
usage:
python test_samtools.py --run_samtools_batch

description:
run samtools for a batch of data/modes
'''
def run_samtools_batch(args):

    #configuration parameters
    reads = [] #list of [readLen, read err, read number, read number description]
    #reads.append([100, 0.00, 100000, '100K'])
    #reads.append([100, 0.01, 100000, '100K'])
    #reads.append([100, 0.00, 1000000, '1m'])
    reads.append([100, 0.01, 1000000, '1m'])
    #reads.append([100, 0.00, 10000000, '10m'])
    #reads.append([100, 0.01, 10000000, '10m'])
    
    N_run_stt = 5
    N_run_stp = 9 #inclusive

    sam2bam = 1
    run_samtools = 1
    convert_samtools = 1

    num_p = 10
    if num_p > 1:
        sam2bam_cmds = []
        run_samtools_cmds = []
        convert_samtools_cmds = []
    
    #default parameters
    SrcDir = '/data1/shunfu1/SNPCalling/data_large_0_idealCov/'
    RootFolder = '/data1/shunfu1/SNPCalling/'

    refGenome = '%s/Chr15.fa'%SrcDir
    
    #auto configuration
    #RootFolder/SimSNPs_MultiRun_<i>/expression & snps & tar genome files & cov files
    #                               /<readsLabel>/
    #                                            /intermediate/ read bed files
    #                                            /read files
    #                                            /rsem/ quantification files
    #                                            /gatk/ gatk files
    #                                                 /2pass/Aligned.out.sorted.bam (*)
    #                                            /snp_res/ snp res files
    #                                            /snp_res_gatk/ snp res (vs gatk) files
    #                                            /snp_res_samtools/ snp res files: snp_samtools.txt

    for ith_run in xrange(N_run_stt, N_run_stp+1):

        #pdb.set_trace()

        RunFolder = '%s/SimSNPs_MultiRun_%d/'%(RootFolder, ith_run)        
        
        for read_config in reads:
            readLen = read_config[0]
            readErr = read_config[1]
            readNum = read_config[2]
            readNumDescription = read_config[3]

            readsLabel = 'reads_N%s_L%s_Err%.2f'%(readNumDescription, readLen, readErr) 

            bamFile = '%s/%s/rsem/Chr15.genome.sorted.bam'%(RunFolder, readsLabel)
            if os.path.exists(bamFile)==False:
                samFile = '%s/%s/rsem/Chr15.genome.bam'%(RunFolder, readsLabel)
                cmd = 'samtools sort -@ 20 -o %s %s'%(bamFile, samFile)
                if sam2bam==1:
                    if num_p==1:
                        run_cmd(cmd)
                    else:
                        sam2bam_cmds.append(cmd)

            resFolder = '%s/%s/snp_res_samtools/'%(RunFolder, readsLabel)
            cmd = 'mkdir -p %s'%resFolder
            run_cmd(cmd)

            #run samtools
            resFile_a = '%s/res.rsem.vcf.gz'%(resFolder)
            resFile_b = '%s/res.rsem.vcf'%(resFolder) 
            resFile_c = '%s/res.rsem.vcf.txt'%(resFolder) 

            cmd = 'samtools mpileup -ugf %s %s | bcftools call -vmO z -o %s'%(refGenome, bamFile, resFile_a)
            if run_samtools==1:
                if num_p==1:
                    run_cmd(cmd)
                else:
                    run_samtools_cmds.append(cmd)

            cmd = 'python test_samtools.py --convert_samtools --snpIn %s'%resFile_a          
            if convert_samtools==1:
                if num_p==1:
                    run_cmd(cmd)
                else:
                    convert_samtools_cmds.append(cmd)

    #pdb.set_trace()
    if sam2bam==1 and num_p>1:
        run_cmds(sam2bam_cmds, num_p)
    if run_samtools==1 and num_p>1:
        run_cmds(run_samtools_cmds, num_p)
    if convert_samtools==1 and num_p>1:
        run_cmds(convert_samtools_cmds, num_p)

    return


if __name__ == "__main__":

    args = sys.argv

    if '--run_samtools_batch' in args:
        run_samtools_batch(args)
    elif '--convert_samtools' in args:
        convert_samtools(args)
    else:
        pdb.set_trace()



