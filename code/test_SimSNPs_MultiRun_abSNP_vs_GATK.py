import sys
from old_code.util import run_cmd
import pdb

#compare abSNP and GATK multiple times
#each time, new snps/tar generated, based on this, different read configured and sampled, and abSNP (test_abSNP_batch) and GATK (test_gatk_batch) run
def main(args):

    #pdb.set_trace()

    #configuration parameters
    reads = [] #list of [readLen, read err, read number, read number description]
    #reads.append([100, 0.00, 100000, '100K'])
    #reads.append([100, 0.01, 100000, '100K'])
    #reads.append([100, 0.00, 1000000, '1m'])
    #reads.append([100, 0.01, 1000000, '1m'])
    #reads.append([100, 0.00, 10000000, '10m'])
    reads.append([100, 0.01, 10000000, '10m'])
    #N_runs = 1
    N_run_snps_generated = [0,1,2,3,4,5,6,7,8,9] #in these folders, no need to generate snps/tar any more
    
    N_run_stt = 5
    N_run_stp = 9 #inclusive
    
    qt = 90

    NumSNP = 1000

    #default parameters
    SrcDir = '/data1/shunfu1/SNPCalling/data_large_0_idealCov/'
    RootFolder = '/data1/shunfu1/SNPCalling/'

    refGenome = '%s/Chr15.fa'%SrcDir
    bedSorted = '%s/hg19_chr15-UCSC-sorted.bed'%SrcDir

    #auto configuration
    #RootFolder/SimSNPs_MultiRun_<i>/expression & snps & tar genome files & cov files
    #                               /<readsLabel>/
    #                                            /intermediate/ read bed files
    #                                            /read files
    #                                            /rsem/ quantification files
    #                                            /gatk/ gatk files
    #                                            /snp_res/ snp res files
    #                                            /snp_res_gatk/ snp res (vs gatk) files

    for ith_run in xrange(N_run_stt, N_run_stp+1):

        RunFolder = '%s/SimSNPs_MultiRun_%d/'%(RootFolder, ith_run)
        cmd = 'mkdir -p %s/'%(RunFolder)
        #pdb.set_trace()
        run_cmd(cmd)

        #SNP & tar generation
        cmd = 'python sim_data_generator.py '+\
                      '--gen_sim_exp_snp_tar '+\
                      '--bedSorted %s '%bedSorted+\
                      '-O %s '%RunFolder+\
                      '--qt %d '%qt+\
                      '--refGenome %s '%refGenome+\
                      '--NumSNP %d'%NumSNP
        if ith_run not in N_run_snps_generated:
            run_cmd(cmd)
            N_run_snps_generated.append(ith_run)
        #pdb.set_trace()

        for read_config in reads:
            readLen = read_config[0]
            readErr = read_config[1]
            readNum = read_config[2]
            readNumDescription = read_config[3]

            readsLabel = 'reads_N%s_L%s_Err%.2f'%(readNumDescription, readLen, readErr)        

            #run abSNP
            cmd = 'python test_abSNP_batch.py '+\
                          '--srcDir %s '%SrcDir+\
                          '--dstDir %s '%RunFolder+\
                          '--numReads %d '%readNum+\
                          '--ReadLen  %d '%readLen+\
                          '--errRate %f '%readErr+\
                          '--readsLabel %s '%readsLabel
            #pdb.set_trace()
            run_cmd(cmd)

            #run GATK
            cmd = 'python test_gatk_batch.py '+\
                          '--srcDir %s '%SrcDir+\
                          '--dstDir %s '%RunFolder+\
                          '--numReads %d '%readNum+\
                          '--ReadLen  %d '%readLen+\
                          '--errRate %f '%readErr+\
                          '--readsLabel %s '%readsLabel
            #pdb.set_trace()
            run_cmd(cmd)

    print('current N_run_snps_generated: %s'%str(N_run_snps_generated))

    return

if __name__ == "__main__":
    
    args = sys.argv

    main(args)
