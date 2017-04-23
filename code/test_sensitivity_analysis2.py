import sys, os, pdb

from old_code.util import run_cmd, run_cmds

def test_sensitivity_analysis_batch(args):

    #configuration (e.g. args)
    readsLabel = 'reads_N100K_L100_Err0.00'
    #readsLabel = 'reads_N100K_L100_Err0.01'
    #readsLabel = 'reads_N1m_L100_Err0.01'
    #readsLabel = 'reads_N10m_L100_Err0.01'
    run_stt = 5
    run_stp = 5 #inclusive

    groupValCalcOption = 2     # 0 (based on m cov/ p cov) 1 (based on ab) 2 (based on gatk/rsem multimapping)
    groupValsStr = '1,2,5,10,100' #'1,2,3,4,5,6,7,8,9,10,20' #mp cov '10,100,500,1000,1e10' #ab '10,100,1000,1e10' #map '1,2,5,10,100'
    
    old_code_dir = 'old_code/'
    srcDir = '/data1/shunfu1/SNPCalling/'

    T = 1 #caller threshold
    rocT=[0.00, 1.00]

    chrom = 'Chr15'
    refGenome = '/data1/shunfu1/SNPCalling/data_large_0_idealCov/%s.fa'%chrom

    num_p = 1 #num parallel

    snpSum3_folder_name = 'snpSum3_SnpLoc0Based'#'snpSum_3_debug_gen_count_sel' #'snpSum3_2' para applied

    #### 0. auto configuration

    run_files = {} #key - ith_run val - {} key - file description val - file path
    # file description & path:
    #
    cmds = []
    for rn in xrange(run_stt, run_stp+1):
        run_files[rn]={}

        runDirName = 'SimSNPs_MultiRun_%d'%rn
        runDir = '%s/%s/'%(srcDir, runDirName)

        run_files[rn]['runDir']=runDir
        run_files[rn]['snpSum3Dir']='%s/%s/%s/'%(runDir, readsLabel, snpSum3_folder_name); run_cmd('mkdir -p %s'%run_files[rn]['snpSum3Dir'])
        run_files[rn]['snpSum3File']='%s/%s/%s/snpSum3.txt'%(runDir, readsLabel, snpSum3_folder_name);

        run_files[rn]['s_m'] = '%s/SNP_m.txt'%(runDir)
        run_files[rn]['s_p'] = '%s/SNP_p.txt'%(runDir)
        run_files[rn]['s_ab_0'] = '%s/%s/snp_res/T%d/caller_output_snp_found_%s.genome.sorted_n_T%d_para_rocT_%.2f.txt'%(runDir, readsLabel, T, chrom, T, rocT[0])
        run_files[rn]['s_ab_1'] = '%s/%s/snp_res/T%d/caller_output_snp_found_%s.genome.sorted_n_T%d_para_rocT_%.2f.txt'%(runDir, readsLabel, T, chrom, T, rocT[1])
        run_files[rn]['s_gatk'] = '%s/%s/gatk/GATK_out/raw_variants.vcf.txt'%(runDir, readsLabel)

        run_files[rn]['sam_gatk'] = '%s/%s/gatk/2pass/dedupped.sorted_n.sam'%(runDir, readsLabel)
        run_files[rn]['bam_gatk'] = '%s/%s/gatk/2pass/dedupped.bam'%(runDir, readsLabel)
        if os.path.exists(run_files[rn]['sam_gatk'])==False:
            cmd = 'samtools sort -n -o %s -@ 20 %s'%(run_files[rn]['sam_gatk'], run_files[rn]['bam_gatk'])
            
            if num_p == 1:
                #pdb.set_trace()
                pass; #run_cmd(cmd)
            else:
                cmds.append(cmd)


        run_files[rn]['sam_rsem'] = '%s/%s/rsem/Chr15.genome.sorted_n.sam'%(runDir, readsLabel)

        run_files[rn]['cov'] = '%s/%s/rsem/rsemCoverage.txt'%(runDir, readsLabel)

        run_files[rn]['cnt_gatk_fn'] = 'count_deduppedSam.txt'
        run_files[rn]['cnt_alt_gatk'] = '%s/%s'%(run_files[rn]['snpSum3Dir'], 'count_deduppedSam_altInfo.txt')
        run_files[rn]['cnt_rsem_fn'] = 'count_rsemSam.txt' 
        run_files[rn]['cnt_alt_rsem'] = '%s/%s'%(run_files[rn]['snpSum3Dir'], 'count_rsemSam_altInfo.txt')

        run_files[rn]['m_rd_bed']='%s/%s/intermediate/reads_m.bed'%(runDir, readsLabel)
        run_files[rn]['p_rd_bed']='%s/%s/intermediate/reads_p.bed'%(runDir, readsLabel)

        run_files[rn]['sens_res']='%s/%s/%s/sensitivity.txt'%(runDir, readsLabel, snpSum3_folder_name);

    if num_p > 1:
        #pdb.set_trace()
        pass; #run_cmds(cmds, num_p)

    #### 1. selective count alt generation
    #pdb.set_trace()
    cmds = []
    for rn in xrange(run_stt, run_stp+1):

        #selective count gatk
        cmd =  'python %s/count_read_lambda.py --gen_count_selectively '%old_code_dir+\
                                       '-s0 %s '%run_files[rn]['s_m']+\
                                       '-s1 %s '%run_files[rn]['s_p']+\
                                       '-s2 %s '%run_files[rn]['s_ab_1']+\
                                       '-s3 %s '%run_files[rn]['s_gatk'] +\
                                       '--sam %s '%run_files[rn]['sam_gatk']+\
                                       '--isRsemSam 0 '+\
                                       '--ref %s '%refGenome+\
                                       '--cov %s '%run_files[rn]['cov']+\
                                       '--countDir %s '%run_files[rn]['snpSum3Dir']+\
                                       '--countFn %s'%run_files[rn]['cnt_gatk_fn']
        #pdb.set_trace()
        if num_p == 1:
            #pdb.set_trace()
            pass; #run_cmd(cmd)
        else:
            cmds.append(cmd)

        #selective count rsem
        cmd =  'python %s/count_read_lambda.py --gen_count_selectively '%old_code_dir+\
                                       '-s0 %s '%run_files[rn]['s_m']+\
                                       '-s1 %s '%run_files[rn]['s_p']+\
                                       '-s2 %s '%run_files[rn]['s_ab_1']+\
                                       '-s3 %s '%run_files[rn]['s_gatk'] +\
                                       '--sam %s '%run_files[rn]['sam_rsem']+\
                                       '--isRsemSam 1 '+\
                                       '--ref %s '%refGenome+\
                                       '--cov %s '%run_files[rn]['cov']+\
                                       '--countDir %s '%run_files[rn]['snpSum3Dir']+\
                                       '--countFn %s'%run_files[rn]['cnt_rsem_fn']
        #pdb.set_trace()
        if num_p == 1:
            #pdb.set_trace()
            pass; #run_cmd(cmd)
        else:
            cmds.append(cmd)
    
    if num_p > 1:
        #pdb.set_trace()
        pass #run_cmds(cmds, num_p)

    #### 2: gen snpSum3
    #pdb.set_trace()
    cmds = []
    for rn in xrange(run_stt, run_stp+1):
        #gen snpSum3       
        cmd = 'python evaluator.py --gen_snpSum3 '+\
                                   '-m %s '%run_files[rn]['s_m']+\
                                   '-p %s '%run_files[rn]['s_p']+\
                                   '--L1 %s --F1 %s '%('abSNP_a0', run_files[rn]['s_ab_0'])+\
                                   '--L2 %s --F2 %s '%('abSNP_a1', run_files[rn]['s_ab_1'])+\
                                   '--L3 %s --F3 %s '%('GATK', run_files[rn]['s_gatk'])+\
                                   '--caG %s '%run_files[rn]['cnt_alt_gatk']+\
                                   '--caR %s '%run_files[rn]['cnt_alt_rsem']+\
                                   '-O %s '%run_files[rn]['snpSum3Dir']+\
                                   '--rm %s '%run_files[rn]['m_rd_bed']+\
                                   '--rp %s '%run_files[rn]['p_rd_bed']+\
                                   '--removeIntermediateFiles'
        #pdb.set_trace()
        if num_p == 1:
            #pdb.set_trace()
            pass; #run_cmd(cmd)
        else:
            cmds.append(cmd)

    if num_p > 1:
        #pdb.set_trace()
        pass #run_cmds(cmds, num_p)

    #### 3: detailed sensitivity
    #pdb.set_trace()
    cmds = []
    for rn in xrange(run_stt, run_stp+1):
        #detailed sensitivity analysis
        cmd = 'python evaluator.py --sensitivity_analysis '+\
                                '--snpSum3 %s '%(run_files[rn]['snpSum3File'])+\
                                '--groupValCalcOption %d '%groupValCalcOption+\
                                '--groupVals %s '%groupValsStr+\
                                '-o %s'%run_files[rn]['sens_res']
        #pdb.set_trace()
        if num_p == 1:
            pdb.set_trace()
            run_cmd(cmd)
        else:
            cmds.append(cmd)

    if num_p > 1:
        #pdb.set_trace()
        run_cmds(cmds, num_p)

    #pooled fig view






'''
usage:

python test_sensitivity_analysis.py --batch
'''
if __name__ == "__main__":

    args = sys.argv

    if '--batch' in args:
        test_sensitivity_analysis_batch(args)