import sys, os, pdb
import numpy as np

from old_code.util import run_cmd, run_cmds

def evaluate_performance_batch(args):

    #configuration (e.g. args)
    #readsLabel = 'reads_N100K_L100_Err0.00'
    #readsLabel = 'reads_N100K_L100_Err0.01'
    readsLabel = 'reads_N1m_L100_Err0.01'
    #readsLabel = 'reads_N10m_L100_Err0.01'

    run_stt = 5
    run_stp = 9 #inclusive

    #opt-0 based on snpSum3 m cov/ p cov
    '''
    groupValCalcOption = 0
    groupValsStr = '10,100,500,1000,1e10'
    '''

    #opt-1 based on snpSum3 ab
    '''
    groupValCalcOption = 1
    groupValsStr = '0,10,100,1e3,1e4,1e6,1e8,1e10'
    '''

    #opt-2 based on snpSum3 gatk/rsem altMap
    #    esp, col-7 for count file in format 1 (gatk) and col-8 for count file in format 1 (rsem), esp:
    #    -1: no snp reads
    #    0:  part or all of snp reads uniq mapped (part, not all of snp reads alt mapped)
    #    1+ (float val): all snp reads alt mapped, # of avg alt mapping (excluding self)
    #'''
    #groupValCalcOption = 2
    #groupValsStr = '-1,0,1,2,4,6,8,10,20'
    #'''
    
    #opt-3 based on snpSum3 ab's percentile [0,1)...[98,99),[99,100) (need to determine the percentile from multi runs first)
    #'''
    #groupValCalcOption = 3
    #tmp_ab_percentile = 'tmp/tmp_ab_percentile.txt'; ab_Str = '--ab 0 '; step_Str = '--step 1 '
    #groupValsStr = '0,25,50,75,100' #'0,20,40,60,80,100' #'0,10,20,30,40,50,60,70,80,90,100'
    #'''

    #opt-4 based on snpSum3 snp reads cov percentile e.g. [0,5)...[95,100) (need to determine the percentile from multi runs first)
    #'''
    groupValCalcOption = 4
    tmp_ab_percentile = 'tmp/tmp_ab_percentile.txt'; ab_Str = '--ab 1 '; step_Str = '--step 5'
    groupValsStr = '0,25,50,75,100' #'0,20,40,60,80,100' #'0,10,20,30,40,50,60,70,80,90,100'
    #'''

    old_code_dir = 'old_code/'
    srcDir = '/data1/shunfu1/SNPCalling/'

    T = 1 #caller threshold
    rocT=[0.00, 1.00]

    chrom = 'Chr15'
    refGenome = '/data1/shunfu1/SNPCalling/data_large_0_idealCov/%s.fa'%chrom

    num_p = 20 #num parallel

    cnt_format = 1 #0-e.g. C,I,lambda_bj,lambda_sum  1-e.g. bj,num_alt_mapping

    #snpSum3_folder_name = 'snpSum3_varscan_filt' #'snpSum3_SnpLoc0Based_CntFormat1'#'snpSum_3_debug_gen_count_sel' #'snpSum3_2' para applied
    snpSum3_folder_name = 'snpSum3_samtools'
    #sensitivity_res_name = 'sensitivity_opt%d.txt'%groupValCalcOption
    #sensitivity_res_name = 'sensitivity_opt3_SnpReadCov_MultiMapSNP_rsem_uniqmap.txt' #'sensitivity_opt2_rsemAltMap_le0.txt' #'sensitivity_opt2.txt' #'sensitivity_opt2_rsemAltMap_ge1.txt'
    #sensitivity_res_name = 'sensitivity_opt3_SnpReadCov_UniqMapSNP.txt'
    #sensitivity_res_name = 'sensitivity_opt3_SnpReadCov_MultSNP_rsem_ge1.txt'
    sensitivity_res_name = 'sensitivity_opt4.txt'

    #### config - teps to run
    do_bam2sam = 1
    do_gen_count_selectively = 1
    do_gen_snpSum3 = 1
    do_detailed_sens_analysis = 0


    #### 0. auto configuration

    run_files = {} #key - ith_run val - {} key - file description val - file path

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
        #run_files[rn]['s_ab_0'] = '%s/%s/snp_res/T%d/caller_output_snp_found_%s.genome.sorted_n_T%d_para_rocT_%.2f.txt'%(runDir, readsLabel, T, chrom, T, rocT[0])
        #run_files[rn]['s_ab_1'] = '%s/%s/snp_res/T%d/caller_output_snp_found_%s.genome.sorted_n_T%d_para_rocT_%.2f.txt'%(runDir, readsLabel, T, chrom, T, rocT[1])
        #run_files[rn]['s_varscan_0'] = '%s/%s/snp_res_varscan/snp_varscan_mode%d.txt'%(runDir, readsLabel, 0)
        #run_files[rn]['s_varscan_1'] = '%s/%s/snp_res_varscan/snp_varscan_mode%d.txt'%(runDir, readsLabel, 1)
        #run_files[rn]['s_varscan_0'] = '%s/%s/snp_res_varscan/snp_varscan_mode%d_filt.txt'%(runDir, readsLabel, 0)
        #run_files[rn]['s_varscan_1'] = '%s/%s/snp_res_varscan/snp_varscan_mode%d_filt.txt'%(runDir, readsLabel, 1)
        run_files[rn]['s_samtools'] = '%s/%s/snp_res_samtools/res.rsem.vcf.txt'%(runDir, readsLabel)
        
        run_files[rn]['s_gatk'] = '%s/%s/gatk/GATK_out/raw_variants.vcf.txt'%(runDir, readsLabel)

        run_files[rn]['sam_gatk'] = '%s/%s/gatk/2pass/dedupped.sorted_n.sam'%(runDir, readsLabel)
        run_files[rn]['bam_gatk'] = '%s/%s/gatk/2pass/dedupped.bam'%(runDir, readsLabel)
        if os.path.exists(run_files[rn]['sam_gatk'])==False:
            cmd = 'samtools sort -n -o %s -@ 20 %s'%(run_files[rn]['sam_gatk'], run_files[rn]['bam_gatk'])
            
            if num_p == 1:
                if do_bam2sam==1: run_cmd(cmd)
            else:
                cmds.append(cmd)


        run_files[rn]['sam_rsem'] = '%s/%s/rsem/Chr15.genome.sorted_n.sam'%(runDir, readsLabel)

        run_files[rn]['cov'] = '%s/%s/rsem/rsemCoverage.txt'%(runDir, readsLabel)

        run_files[rn]['cnt_gatk_fn'] = 'count_deduppedSam.txt'
        run_files[rn]['cnt_gatk'] = '%s/%s'%(run_files[rn]['snpSum3Dir'], 'count_deduppedSam.txt')
        run_files[rn]['cnt_alt_gatk'] = '%s/%s'%(run_files[rn]['snpSum3Dir'], 'count_deduppedSam_altInfo.txt')
        run_files[rn]['cnt_rsem_fn'] = 'count_rsemSam.txt' 
        run_files[rn]['cnt_rsem'] = '%s/%s'%(run_files[rn]['snpSum3Dir'], 'count_rsemSam.txt')
        run_files[rn]['cnt_alt_rsem'] = '%s/%s'%(run_files[rn]['snpSum3Dir'], 'count_rsemSam_altInfo.txt')

        run_files[rn]['m_rd_bed']='%s/%s/intermediate/reads_m.bed'%(runDir, readsLabel)
        run_files[rn]['p_rd_bed']='%s/%s/intermediate/reads_p.bed'%(runDir, readsLabel)

        run_files[rn]['sens_res']='%s/%s/%s/%s'%(runDir, readsLabel, snpSum3_folder_name, sensitivity_res_name);

    if num_p > 1:
        if do_bam2sam==1: run_cmds(cmds, num_p) #pass; #run_cmds(cmds, num_p)

    #### 1. selective count alt generation
    cmds = []
    for rn in xrange(run_stt, run_stp+1):

        #selective count gatk
        cmd =  'python %s/count_read_lambda.py --gen_count_selectively '%old_code_dir+\
                                       '--format %d '%(cnt_format)+\
                                       '-s0 %s '%run_files[rn]['s_m']+\
                                       '-s1 %s '%run_files[rn]['s_p']+\
                                       '-s2 %s '%run_files[rn]['s_samtools']+\
                                       '-s3 %s '%run_files[rn]['s_gatk'] +\
                                       '--sam %s '%run_files[rn]['sam_gatk']+\
                                       '--isRsemSam 0 '+\
                                       '--ref %s '%refGenome+\
                                       '--cov %s '%run_files[rn]['cov']+\
                                       '--countDir %s '%run_files[rn]['snpSum3Dir']+\
                                       '--countFn %s'%run_files[rn]['cnt_gatk_fn']
        if num_p == 1:
            if do_gen_count_selectively==1: run_cmd(cmd) #pass; #run_cmd(cmd)
        else:
            cmds.append(cmd)

        #selective count rsem
        cmd =  'python %s/count_read_lambda.py --gen_count_selectively '%old_code_dir+\
                                       '--format %d '%(cnt_format)+\
                                       '-s0 %s '%run_files[rn]['s_m']+\
                                       '-s1 %s '%run_files[rn]['s_p']+\
                                       '-s2 %s '%run_files[rn]['s_samtools']+\
                                       '-s3 %s '%run_files[rn]['s_gatk'] +\
                                       '--sam %s '%run_files[rn]['sam_rsem']+\
                                       '--isRsemSam 1 '+\
                                       '--ref %s '%refGenome+\
                                       '--cov %s '%run_files[rn]['cov']+\
                                       '--countDir %s '%run_files[rn]['snpSum3Dir']+\
                                       '--countFn %s'%run_files[rn]['cnt_rsem_fn']
        if num_p == 1:
            if do_gen_count_selectively==1: run_cmd(cmd) #pass; #run_cmd(cmd)
        else:
            cmds.append(cmd)
    
    if num_p > 1:
        if do_gen_count_selectively==1: run_cmds(cmds, num_p) #pass #run_cmds(cmds, num_p)

    #### 2: gen snpSum3
    cmds = []
    for rn in xrange(run_stt, run_stp+1):
        #gen snpSum3 use --caG --caR (count alt file)
        '''
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
        '''
        #gen snpSum3 use --cG1 --cR1 (count file fmt 1 to check snp reads alt map)
        '''
        cmd = 'python evaluator.py --gen_snpSum3 '+\
                                   '-m %s '%run_files[rn]['s_m']+\
                                   '-p %s '%run_files[rn]['s_p']+\
                                   '--L1 %s --F1 %s '%('varscan_0', run_files[rn]['s_varscan_0'])+\
                                   '--L2 %s --F2 %s '%('varscan_1', run_files[rn]['s_varscan_1'])+\
                                   '--L3 %s --F3 %s '%('GATK', run_files[rn]['s_gatk'])+\
                                   '--cG1 %s '%run_files[rn]['cnt_gatk']+\
                                   '--cR1 %s '%run_files[rn]['cnt_rsem']+\
                                   '-O %s '%run_files[rn]['snpSum3Dir']+\
                                   '--rm %s '%run_files[rn]['m_rd_bed']+\
                                   '--rp %s '%run_files[rn]['p_rd_bed']+\
                                   '--removeIntermediateFiles '+\
                                   '--md_fp'
        '''

        cmd = 'python evaluator.py --gen_snpSum3 '+\
                                   '-m %s '%run_files[rn]['s_m']+\
                                   '-p %s '%run_files[rn]['s_p']+\
                                   '--L1 %s --F1 %s '%('samtools', run_files[rn]['s_samtools'])+\
                                   '--L2 %s --F2 %s '%('GATK', run_files[rn]['s_gatk'])+\
                                   '--cG1 %s '%run_files[rn]['cnt_gatk']+\
                                   '--cR1 %s '%run_files[rn]['cnt_rsem']+\
                                   '-O %s '%run_files[rn]['snpSum3Dir']+\
                                   '--rm %s '%run_files[rn]['m_rd_bed']+\
                                   '--rp %s '%run_files[rn]['p_rd_bed']+\
                                   '--removeIntermediateFiles '+\
                                   '--md_fp'

        if num_p == 1:
            if do_gen_snpSum3==1:
                #pdb.set_trace()
                run_cmd(cmd) #pass; #run_cmd(cmd)
        else:
            cmds.append(cmd)

    if num_p > 1:
        if do_gen_snpSum3==1: run_cmds(cmds, num_p) #pass #run_cmds(cmds, num_p)

    
    #### 3: in case percentile group, check ab percentile from multi runs
    if groupValCalcOption==3 or groupValCalcOption==4:
        cmd = 'python evaluator.py --calcAbPercentile '
        snpSum3_idx = 0
        for rn in xrange(run_stt, run_stp+1):
            cmd += '-%d %s '%(snpSum3_idx, run_files[rn]['snpSum3File'])
            snpSum3_idx+=1
        cmd += '-o %s '%tmp_ab_percentile
        cmd += '%s %s'%(ab_Str, step_Str)
        #pdb.set_trace()
        run_cmd(cmd)
    #pdb.set_trace()

    #### 4: detailed sensitivity
    cmds = []
    for rn in xrange(run_stt, run_stp+1):
        #detailed sensitivity analysis
        if groupValCalcOption==3 or groupValCalcOption==4:
            ab_per_str = '--ab_percentile %s '%tmp_ab_percentile
        else:
            ab_per_str = ''

        cmd = 'python evaluator.py --sensitivity_analysis '+\
                                '--snpSum3 %s '%(run_files[rn]['snpSum3File'])+\
                                '--groupValCalcOption %d %s '%(groupValCalcOption, ab_per_str)+\
                                '--groupVals %s '%groupValsStr+\
                                '-o %s'%run_files[rn]['sens_res']
        if num_p == 1:
            if do_detailed_sens_analysis==1: run_cmd(cmd)
        else:
            cmds.append(cmd)

    if num_p > 1:
        if do_detailed_sens_analysis==1: run_cmds(cmds, num_p)

    return

'''
calc avg/stdev of a set of sens files. and output outDir/outFn.txt & outFn_stdev.txt

sens file format:
line-0: #snpSum3File
line-1: #groupValCalcOption
line-2: #caller   [g0,g1)  [g1,g2) ... [gN-1, gN)
line-3: tot       cnts    cnts        cnts
line-4: caller0   cnts    cnts        cnts
[line-5 etc]

python evaluate_performance_othertools.py --calc_sens_avg_stdev -O outDir --fn outFn (e.g. *.txt) [--normalize ref_c_id]

ref_c_id: reference caller id for normalization, default 0 ('Tot')
'''
def calc_sens_avg_stdev(args):

    #configuration (e.g. args)
    
    run_stt = 5
    run_stp = 9 #inclusive

    run_files = {} #key - run id val - file path
    for rf in xrange(run_stt, run_stp+1):
        #run_files[rf] = 'tmp/sens%d.txt'%rf
        #run_files[rf] = '/data1/shunfu1/SNPCalling//SimSNPs_MultiRun_%d/reads_N1m_L100_Err0.01/snpSum3_SnpLoc0Based_CntFormat1/sensitivity_opt2.txt'%(rf)
        #run_files[rf] = '/data1/shunfu1/SNPCalling//SimSNPs_MultiRun_%d/reads_N1m_L100_Err0.01/snpSum3_SnpLoc0Based_CntFormat1/sensitivity_opt2_rsemAltMap_le0.txt'%(rf)
        #run_files[rf] = '/data1/shunfu1/SNPCalling//SimSNPs_MultiRun_%d/reads_N1m_L100_Err0.01/snpSum3_SnpLoc0Based_CntFormat1/sensitivity_opt2_rsemAltMap_ge1.txt'%(rf)

        #run_files[rf] = '/data1/shunfu1/SNPCalling//SimSNPs_MultiRun_%d/reads_N1m_L100_Err0.01/snpSum3_SnpLoc0Based_CntFormat1/sensitivity_opt1.txt'%(rf)

        #run_files[rf] = '/data1/shunfu1/SNPCalling//SimSNPs_MultiRun_%d/reads_N1m_L100_Err0.01/snpSum3_SnpLoc0Based_CntFormat1/sensitivity_opt3.txt'%(rf)
        #run_files[rf] = '/data1/shunfu1/SNPCalling//SimSNPs_MultiRun_%d/reads_N1m_L100_Err0.01/snpSum3_SnpLoc0Based_CntFormat1/sensitivity_opt3_SnpReadCov.txt'%(rf) sensitivity_opt3_SnpReadCov_MultiMapSNP
        #run_files[rf] = '/data1/shunfu1/SNPCalling//SimSNPs_MultiRun_%d/reads_N1m_L100_Err0.01/snpSum3_SnpLoc0Based_CntFormat1/sensitivity_opt3_SnpReadCov_MultSNP_rsem_ge1.txt'%(rf) 
        run_files[rf] = '/data1/shunfu1/SNPCalling//SimSNPs_MultiRun_%d/reads_N1m_L100_Err0.01/snpSum3_SnpLoc0Based_CntFormat1/sensitivity_opt4.txt'%(rf) 
        #run_cmd('cp %s tmp/sensitivity_opt3_SnpReadCov_MultMapSNP_%d.txt'%(run_files[rf], rf))

    #parse structure
    with open(run_files[run_stt], 'r') as f:

        f.readline()
        f.readline()
        groupInfo = f.readline()
        groupDescriptions = [gd for gd in groupInfo.strip().split('\t')[1:] if gd != ''] #['[%s'%str(v.strip()) for v in groupInfo.split('[') if v[0] != '#' and v != '']
        callerDescriptions = []
        for line in f:
            if line[0]=='#' or line.strip()=='': continue
            callerDescriptions.append(line.split()[0])
        #pdb.set_trace()

    #prepare data
    #pdb.set_trace()
    data = {} #key: (caller_id,group_id) val: list of entry values (per run)
    for g in range(len(groupDescriptions)):
        for c in range(len(callerDescriptions)):
            data[(c,g)]=[] #list of N_runs

    for r in xrange(run_stt, run_stp+1):
        run_file = run_files[r]
        with open(run_file, 'r') as f:
            #skip header lines
            for _ in range(3):
                f.readline()

            c_id = 0
            for line in f:
                if line[0]=='#' or line.strip()=='': continue
                vals = [float(v) for v in line.split()[1:] if v != []]
                for g_id in range(len(vals)):
                    data[(c_id, g_id)].append(vals[g_id])
                c_id += 1
    #pdb.set_trace()

    if '--normalize' in args:
        #pdb.set_trace()
        ref_c_id = int(args[args.index('--normalize')+1])
        for g in range(len(groupDescriptions)):
            #for non ref_c_id entries
            for c in range(len(callerDescriptions)):
                if c==ref_c_id: continue

                for r in range(len(data[(c,g)])):
                    if data[(ref_c_id,g)][r]!=0:
                        data[(c,g)][r] = float(data[(c,g)][r])/data[(ref_c_id,g)][r]
                    else:
                        pass #pdb.set_trace()

            #for ref_c_id entries
            for r in range(len(data[(ref_c_id, g)])):
                if data[(ref_c_id, g)][r]!=0:
                    data[(ref_c_id, g)][r] = 1.0

    #output stat
    #pdb.set_trace()
    outDir = args[args.index('-O')+1]
    outFn = args[args.index('--fn')+1]
    outf = '%s/%s'%(outDir, outFn)
    outf_stdev = outf[:-4]+'_stdev.txt'
    with open(outf,'w') as f, open(outf_stdev, 'w') as f_stdev:
        st = '%s'%groupInfo
        f.write('%s'%st); f_stdev.write('%s'%st);

        for c_id in range(len(callerDescriptions)):
            st1 = '%10s\t'%callerDescriptions[c_id]
            st2 = '%10s\t'%callerDescriptions[c_id]

            for g_id in range(len(groupDescriptions)):
                vals = data[(c_id, g_id)]
                #pdb.set_trace()
                st1 += '%5.2f\t'%(np.mean(vals))
                st2 += '%5.2f\t'%(np.std(vals))

            st1 += '\n'
            st2 += '\n'

            f.write(st1)
            f_stdev.write(st2)

    print('%s written'%outf)
    print('%s written'%outf_stdev)

    return

'''
read a set of roc files (from batch of runs), collapse mean and stdev into output: outDir/outFn.txt & outFn_stdev.txt

roc file format:
line-0: caller0     caller1     etc
line-1(e.g. md): v0     v1      etc
line-2(e.g. fp): w0     w1      etc
etc lines (e.g. tot)

python evaluate_performance_othertools.py --calc_roc_avg_stdev -O outDir --fn outFn (e.g. *.txt)

'''
def calc_roc_avg_stdev(args):

    run_stt = 5
    run_stp = 9

    run_files = []
    for r in xrange(run_stt, run_stp+1):
        #f = '/data1/shunfu1/SNPCalling//SimSNPs_MultiRun_%d//reads_N10m_L100_Err0.01/snpSum3_SnpLoc0Based_CntFormat1//snpSum3_md_fp.txt'%(r)
        #f = '/data1/shunfu1/SNPCalling//SimSNPs_MultiRun_%d//reads_N100K_L100_Err0.00/snpSum3_SnpLoc0Based_CntFormat1//snpSum3_md_fp.txt'%(r)
        #f = '/data1/shunfu1/SNPCalling//SimSNPs_MultiRun_%d//reads_N100K_L100_Err0.00/snpSum3_varscan//snpSum3_md_fp.txt'%(r)
        #f = '/data1/shunfu1/SNPCalling//SimSNPs_MultiRun_%d//reads_N100K_L100_Err0.00/snpSum3_varscan_filt//snpSum3_md_fp.txt'%(r)
        #f = '/data1/shunfu1/SNPCalling//SimSNPs_MultiRun_%d//reads_N100K_L100_Err0.00/snpSum3_samtools/snpSum3_md_fp.txt'%(r)

        #f = '/data1/shunfu1/SNPCalling//SimSNPs_MultiRun_%d//reads_N1m_L100_Err0.01/snpSum3_varscan//snpSum3_md_fp.txt'%(r)
        #f = '/data1/shunfu1/SNPCalling//SimSNPs_MultiRun_%d//reads_N1m_L100_Err0.01/snpSum3_varscan_filt//snpSum3_md_fp.txt'%(r)
        f = '/data1/shunfu1/SNPCalling//SimSNPs_MultiRun_%d//reads_N1m_L100_Err0.01/snpSum3_samtools/snpSum3_md_fp.txt'%(r)
        run_files.append(f)
    
    roc = {} #key - line id (excluding 1st line), val - description of line
    roc[0]='md'
    roc[1]='fp'

    #pdb.set_trace()

    callers = []
    with open(run_files[0], 'r') as f:
        callers = f.readline().split()
    print('callers: %s'%str(callers))

    data = {} #key: (caller id, line id (excluding 1st line)), val: list of caller's val in the line across multi runs
    for r in xrange(run_stt, run_stp+1):
        #pdb.set_trace()
        with open(run_files[r-run_stt], 'r') as f:
            f.readline()
            line_id = 0
            for line in f:
                if line[0]=='#' or line.strip()=='': continue
                if line_id in roc.keys():
                    line_vals = [int(v) for v in line.split()]
                    for c_id in range(len(callers)):
                        key = (c_id, line_id)
                        if key in data:
                            data[key].append(line_vals[c_id])
                        else:
                            data[key]=[line_vals[c_id]]
                    line_id += 1
                else:
                    line_id += 1
    #pdb.set_trace()

    outDir = args[args.index('-O')+1]
    outFn = args[args.index('--fn')+1]
    out_avg = '%s/%s'%(outDir, outFn)
    out_stdev = '%s/%s'%(outDir, outFn[:-4]+'_stdev.txt')
    with open(out_avg, 'w') as oa, open(out_stdev, 'w') as ost:
        header = '\t'.join(callers)
        oa.write(header+'\n')
        ost.write(header+'\n')
        for line_id in roc.keys():
            st = ''
            st1 = ''
            for c_id in range(len(callers)):
                st += '%5.2f\t'%np.mean(data[(c_id, line_id)])
                st1 += '%5.2f\t'%np.std(data[(c_id, line_id)])
            st += '\n'
            st1 += '\n'
            oa.write(st)
            ost.write(st1)
    print('%s written'%out_avg)
    print('%s written'%out_stdev)
    #pdb.set_trace()

    return

def parse_sens(sensFile, num_header_lines=3):

    with open(sensFile, 'r') as f:

        for i in range(num_header_lines-1):
            f.readline()
        groupInfo = f.readline()
        groupDescriptions = [gd for gd in groupInfo.strip().split('\t')[1:] if gd != ''] #['[%s'%str(v.strip()) for v in groupInfo.split('[') if v[0] != '#' and v != '']
        callerDescriptions = []
        for line in f:
            if line[0]=='#' or line.strip()=='': continue
            callerDescriptions.append(line.split()[0])
        #pdb.set_trace()

    return [callerDescriptions, groupDescriptions]

def read_sens_data(sensFile, callerDescriptions, groupDescriptions, num_header_lines=3):

    if sensFile == '':
        return None

    data = {} #key: (caller_id,group_id) val: avg val
    for g in range(len(groupDescriptions)):
        for c in range(len(callerDescriptions)):
            data[(c,g)]=0 #list of N_runs

    with open(sensFile, 'r') as f:
        #skip header lines
        for _ in range(num_header_lines):
            f.readline()

        c_id = 0
        for line in f:
            if line[0]=='#' or line.strip()=='': continue
            vals = [float(v) for v in line.split()[1:] if v != []]
            for g_id in range(len(vals)):
                data[(c_id, g_id)] = vals[g_id]
            c_id += 1

    return data

def barplot_with_errorbar(groups, #e.g. ['100K', '1M']
                          bars,   #e.g. ['abSNP (a=0)', 'abSNP (a=1)', 'GATK']
                          bars_color, #e.g. ['red', 'gold', 'blue']
                          bar_vals, #e.g. #per element (caller) is a list of vals (mean) of a bar/caller versus groups
                          bar_stds,
                          x_lab, #e.g. 'Number of reads'
                          y_lab, #e.g. 'Number of total errors'
                          tit,#e.g. 'Sim 1K SNPs per allele, Read error rate 0.00'
                          bars2=[],
                          bars_color2=[],
                          bar_vals2=[],
                          bar_stds2=[],
                          ):

    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()

    #pdb.set_trace()

    #groups = ['100K', '1M']
    #bars = ['abSNP (a=0)', 'abSNP (a=1)', 'GATK']
    #bars_color = ['red', 'gold', 'blue']
    
    n_groups = len(groups)

    #bottom, mis-detection
    #bar_vals = [] #per element is a list of vals of a bar versus groups
    #bar_stds = []

    #bars[0]
    #bar_vals.append([639.1, 197.375])
    #bar_stds.append([356.86, 137.01])

    #bars[1]
    #bar_vals.append([617.8, 156])
    #bar_stds.append([360.19, 149.09])

    #bars[2]
    #bar_vals.append([729.2, 268.75])
    #bar_stds.append([343.79, 147.04])

    fig, ax = plt.subplots()

    index = np.arange(n_groups)+0.05
    bar_width = 0.1
    index = bar_width/2+index

    for i in range(len(bars)):
        #pdb.set_trace()
        plt.bar(left=index+bar_width*i,
                height=bar_vals[i], 
                width=bar_width,
                color=bars_color[i], #yerr=bar_stds[i],
                label=bars[i],
                error_kw={'ecolor':'black'})

        if bar_vals2 != []:
            #no label
            plt.bar(left=index+bar_width*i,
                height=bar_vals2[i], 
                width=bar_width,
                color=bars_color2[i], #yerr=bar_stds2[i],
                error_kw={'ecolor':'black'},
                bottom=bar_vals[i])

    plt.xlabel(x_lab)
    plt.ylabel(y_lab)
    plt.title(tit)
    plt.xticks(index+len(bars)/2.0*bar_width, groups)
    plt.legend(loc='upper left')
    #plt.legend(loc='upper right')

    ax.set_ylim([0.00,1.1])
    #ax.set_yscale('log')

    #plt.tight_layout()
    plt.show()
    pdb.set_trace()

    return

'''
at local pc
python evaluate_performance_othertools.py --plot_sens                                     

                                     -n num_header_lines

                                     [--skip_groups gid0,gid1,...]
                                     [--groupLabels gLa,gLb,...]
                                     [--skip_callers cid0,cid1,...]

                                     -a0 sensFile [-d0 sensFile_stdev]                                                                          
                                     [--callerLabels0 cLa,cLb,...] #a,b... should correspond to original callers minus skip callers
                                     --callerColors0 cla,clb,...

                                     [-a1 sensFile [-d1 sensFile_stdev]]
                                     [--callerLabels1 cLa,cLb,...] #won't show legend if starts with '_'
                                     [--callerColors1 cla,clb,...] #need to appear if -a1 is there

                                     [-a2 -d2 --callerLabels2 --callerColors2 etc]

                                     --x_lab x_lab
                                     --y_lab y_lab [--y_lim low,high]
                                     [--legend_loc loc]

                                     --title tit

-a<i> and -d<i> i>=1 represent stacked bars

sens file format:
line-0: #snpSum3File
line-1: #groupValCalcOption
line-2: #caller             [g0,g1)     [g1,g2) ... [gN-1, gN)
line-3: caller0 (e.g. Tot)  cnts        cnts        cnts
line-4: caller1             cnts        cnts        cnts
[line-5 etc]

fig will be saved to tmp/figure.png

legend locations
- default: outside fig
- non-default: put in 'upper left', 'upper right' etc of the figure
- see http://matplotlib.org/users/legend_guide.html#plotting-guide-legend
'''
def plot_sens(args):

    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()

    num_header_lines = int(args[args.index('-n')+1])

    x_lab = args[args.index('--x_lab')+1]
    y_lab = args[args.index('--y_lab')+1]
    if '--y_lim' in args:
        y_lim = args[args.index('--y_lim')+1]
        y_lim = [float(v) for v in y_lim.split(',') if v != '']
    else:
        y_lim = [0.0,1.0]
    tit = args[args.index('--title')+1]
    tit = tit.decode('string_escape')

    #pdb.set_trace()

    if '--legend_loc' in args:
        loc = args[args.index('--legend_loc')+1]
    else:
        loc = ''

    if '--skip_groups' in args:
        skip_groups = args[args.index('--skip_groups')+1]
        skip_group_list = [int(v) for v in skip_groups.split(',') if v != '']
    else:
        skip_group_list = []

    groupLabels = [] #list of two lists (one for bottom callers, one for top stacked callers)
    if '--groupLabels' in args:
        groupLabels = args[args.index('--groupLabels')+1]
        groupLabels = [gl for gl in groupLabels.split(',') if gl != '']
    else:
        groupLabels = []

    if '--skip_callers' in args:
        skip_callers = args[args.index('--skip_callers')+1]
        skip_caller_list = [int(v) for v in skip_callers.split(',') if v != '']
    else:
        skip_caller_list = []
    
    f_idx = 0
    prev_callers_vals = []
    while '-a%d'%f_idx in args:

        sensFile = args[args.index('-a%d'%f_idx)+1]
        if '-d%d'%f_idx in args:
            sensFile_stdev = args[args.index('-d%d'%f_idx)+1]
        else:
            sensFile_stdev = ''
        
        #parse structure
        [callerDescriptions, groupDescriptions] = parse_sens(sensFile, num_header_lines=num_header_lines)
        #key - (c_id, g_id) val - val
        data = read_sens_data(sensFile, callerDescriptions, groupDescriptions, num_header_lines=num_header_lines)
        data_stdev = read_sens_data(sensFile_stdev, callerDescriptions, groupDescriptions, num_header_lines=num_header_lines)

        if groupLabels == [] and f_idx==0: #do this only once
            groupLabels = [groupDescriptions[g] for g in list(np.arange(0,len(groupDescriptions))) if g not in skip_group_list]

        if '--callerLabels%d'%f_idx in args:
            callerLabels = args[args.index('--callerLabels%d'%f_idx)+1]
            callerLabels = [gl for gl in callerLabels.split(',') if gl != '']
        else:
            callerLabels = [callerDescriptions[c] for c in list(np.arange(0,len(callerDescriptions))) if c not in skip_caller_list]
        
        callerColors = args[args.index('--callerColors%d'%f_idx)+1]
        callerColors = [gl for gl in callerColors.split(',') if gl != '']
        
        callers_vals = [] #row wrt caller col wrt grp
        callers_stds = []
        for c_id in range(len(callerDescriptions)):
            if c_id in skip_caller_list: continue

            per_caller_vals = []
            per_caller_stds = []
        
            for g_id in range(len(groupDescriptions)):
                if g_id in skip_group_list: continue

                per_caller_vals.append(data[(c_id, g_id)])
                if data_stdev is not None:
                    per_caller_stds.append(data_stdev[(c_id, g_id)])

            callers_vals.append(per_caller_vals)
            if data_stdev is not None:
                callers_stds.append(per_caller_stds)

        ### plot
        bar_width = 0.1
        index = float(bar_width*2)/3+np.arange(len(groupLabels))

        for i in range(len(callerLabels)):
            if data_stdev is None:
                yerr = None
            else:
                yerr = callers_stds[i]
            if f_idx == 0:
                bottom = None
            else:
                bottom = prev_callers_vals[i]

            plt.bar(left=index+bar_width*i,
                    height=callers_vals[i], 
                    width=bar_width,
                    color=callerColors[i],
                    edgecolor=callerColors[i],
                    yerr=yerr,
                    bottom=bottom,
                    label=callerLabels[i],
                    error_kw={'ecolor':'black'})

        ### next iteration
        f_idx += 1
        prev_callers_vals = callers_vals

    ####
    plt.xlabel(x_lab)
    plt.ylabel(y_lab)
    plt.title(tit)

    plt.xticks(index+len(groupLabels)/2.0*bar_width, groupLabels)
    
    if loc=='': #default
        #lgd = plt.legend(loc='upper left', bbox_to_anchor=(1, 1), ncol=1)
        lgd = plt.legend(loc='upper left', bbox_to_anchor=(0, -0.3, 1, 0.2), ncol=len(callerLabels), mode='expand')
    else:
        lgd = plt.legend(loc=loc)

    ax.set_ylim(y_lim)
    ax.yaxis.grid(True)

    fig.savefig('tmp/figure', bbox_extra_artists=(lgd,), bbox_inches='tight') #to put legend outside fig
    print('fig saved to tmp/figure')
    #pdb.set_trace()

    return

'''
at local pc
python evaluate_performance_othertools.py --plot_roc
                                     -x x_line_id(2nd line has id 0) --x_lab x_lab -y y_line_id --y_lab y_lab
                                     --title tit

                                     -i rocFile [-s rocFile_stdev]
                                     -c1 c1_1,c1_2,... -c1_lab c1_lab -c1_color c1_color [-c1_fmt c1_fmt]
                                     -c2 c2_1,c2_2,... -c2_lab c2_lab -c2_color c2_color [-c2_fmt c2_fmt]
                                     [-c3 -c3_lab -c3_color etc]

                                    [--dummy_tot_err_line]

roc file format:
line-0: caller0     caller1     etc
line-1(e.g. md): v0     v1      etc
line-2(e.g. fp): w0     w1      etc
etc lines (e.g. tot)

Options:

-c<i> c<i>_1,c<i>_2... means choose columns c<i>_1,c<i>_2 etc as curve c<i>

-for multiple roc files, we may merge contents into one and select related columns to plot curves

--dummy_tot_err_line: add a -45-degree line cross each curve points to show total error

Misc:

named colors:
http://matplotlib.org/examples/color/named_colors.html

marker reference:
http://matplotlib.org/examples/lines_bars_and_markers/marker_reference.html
'''
def plot_roc(args):    

    rocFile = args[args.index('-i')+1]

    with open(rocFile, 'r') as rocF:
        header = rocF.readline()
        callers = header.split()

        line_vals = [] #each entry is a list (ith line), which is a list of vals vs callers
        for line in rocF:
            line_vals.append([float(v) for v in line.split()])
    #pdb.set_trace()

    x_line_id = int(args[args.index('-x')+1])
    x_lab = args[args.index('--x_lab')+1]
    y_line_id = int(args[args.index('-y')+1])
    y_lab = args[args.index('--y_lab')+1]
    tit = args[args.index('--title')+1]
    tit = tit.decode('string_escape')
    #pdb.set_trace()

    if '--dummy_tot_err_line' in args:
        dummy_line = 1
    else:
        dummy_line = 0

    curve_idx = 1

    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()

    while '-c%d'%curve_idx in args:

        c_ids = args[args.index('-c%d'%curve_idx)+1] #curve points (e.g. caller0 and caller1 belong to the same curve)
        c_ids = [int(c) for c in c_ids.split(',') if c!='']

        x_vals = [line_vals[x_line_id][c] for c in c_ids]
        y_vals = [line_vals[y_line_id][c] for c in c_ids]

        c_lab = args[args.index('-c%d_lab'%curve_idx)+1]
        c_color = args[args.index('-c%d_color'%curve_idx)+1]
        
        if '-c%d_fmt'%curve_idx in args:
            c_fmt = args[args.index('-c%d_fmt'%curve_idx)+1]
        else:
            c_fmt = '-o'


        ax.errorbar(x_vals,
                    y_vals,
                    fmt=c_fmt,
                    label=c_lab,
                    color=c_color)

        '''
        for i in range(len(x_vals)):
            x = int(x_vals[i])
            y = int(y_vals[i])
            ax.annotate('(%d, %d)'%(x,y), xy=(x+1,y+1), textcoords='data')
        
        ss = 65
        tt = 189
        cl = 'red'
        ax.annotate('m = %d - f'%tt, xy=(ss,tt-ss), xytext=(ss+15,tt-ss+10), arrowprops=dict(facecolor=cl, shrink=0.05)) #, textcoords='data')

        ss = 55
        tt = 127
        cl = 'blue'
        ax.annotate('m = %d - f'%tt, xy=(ss,tt-ss), xytext=(ss+15,tt-ss+10), arrowprops=dict(facecolor=cl, shrink=0.05)) #, textcoords='data')

        ss = 75
        tt = 102
        cl = 'blue'
        ax.annotate('m = %d - f'%tt, xy=(ss,tt-ss), xytext=(ss+20,tt-ss+15), arrowprops=dict(facecolor=cl, shrink=0.05)) #, textcoords='data')
        '''

        if dummy_line==1:
            for i in range(len(x_vals)): #each point has a dummy line
                tot = x_vals[i]+y_vals[i]
                x_dummy_vals = [0, tot]
                y_dummy_vals = [tot, 0]
                ax.errorbar(x_dummy_vals,
                            y_dummy_vals,
                            fmt='--',
                            color=c_color)       

        curve_idx += 1
    #pdb.set_trace()

    plt.xlabel(x_lab)
    plt.ylabel(y_lab)
    plt.title(tit)
    lgd = plt.legend(loc='upper right')
    ax.grid(True)

    fig.savefig('tmp/figure')
    #fig.savefig('tmp/figure', bbox_extra_artists=(lgd,), bbox_inches='tight') #to put legend outside fig
    print('fig saved to tmp/figure')

    return

'''
usage:

#selective count alt generation --> gen snpSum3 --> detailed sensitivity for different runs of certain case (e.g. fixed read num, len and err rate)
#generate: sel count, count alt, snpSum3 and roc & sensitivity files

python evaluate_performance_othertools.py --batch

#calc avg/stdev of a set of sens files. and output outDir/outFn.txt & outFn_stdev.txt
#
#sens file format:
#line-0: #snpSum3File
#line-1: #groupValCalcOption
#line-2: #caller   [g0,g1)  [g1,g2) ... [gN-1, gN)
#line-3: tot       cnts    cnts        cnts
#line-4: caller0   cnts    cnts        cnts
#[line-5 etc]

python evaluate_performance_othertools.py --calc_sens_avg_stdev -O outDir --fn outFn (e.g. *.txt) [--normalize ref_c_id]

#at local pc
python evaluate_performance_othertools.py --plot_sens -i sensFile [-i2 sensFile_stdev]
'''
if __name__ == "__main__":

    args = sys.argv

    if '--batch' in args:
        evaluate_performance_batch(args)
    elif '--calc_sens_avg_stdev' in args:
        calc_sens_avg_stdev(args)
    elif '--calc_roc_avg_stdev' in args:
        calc_roc_avg_stdev(args)
    elif '--plot_sens' in args:
        plot_sens(args)
    elif '--plot_roc' in args:
        plot_roc(args)
