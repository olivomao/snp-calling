import sys, os, pdb
import numpy as np

from old_code.util import run_cmd, run_cmds

def test_sensitivity_analysis_batch(args):

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
    #'''
    groupValCalcOption = 1
    groupValsStr = '0,10,100,1e3,1e4,1e6,1e8,1e10'
    #'''

    #opt-2 based on snpSum3 gatk/rsem altMap
    #    esp, col-7 for count file in format 1 (gatk) and col-8 for count file in format 1 (rsem), esp:
    #    -1: no snp reads
    #    0:  part or all of snp reads uniq mapped (part, not all of snp reads alt mapped)
    #    1+ (float val): all snp reads alt mapped, # of avg alt mapping (excluding self)
    '''
    groupValCalcOption = 2
    groupValsStr = '-1,0,1,2,4,6,8,10,20'
    '''
    
    old_code_dir = 'old_code/'
    srcDir = '/data1/shunfu1/SNPCalling/'

    T = 1 #caller threshold
    rocT=[0.00, 1.00]

    chrom = 'Chr15'
    refGenome = '/data1/shunfu1/SNPCalling/data_large_0_idealCov/%s.fa'%chrom

    num_p = 20 #num parallel

    cnt_format = 1 #0-e.g. C,I,lambda_bj,lambda_sum  1-e.g. bj,num_alt_mapping

    snpSum3_folder_name = 'snpSum3_SnpLoc0Based_CntFormat1'#'snpSum_3_debug_gen_count_sel' #'snpSum3_2' para applied
    #sensitivity_res_name = 'sensitivity_opt%d.txt'%groupValCalcOption
    sensitivity_res_name = 'sensitivity_opt1.txt' #'sensitivity_opt2_rsemAltMap_le0.txt' #'sensitivity_opt2.txt' #'sensitivity_opt2_rsemAltMap_ge1.txt'

    #### config - steps to run
    do_bam2sam = 0
    do_gen_count_selectively = 0
    do_gen_snpSum3 = 0
    do_detailed_sens_analysis = 1


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
        run_files[rn]['s_ab_0'] = '%s/%s/snp_res/T%d/caller_output_snp_found_%s.genome.sorted_n_T%d_para_rocT_%.2f.txt'%(runDir, readsLabel, T, chrom, T, rocT[0])
        run_files[rn]['s_ab_1'] = '%s/%s/snp_res/T%d/caller_output_snp_found_%s.genome.sorted_n_T%d_para_rocT_%.2f.txt'%(runDir, readsLabel, T, chrom, T, rocT[1])
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
                                       '-s2 %s '%run_files[rn]['s_ab_1']+\
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
                                       '-s2 %s '%run_files[rn]['s_ab_1']+\
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
        cmd = 'python evaluator.py --gen_snpSum3 '+\
                                   '-m %s '%run_files[rn]['s_m']+\
                                   '-p %s '%run_files[rn]['s_p']+\
                                   '--L1 %s --F1 %s '%('abSNP_a0', run_files[rn]['s_ab_0'])+\
                                   '--L2 %s --F2 %s '%('abSNP_a1', run_files[rn]['s_ab_1'])+\
                                   '--L3 %s --F3 %s '%('GATK', run_files[rn]['s_gatk'])+\
                                   '--cG1 %s '%run_files[rn]['cnt_gatk']+\
                                   '--cR1 %s '%run_files[rn]['cnt_rsem']+\
                                   '-O %s '%run_files[rn]['snpSum3Dir']+\
                                   '--rm %s '%run_files[rn]['m_rd_bed']+\
                                   '--rp %s '%run_files[rn]['p_rd_bed']+\
                                   '--removeIntermediateFiles'
        if num_p == 1:
            if do_gen_snpSum3==1: run_cmd(cmd) #pass; #run_cmd(cmd)
        else:
            cmds.append(cmd)

    if num_p > 1:
        if do_gen_snpSum3==1: run_cmds(cmds, num_p) #pass #run_cmds(cmds, num_p)

    #### 3: detailed sensitivity
    cmds = []
    for rn in xrange(run_stt, run_stp+1):
        #detailed sensitivity analysis
        cmd = 'python evaluator.py --sensitivity_analysis '+\
                                '--snpSum3 %s '%(run_files[rn]['snpSum3File'])+\
                                '--groupValCalcOption %d '%groupValCalcOption+\
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

python test_sensitivity_analysis2.py --calc_sens_avg_stdev -O outDir --fn outFn (e.g. *.txt)
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

        run_files[rf] = '/data1/shunfu1/SNPCalling//SimSNPs_MultiRun_%d/reads_N1m_L100_Err0.01/snpSum3_SnpLoc0Based_CntFormat1/sensitivity_opt1.txt'%(rf)
        

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
                st1 += '%5.2f\t'%(np.mean(vals))
                st2 += '%5.2f\t'%(np.std(vals))

            st1 += '\n'
            st2 += '\n'

            f.write(st1)
            f_stdev.write(st2)

    print('%s written'%outf)
    print('%s written'%outf_stdev)

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

    index = np.arange(n_groups)
    bar_width = 0.1
    index = bar_width/2+index

    for i in range(len(bars)):
        #pdb.set_trace()
        plt.bar(left=index+bar_width*i,
                height=bar_vals[i], 
                width=bar_width,
                color=bars_color[i],
                yerr=bar_stds[i],
                label=bars[i],
                error_kw={'ecolor':'black'})

        if bar_vals2 != []:
            #no label
            plt.bar(left=index+bar_width*i,
                height=bar_vals2[i], 
                width=bar_width,
                color=bars_color2[i],
                yerr=bar_stds2[i],
                error_kw={'ecolor':'black'},
                bottom=bar_vals[i])

    plt.xlabel(x_lab)
    plt.ylabel(y_lab)
    plt.title(tit)
    plt.xticks(index+len(bars)/2.0*bar_width, groups)
    plt.legend(loc='upper right')

    ax.set_ylim([-5,ax.get_ylim()[1]])
    #ax.set_yscale('log')

    #plt.tight_layout()
    plt.show()
    pdb.set_trace()

    return

'''
at local pc
python test_sensitivity_analysis2.py --plot -i sensFile [-i2 sensFile_stdev] -n num_header_lines
                                            [-s sensFile2 [-s2 sensFile2_stdev]]

-s and -s2 represent stacked bars
'''
def plot_sens(args):

    x_lab = 'Num of Alternative Mapping'
    y_lab = 'SNPs correctly detected'
    tit = 'Sim 1K SNPs per allele, 1M 100-bp Reads with error rate 0.01'

    colors = ['red', 'green', 'blue', 'gold', 'black']
    colors2 = ['orangered', 'lime', 'lightskyblue', 'khaki', 'gray']

    #plot
    print('skip_group_list:');pdb.set_trace()
    skip_group_list = [5,6]
    skip_caller_list = []

    sensFile = args[args.index('-i')+1]
    if '-i2' in args:
        sensFile_stdev = args[args.index('-i2')+1]
    else:
        sensFile_stdev = ''

    if '-s' in args:
        sensFile2 = args[args.index('-s')+1]
    else:
        sensFile2 = ''

    if '-s2' in args:
        sensFile2_stdev = args[args.index('-s2')+1]
    else:
        sensFile2_stdev = ''

    num_header_lines = int(args[args.index('-n')+1])

    #parse structure
    [callerDescriptions, groupDescriptions] = parse_sens(sensFile, num_header_lines=num_header_lines)
    #pdb.set_trace()

    #key - (c_id, g_id) val - val
    data = read_sens_data(sensFile, callerDescriptions, groupDescriptions, num_header_lines=num_header_lines)
    #pdb.set_trace()
    
    if sensFile_stdev != '' and os.path.exists(sensFile_stdev)==True:
        data_stdev = read_sens_data(sensFile_stdev, callerDescriptions, groupDescriptions, num_header_lines=num_header_lines)
        #pdb.set_trace()
    else:
        data_stdev = None

    if sensFile2 != '':
        data2 = read_sens_data(sensFile2, callerDescriptions, groupDescriptions, num_header_lines=num_header_lines)
    else:
        data2 = None

    if sensFile2_stdev != '':
        data2_stdev = read_sens_data(sensFile2_stdev, callerDescriptions, groupDescriptions, num_header_lines=num_header_lines)
    else:
        data2_stdev = None    

    groups = [groupDescriptions[g] for g in list(np.arange(0,len(groupDescriptions))) if g not in skip_group_list]
    callers = [callerDescriptions[c] for c in list(np.arange(0,len(callerDescriptions))) if c not in skip_caller_list]
    if data2 is not None:
        callers2 = [callerDescriptions[c]+'(2)' for c in list(np.arange(0,len(callerDescriptions))) if c not in skip_caller_list]
    else:
        callers2 = []

    
    bar_vals = []
    bar_stds = []

    bar_vals2 = [] #to be stacked
    bar_stds2 = []

    for c_id in range(len(callerDescriptions)):
        if c_id in skip_caller_list: continue

        per_bar_vals = []
        per_bar_stds = []

        per_bar_vals2 = []
        per_bar_stds2 = []
        
        for g_id in range(len(groupDescriptions)):
            if g_id in skip_group_list: continue

            per_bar_vals.append(data[(c_id, g_id)])
            per_bar_stds.append(data_stdev[(c_id, g_id)])

            if data2 is not None:
                if c_id==0:
                    per_bar_vals2.append(0)
                else:
                    per_bar_vals2.append(data2[(c_id, g_id)])
            if data2_stdev is not None:
                if c_id==0:
                    per_bar_stds2.append(0)
                else:
                    per_bar_stds2.append(data2_stdev[(c_id, g_id)])

        bar_vals.append(per_bar_vals)
        bar_stds.append(per_bar_stds)

        if data2 is not None:
            bar_vals2.append(per_bar_vals2)
        if data2_stdev is not None:
            bar_stds2.append(per_bar_stds2)
    pdb.set_trace()

    barplot_with_errorbar(groups=groups, #e.g. ['100K', '1M']
                          bars=callers,   #e.g. ['abSNP (a=0)', 'abSNP (a=1)', 'GATK']
                          bars_color=colors[:len(callerDescriptions)], #e.g. ['red', 'gold', 'blue']
                          bar_vals=bar_vals, #e.g. #per element (caller) is a list of vals (mean) of a bar/caller versus groups
                          bar_stds=bar_stds,
                          x_lab=x_lab, #e.g. 'Number of reads'
                          y_lab=y_lab, #e.g. 'Number of total errors'
                          tit=tit, #e.g. 'Sim 1K SNPs per allele, Read error rate 0.00'
                          bars2 = callers2,
                          bars_color2=colors2[:len(callerDescriptions)],
                          bar_vals2=bar_vals2,
                          bar_stds2=bar_stds2
                          )
    return

'''
usage:

#selective count alt generation --> gen snpSum3 --> detailed sensitivity for different runs of certain case (e.g. fixed read num, len and err rate)
#generate: sel count, count alt, snpSum3 and sensitivity files

python test_sensitivity_analysis2.py --batch

#calc avg/stdev of a set of sens files. and output outDir/outFn.txt & outFn_stdev.txt
#
#sens file format:
#line-0: #snpSum3File
#line-1: #groupValCalcOption
#line-2: #caller   [g0,g1)  [g1,g2) ... [gN-1, gN)
#line-3: tot       cnts    cnts        cnts
#line-4: caller0   cnts    cnts        cnts
#[line-5 etc]

python test_sensitivity_analysis2.py --calc_sens_avg_stdev -O outDir -fn outFn

#at local pc
python test_sensitivity_analysis2.py --plot -i sensFile [-i2 sensFile_stdev]
'''
if __name__ == "__main__":

    args = sys.argv

    if '--batch' in args:
        test_sensitivity_analysis_batch(args)
    elif '--calc_sens_avg_stdev' in args:
        calc_sens_avg_stdev(args)
    elif '--plot' in args:
        plot_sens(args)