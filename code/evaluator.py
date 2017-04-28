import sys, pdb
from intervaltree import Interval, IntervalTree 
import numpy as np

from sim_data_generator import snp_read_cov1

from old_code.util import run_cmd
from old_code.snp_analysis import SnpNodeDic
from old_code.snp_oper import load_snps
from old_code.count_read_lambda import check_lambda_multimapping, check_count_fmt1_lam_altMap


'''
python evaluator.py    --gen_snpSum3
                       -m true_m_snp
                       -p true_p_snp
                       --L1 snp_res_label_1 --F1 snp_res_file_1 [--L2 L2 --F2 F2 etc] (s.t. we can have flexible number of snp files)
                       (--caG count_alt_gatk --caR count_alt_rsem) or (--cG1 count_gatk_fmt1 --cR1 count_rsem_fmt1)
                       -O outDir
                       [--rm true_read_m_bed --rp true_read_p_bed]
                       [--removeIntermediateFiles]
                       [--md_fp]
#output:
outDir/snpSum3.txt
       snps_mp.txt (intermediate file; pooled snps)
       snps_rd_cov.txt (intermediate file)
       snpSum3_md_fp.txt (if --md_fp specified)

snpSum3.txt:
line-0: header line describe following cols
line-1 and so on:
col-0: gPos (1-based)
col-1: rB
col-2: tB
col-3: allele (m or p or b)
col-4: # of m reads covering it (-1 if info unavailable)
col-5: # of p reads covering it (-1 if info unavailable)
col-6: estimated abundance
col-7 & 8:
    if --caG and --caR used:
        col-7: # of multi mapped reads (gatk), excluding self (-1 if info unavailable)
        col-8: # of multi mapped reads (rsem), excluding self (-1 if info unavailable)
    if --cG1 and --cR1 used:
        col-7 for count file in format 1 (gatk) and col-8 for count file in format 1 (rsem), esp:
               -1: no snp reads
               0:  part or all of snp reads uniq mapped (part, not all of snp reads alt mapped)
               1+ (float val): all snp reads alt mapped, # of avg alt mapping (excluding self)
col-9: is detected by snp_res_file_1 (e.g. abSNP a=0) (1 or 0 or 0.5/pos right, but different tB)
[col-10: is detected by snp_res_file_2 (e.g. abSNP a=1) (1 or 0 or 0.5/pos right, but different tB)]
[col-11: is detected by snp_res_file_3 (e.g. GATK) (1 or 0 or 0.5/pos right, but different tB)]

snpSum3_md_fp.txt
- Note: SNPs w/o SNP reads are excluded from md calculation
line-0: caller0     caller1     etc
line-1: md0         md1         etc
line-2: fp0         fp1         etc
'''
def gen_snpSum3(args):

    #pdb.set_trace()

    outDir = args[args.index('-O')+1]
    intermediate_files = []

    snpSum3 = {} #key: gPos 1-based val: list of vals for col-1 ~ col-9+ of snpSum3

    true_m_snp = load_snps(args[args.index('-m')+1])
    true_p_snp = load_snps(args[args.index('-p')+1])

    i=0 #num of additional snp res files
    op_L='--L%d'%(i+1) 
    op_F='--F%d'%(i+1) 
    snp_res_files = {} #key - file description val - file path
    while op_L in args:
        f_description = args[args.index(op_L)+1]
        f_path = args[args.index(op_F)+1]
        snp_res_files[f_description] = f_path
        i+=1
        op_L='--L%d'%(i+1) 
        op_F='--F%d'%(i+1)

    snp_res_files = snp_res_files.items()
    snp_res_files = sorted(snp_res_files, key=lambda x:x[0])
    print('snp_res_files: %s'%str(snp_res_files))
    #pdb.set_trace()

    #col-1 ~ col-3
    for gPos, rB_tB in true_m_snp.items():
        snpSum3[gPos] = [rB_tB[0], rB_tB[1], 'm', -1, -1, 0.0, -1, -1]+[0]*i 

    for gPos, rB_tB in true_p_snp.items():
        if gPos in snpSum3:
            snpSum3[gPos][2] = 'b'
        else:
            snpSum3[gPos] = [rB_tB[0], rB_tB[1], 'p', -1, -1, 0.0, -1, -1]+[0]*i

    int_snps_mn = '%s/snps_mn.txt'%(outDir); intermediate_files.append(int_snps_mn)
    with open(int_snps_mn, 'w') as int_snps_mn_f:
        for gPos, vals in snpSum3.items():
            int_snps_mn_f.write('%d\t%s\t-->\t%s\n'%(gPos, vals[0], vals[1]))
    #pdb.set_trace()

    #col-9+
    #pdb.set_trace()
    for i in range(len(snp_res_files)):
        f_description = snp_res_files[i][0]
        f_path = snp_res_files[i][1]

        detected_snps = load_snps(f_path)
        
        for gPos, rB_tB in detected_snps.items():
            if gPos in snpSum3:
                if rB_tB[1]==snpSum3[gPos][1]:
                    snpSum3[gPos][8+i]=1
                else:
                    snpSum3[gPos][8+i]=0.5

    #col-4 ~ col-5
    if '--rm' in args and '--rp' in args:
        m_read_bed = args[args.index('--rm')+1]
        p_read_bed = args[args.index('--rp')+1]
        int_snp_rd_cov = '%s/snps_rd_cov.txt'%outDir; intermediate_files.append(int_snp_rd_cov)
        curr_args = '-s %s -m %s -p %s -o %s'%(int_snps_mn, m_read_bed, p_read_bed, int_snp_rd_cov)
        snp_read_cov1(curr_args.split())
        with open(int_snp_rd_cov, 'r') as int_snp_rd_cov_f:
            for line in int_snp_rd_cov_f:
                if line[0]!='#' and len(line.strip().split())>=3:
                    tokens = line.strip().split()
                    gPos = int(tokens[0])
                    mCov = int(tokens[1])
                    pCov = int(tokens[2])
                    if gPos not in snpSum3:
                        print('unexpected'); pdb.set_trace()
                    snpSum3[gPos][3]=mCov
                    snpSum3[gPos][4]=pCov
    #pdb.set_trace()

    #col-6 ~ col-8
    if '--caG' in args and '--caR' in args:
        if '--caG' in args:
            lam_multimapping_count_alt_gatk = check_lambda_multimapping(args[args.index('--caG')+1])
            for k in snpSum3.keys():
                k2 = k#-1 #0-based
                if k2 in lam_multimapping_count_alt_gatk:
                    snpSum3[k][5] = lam_multimapping_count_alt_gatk[k2][0]
                    snpSum3[k][6] = lam_multimapping_count_alt_gatk[k2][1]
                else:
                    print('unexpected: k2=%d not in count_alt_gatk'%k2)
                    pdb.set_trace()
        else:
            lam_multimapping_count_alt_gatk = None

        if '--caR' in args:
            lam_multimapping_count_alt_rsem = check_lambda_multimapping(args[args.index('--caR')+1])
            for k in snpSum3.keys():
                k2 = k#-1 #0-based
                if k2 in lam_multimapping_count_alt_rsem:
                    snpSum3[k][5] = lam_multimapping_count_alt_rsem[k2][0]
                    snpSum3[k][7] = lam_multimapping_count_alt_rsem[k2][1]
                else:
                    print('unexpected: k2=%d not in count_alt_rsem'%k2)
                    pdb.set_trace()
        else:
            lam_multimapping_count_alt_rsem = None
        #pdb.set_trace()
    elif '--cG1' in args and '--cR1' in args:
        c_cases = []
        c_cases.append(['--cG1', 6]) #col-7
        c_cases.append(['--cR1', 7]) #col-8

        for c_case in c_cases:
            c_op = c_case[0] #count file operator (indicator)
            snpSum3_idx = c_case[1]
            #pdb.set_trace()

            lam_altMap_count = check_count_fmt1_lam_altMap(args[args.index(c_op)+1])
            for k in snpSum3.keys():
                k2 = k
                if k2 in lam_altMap_count:
                    snpSum3[k][5] = lam_altMap_count[k2][0]
                    snpSum3[k][snpSum3_idx] = lam_altMap_count[k2][1]
                else:
                    print('unexpected: k2=%d not in %s'%(k2, c_op))
                    pdb.set_trace()

            #pdb.set_trace()

    else:
        print('need (--caG --caR) or (--cG1 --cR1)'); pdb.set_trace()

    if '--md_fp' in args:
        #pdb.set_trace()
        md = [0]*len(snp_res_files)
        fp = [0]*len(snp_res_files)
        #check md
        for gPos, v_list in snpSum3.items():
            if float(v_list[7-1])==-1.0:
                #pdb.set_trace()
                continue
            for c_id in range(len(snp_res_files)):
                if v_list[9-1+c_id]!=1:
                    #pdb.set_trace()
                    md[c_id]+=1
        #check_fp
        for i in range(len(snp_res_files)):
            f_description = snp_res_files[i][0]
            f_path = snp_res_files[i][1]

            detected_snps = load_snps(f_path)
        
            for gPos, rB_tB in detected_snps.items():
                if gPos not in snpSum3:
                    #pdb.set_trace()
                    fp[i]+=1
                elif rB_tB[1]!=snpSum3[gPos][1]:
                    pdb.set_trace()
                    fp[i]+=1
        #pdb.set_trace()
        out_file_md_fp = '%s/snpSum3_md_fp.txt'%(outDir)
        with open(out_file_md_fp, 'w') as of:
            st = ''
            for i in range(len(snp_res_files)):
                st += snp_res_files[i][0]+'\t'
            st += '\n'
            of.write(st)

            of.write('\t'.join([str(md_val) for md_val in md])+'\n')
            of.write('\t'.join([str(fp_val) for fp_val in fp])+'\n')

            print('%s written'%(out_file_md_fp))

    #output
    out_file = '%s/snpSum3.txt'%outDir
    with open(out_file, 'w') as of:

        #headline
        #st = '#gPos1\trB\ttB\tallele\tnum_m_rds\tnum_p_rds\tab\tmultMap_gatk\tmultMap_rsem\t'
        st = '#gPos1\trB\ttB\tallele\tnum_m_rds\tnum_p_rds\tab\taltMap_gatk\taltMap_rsem\t'
        for i in range(len(snp_res_files)):
            f_description = snp_res_files[i][0]
            st += f_description + '\t'
        st += '\n'
        of.write(st)

        snpSum3 = snpSum3.items()
        snpSum3 = sorted(snpSum3, key=lambda x:x[0])
        for k, v_list in snpSum3:
            st = '%d\t'%k
            st += '\t'.join([str(v) for v in v_list])
            st += '\n'
            of.write(st)
    #pdb.set_trace()

    if '--removeIntermediateFiles' in args:
        for int_f in intermediate_files:
            run_cmd('rm %s'%int_f)
    #pdb.set_trace()

    return snpSum3

#return: snpSum3 object {} (key - gPos 1 based val - list of col1 ~ col9+)
#        caller_descriptions: list of caller names embedded in snpSum3 file
def load_snpSum3(file):
    snpSum3 = {}
    with open(file, 'r') as f:

        header_line = f.readline()
        caller_descriptions = [itm for itm in header_line[1:].split() if itm != ''][9:] #tokens[9:] describe callers
        print('caller descriptions: %s'%str(caller_descriptions))

        for line in f:
            if line[0]=='#' or line.strip()=='': continue
            tokens = line.split()
            gPos = int(tokens[0]) #1-based
            rB = tokens[1] #col-1: rB
            tB = tokens[2] #col-2: tB
            allele = tokens[3] #col-3: allele (m or p or b)
            m_cov = int(tokens[4]) #col-4: # of m reads covering it (-1 if info unavailable)
            p_cov = int(tokens[5]) #col-5: # of p reads covering it (-1 if info unavailable)
            ab = float(tokens[6]) #col-6: estimated abundance
            mult_map_gatk = float(tokens[7]) #col-7: # of multi mapped reads (gatk), excluding self (-1 if info unavailable)
            mult_map_rsem = float(tokens[8]) #col-8: # of multi mapped reads (rsem), excluding self (-1 if info unavailable)
            snpSum3[gPos]=[rB, tB, allele, m_cov, p_cov, ab, mult_map_gatk, mult_map_rsem]

            for c in tokens[9:]:
                snpSum3[gPos].append(int(c))

    return [caller_descriptions, snpSum3]

#v_list: list corresponding to snpSum3 file's data line of col1 ~ col9+
#groupValCalcOption 0~3, 0 (based on m cov/ p cov) 1 (based on ab) 2 (based on gatk/rsem multimapping) 3 (based on ab perc, need ab_perc)
#
#ab_perc is an interval tree with [ab_stt, ab_stp):perc
#we can use it to decide which target perc a SNP falls into
#
#return the val (indicates which value the snp line needs to have in order to be assigned to some group); None returned in case of any exception
def calc_curr_group_val(v_list, groupValCalcOption, ab_perc=None):
    val = None
    if groupValCalcOption==0:
        #pdb.set_trace()
        val = v_list[4-1]+v_list[5-1] # sum of m_cov and p_cov
    elif groupValCalcOption==1:
        #pdb.set_trace()
        val = v_list[6-1]
    elif groupValCalcOption==2:
        #pdb.set_trace()
        val = v_list[7-1]#7-1: gatk 8-1: rsem
        #val = v_list[8-1]
        #val = min(v_list[7-1], v_list[8-1])
        #val = max(v_list[7-1], v_list[8-1])
    elif groupValCalcOption==3:
        if ab_perc is None:
            print('unexpected ab_perc')
            pdb.set_trace()
        else:
            #pdb.set_trace()
            ab = float(v_list[6-1])
            #if v_list[3-1]=='m' or v_list[3-1]=='b': ab = float(v_list[4-1])
            #if v_list[3-1]=='p': ab = float(v_list[5-1])

            perc_query = ab_perc[ab]
            if len(perc_query)!=0:
                val =  sorted(perc_query)[0].data #perc
            else:
                print('cant find per from perc_query, ab=%f'%ab)
                pdb.set_trace()

    return val

#groupValCalcOption 0~2, 0 (based on m cov/ p cov) 1 (based on ab) 2 (based on gatk/rsem multimapping)
#group_vals: a list of float values corresponding to g0,g1 to gN
#snpSum3: snpSum3 object {} (key - gPos 1 based val - list of col1 ~ col9+)
#caller_descriptions: list of caller names embedded in snpSum3 file / ordered
#
#output:
#resFile & resFile.log
#
#resFile:
#line-0: #snpSum3File
#line-1: #groupValCalcOption
#line-2: #caller   [g0,g1)  [g1,g2) ... [gN-1, gN)
#line-3: tot       cnts    cnts        cnts
#line-4: caller0   cnts    cnts        cnts
#[line-5 etc]
#
#resFile.log
#(caller0, 0-th grp): 
#list of snp locations (1-based)
#...
#(caller0, N-1-th grp):
#list of snp locations
#(caller1, 0-th grp):
#list of snp locations
#...
#(caller1, N-1-th grp):
#list of snp locations
#...
def count_by_group(groupValCalcOption, group_vals, snpSum3File, snpSum3, caller_descriptions, resFile, ab_perc=None, perc_ab=None):

    N = len(group_vals)-1
    
    treeGroup = IntervalTree() #data is group index
    for i in range(N):
        g_stt = group_vals[i]
        g_stp = group_vals[i+1]
        treeGroup.add(Interval(g_stt, g_stp, i))

    Counts = {} #key-caller/row val-list of cnts falling into per group
    for c in caller_descriptions:
        Counts[c] = [0]*N
    Counts['Tot'] = [0]*N 

    Counts2 = {} #key-(caller/excluding 'Tot', group_id) val-list of snp locations (1 based) falling into group_id called by caller
    for c in caller_descriptions:
        for i in range(N):
            key = (c, i)
            Counts2[key]=[]

    for gPos, v_list in snpSum3.items():
        val = calc_curr_group_val(v_list, groupValCalcOption, ab_perc=ab_perc) ##
        if val is None:
            print('unexpected group val'); pdb.set_trace()
            continue
        #pdb.set_trace()
        query_res = treeGroup[val]
        if len(query_res)!=0:
            index =  sorted(query_res)[0].data #group id

            Counts['Tot'][index] += 1

            for c_idx in range(len(caller_descriptions)):
                if v_list[8+c_idx]==1:
                    c = caller_descriptions[c_idx]
                    #print('modify here for sub grouping'); #pdb.set_trace()
                    if True: #v_list[8-1]>=1:
                        Counts[c][index] += 1
                        Counts2[(c,index)].append(gPos)
                    #pdb.set_trace()
        else:
            print('find no group for val %s (gPos: %d v_list: %s groupValCalcOption: %d)'%(str(val), gPos, str(v_list), groupValCalcOption))
            pdb.set_trace()

    #print and output
    print('ready to output sensitivity res and log.')
    #pdb.set_trace()
    resFileLog = resFile + '.log'
    with open(resFile, 'w') as rF, open(resFileLog, 'w') as rF_log:

        #rF and print
        st = '#snpSum3File: %s\n'%snpSum3File
        st += '#groupValCalcOption: %d\n'%groupValCalcOption
        st += '%10s\t'%'#caller   '
        #st += '%10s\t'%'#caller   '
        for i in range(N):
            g_stt = group_vals[i]
            g_stp = group_vals[i+1]
            if groupValCalcOption==3:
                #change perc to ab
                g_stt = perc_ab[g_stt]
                if g_stp == 100:
                    g_stp = np.inf 
                else:
                    g_stp = perc_ab[g_stp]
            st += '%10s\t'%('[%.2f, %.2f)'%(g_stt, g_stp))
            #st += '%s\t'%('%.2f'%(g_stt))
        st += '\n'
        print(st)
        rF.write(st) #header lines info

        #Tot
        st = '%10s\t'%('Tot')
        st += '\t'.join(['%10s'%str(v) for v in Counts['Tot']])
        #st = '%10s\t'%('Tot')
        #st += '\t'.join(['%s'%str(v) for v in Counts['Tot']])
        st += '\n'
        print(st)
        rF.write(st)

        for caller in caller_descriptions:
            st = '%10s\t'%caller
            st += '\t'.join(['%10s'%str(v) for v in Counts[caller]])
            #st = '%10s\t'%caller
            #st += '\t'.join(['%s'%str(v) for v in Counts[caller]])
            st += '\n'
            print(st)
            rF.write(st)

        #rF.log
        for caller in caller_descriptions:
            for i in range(N):
                key = (caller, i)
                st = '%s:\n'%str(key)
                st += ','.join([str(v) for v in Counts2[key]])
                st += '\n'
                rF_log.write(st)

    return

#return an interval tree and perc_ab dic (key-perc, val-ab)
#[ab0, ab1): percentile 0
#...
#[ab99,np.inf): percentile 99
def load_ab_perc(ab_perc_file):
    perc_ab = []
    #pdb.set_trace()
    with open(ab_perc_file, 'r') as f:
        for line in f:
            perc = int(line.split()[0])
            ab = float(line.split()[1])
            perc_ab.append([perc, ab])

    #pdb.set_trace()
    res = IntervalTree() #data is group index
    for i in range(len(perc_ab)):
        if i<len(perc_ab)-1:
            ab_stt = perc_ab[i][1]
            ab_stp = perc_ab[i+1][1]
        else:
            ab_stt = perc_ab[i][1]
            ab_stp = np.inf
        res.add(Interval(ab_stt, ab_stp, perc_ab[i][0]))
    #pdb.set_trace()
    res2 = {}
    for perc, ab in perc_ab:
        res2[perc]=ab
    return res, res2


#load snpSum3 file and analyze sensitivity
#usage:
#python evaluator --sensitivity_analysis
#                 --snpSum3 snpSum3File
#                 --groupValCalcOption 0~3     # 0 (based on m cov/ p cov) 1 (based on ab) 2 (based on gatk/rsem multimapping) 3 (based on ab percentile)
#                 --ab_percentile ab_percentile_file
#                 --groupVals g0,g1,g2,...,gN     # [g0,g1) ... [gN-1, gN)
#                 -o resFile
#output:
#resFile & resFile.log
#
#resFile:
#line-0: #snpSum3File
#line-1: #groupValCalcOption
#line-2: #caller   [g0,g1)  [g1,g2) ... [gN-1, gN)
#line-3: tot       cnts    cnts        cnts
#line-4: caller0   cnts    cnts        cnts
#[line-5 etc]
#
#resFile.log
#(caller0, 0-th grp): list of snp locations (1-based)
#...
#(caller0, N-1-th grp): list of snp locations
#(caller1, 0-th grp): list of snp locations
#...
#(caller1, N-1-th grp): list of snp locations
#...
def sensitivity_analysis(args):

    snpSum3File = args[args.index('--snpSum3')+1]
    groupValCalcOption = int(args[args.index('--groupValCalcOption')+1])
    groupValsStr = args[args.index('--groupVals')+1]
    groupVals = [float(v) for v in groupValsStr.split(',') if v != '']
    resFile = args[args.index('-o')+1]
    #pdb.set_trace()

    caller_descriptions, snpSum3 = load_snpSum3(snpSum3File)
    #pdb.set_trace()

    if '--ab_percentile' in args:
        ab_perc, perc_ab = load_ab_perc(args[args.index('--ab_percentile')+1])
    else:
        ab_perc = None
        perc_ab = None
    #pdb.set_trace()

    count_by_group(groupValCalcOption, groupVals, snpSum3File, snpSum3, caller_descriptions, resFile, ab_perc=ab_perc, perc_ab=perc_ab)
    #pdb.set_trace()

    return

'''
python evaluator.py --calcAbPercentile -1 snpSum3_1 [-2 snpSum3_2 etc] -o ab_percentile_file

ab_percentile_file format:
col-0: percentile e.g. 0~99
col-1: related ab val
'''

def calcAbPercentile(args):
    #pdb.set_trace()

    idx = 0
    files = []
    while '-%d'%idx in args:
        files.append(args[args.index('-%d'%idx)+1])
        idx += 1
    #pdb.set_trace()

    ab_values = []
    for file in files:
        with open(file, 'r') as f:
            for line in f:
                if line[0]=='#' or line.strip()=='': continue
                ab_values.append(float(line.split()[6]))
                
                #if line.split()[3]=='m': ab_values.append(float(line.split()[4]))
                #if line.split()[3]=='p': ab_values.append(float(line.split()[5]))

    #pdb.set_trace()

    ab_percentile_file = args[args.index('-o')+1]
    with open(ab_percentile_file, 'w') as f:
        for qt in xrange(0,100,1):
            thre = np.percentile(ab_values, qt)
            f.write('%d\t%.2f\n'%(qt, thre))
    #pdb.set_trace()

    print('%s written'%ab_percentile_file)

    return

'''
Description:
- use snp_analysis.py
- compare snp res among different callers and true snps
- loadSnpInfo: read snp info from true snps or called snps, add count/countAlt info if related count/countAlt files supplied (to check uniqueness)

Usage:

python evaluator.py    --loadSnpInfo 
                       -L1 label1 -F1 file1 [-L2 label2 -F2 file2 ...]
                       [-C1 countFile] [-C2 countAltFile]
                       [--snpLog outFile]
                       [--snpSum sumFile]
                       [--snpSum2 sumFile2]

python evaluator.py    --gen_snpSum3
                       -m true_m_snp
                       -p true_p_snp
                       -o abSNP_snp
                       -g gatk_snp
                       -c count_alt
                       -o snpSum3_file
                       [--rm true_read_m_bed]
                       [--rp true_read_p_bed]
'''

if __name__ == "__main__":

    args = sys.argv

    if '--loadSnpInfo' in args:

        snDic = SnpNodeDic()

        snDic.loadSnpInfo(args)

        if '--snpLog' in args:
            snDic.writeSnpLog(args)

        if '--snpSum' in args and '--snpSum2' in args:
            snDic.writeSnpSummary2(args) #snDic.writeSnpSummary(args)
    elif '--gen_snpSum3' in args:
        gen_snpSum3(args)
    elif '--sensitivity_analysis' in args:
        sensitivity_analysis(args)
    elif '--calcAbPercentile' in args:
        calcAbPercentile(args)
