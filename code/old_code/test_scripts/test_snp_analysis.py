import pdb, sys
from util import *

'''
test functions in snp_analysis.py

usage:

python test_snp_analysis.py --condition conditionStr 
                            --copyRes(make a copy to code/tmp/<note>/ or not) 
                            [--round n (not to be used)]
                            [--thre thre] 
                            --para
                            --compare 0(non-filt)/1(GATK)/2... 
                            --sam samfile
note:
compare: 0 - compare with non-filtered; 1 - compare with GATK; 2 etc - TBD -- info included as c0/c1 ... in snpLog, snpSum, snpSum2
samfile: helps decide count/countAlt/res naming -- info included as [sam_fn] in snpLog, snpSum, snpSum2

compare our called snp (filtered) with non-filtered or GATK
generate snpLog (count/countAlt Str etc for debug purpose), snpSum (md and fp), snpSum2 (related gPos of snpSum)

copy true snps, called snps, [count, countAlt not copied] files into code/tmp/<note: related to conditionStr [and round]>/
'''

def test_loadSnpFiles(args):

    #pdb.set_trace()

    if '--condition' in args:
        conditionStr = args[args.index('--condition')+1]
    else:
        conditionStr = ''

    if '--round' in args:
        n_th_round = args[args.index('--round')+1]
        note = conditionStr + '_round_' + n_th_round
    else:
        n_th_round = ''
        note = conditionStr

    if '--thre' in args:
        thre = int(args[args.index('--thre')+1])
        threStr = '_T%d'%thre
    else:
        thre = 1
        threStr = '_T1'

    #pdb.set_trace()
    samfile = args[args.index('--sam')+1]
    #sam_dir = '/'+'/'.join([itm for itm in samfile.split('/')[:-1] if itm != ''])+'/'
    sam_fn = samfile.split('/')[-1][:-4] #no .txt file type
    sam_dir = samfile[0:len(samfile)-len(sam_fn)-4]
    

    if '--para' in args:
        paraStr = '_para'
        #C1 = ['/data1/shunfu1/SNPCalling/data/data_GATK/2pass/split_sorted_sam_dedupped/count_y/count_dedupped.txt', 'count_dedupped.txt']
        #C2 = ['/data1/shunfu1/SNPCalling/data/data_GATK/2pass/split_sorted_sam_dedupped/count_y_altInfo/count_dedupped_altInfo.txt', 'count_dedupped_altInfo.txt']
        C1 = ['%s/split_sorted_sam_%s/count_y/count_%s.txt'%(sam_dir, sam_fn, sam_fn), 'count_%s.txt'%(sam_fn)]
        C2 = ['%s/split_sorted_sam_%s/count_y_altInfo/count_%s_altInfo.txt'%(sam_dir, sam_fn, sam_fn), 'count_%s_altInfo.txt'%(sam_fn)]
    else:
        paraStr = ''
        #C1 = ['/data1/shunfu1/SNPCalling/data/data_GATK/2pass/count_dedupped.txt', 'count_dedupped.txt']
        #C2 = ['/data1/shunfu1/SNPCalling/data/data_GATK/2pass/count_dedupped_altInfo.txt', 'count_dedupped_altInfo.txt']
        C1 = ['%s/count_%s.txt'%(sam_dir, sam_fn), 'count_%s.txt'%(sam_fn)]
        C2 = ['%s/count_%s_altInfo.txt'%(sam_dir, sam_fn), 'count_%s_altInfo.txt'%(sam_fn)]

    compare = int(args[args.index('--compare')+1])

    files = []
    #files.append(['O','/data1/shunfu1/SNPCalling/data/caller_output_snp_found_dedupped%s%s_filt.txt'%(threStr, paraStr), 'caller_output_snp_found_dedupped%s%s_filt.txt'%(threStr, paraStr)])
    files.append(['O','/data1/shunfu1/SNPCalling/data/caller_output_snp_found_%s%s%s_filt.txt'%(sam_fn, threStr, paraStr), 'caller_output_snp_found_%s%s%s_filt.txt'%(sam_fn, threStr, paraStr)])
    files.append(['m','/data1/shunfu1/SNPCalling/data/SNP_m.txt', 'SNP_m.txt'])
    files.append(['p','/data1/shunfu1/SNPCalling/data/SNP_p.txt', 'SNP_p.txt'])
    if compare==0:
        #files.append(['G','/data1/shunfu1/SNPCalling/data/caller_output_snp_found_dedupped%s%s.txt'%(threStr, paraStr), 'caller_output_snp_found_dedupped%s%s.txt'%(threStr, paraStr)])
        files.append(['G','/data1/shunfu1/SNPCalling/data/caller_output_snp_found_%s%s%s.txt'%(sam_fn, threStr, paraStr), 'caller_output_snp_found_%s%s%s.txt'%(sam_fn, threStr, paraStr)])
    elif compare==1:
        files.append(['G','/data1/shunfu1/SNPCalling/data/data_GATK/GATK_out/raw_variants.vcf.txt','raw_variants.vcf.txt'])
    else:
        print('unknown compare mode: %d'%compare)
        pdb.set_trace()
    
    snpLog = '/data1/shunfu1/SNPCalling/data/snpLog_[%s]%s_c%d.txt'%(sam_fn, threStr, compare)
    snpLog_fn = snpLog.split('/')[-1]
    snpSum = '/data1/shunfu1/SNPCalling/data/snpSum_[%s]%s_c%d.txt'%(sam_fn, threStr, compare)
    snpSum_fn = snpSum.split('/')[-1]
    snpSum2 = '/data1/shunfu1/SNPCalling/data/snpSum2_[%s]%s_c%d.txt'%(sam_fn, threStr, compare)
    snpSum2_fn = snpSum2.split('/')[-1]

    #run snp_analysis
    nFiles = len(files)

    cmd = 'python snp_analysis.py --loadSnpInfo '

    for i in range(nFiles):

        cmd += ' -L%d %s -F%d %s'%(i+1, files[i][0], i+1, files[i][1])

    if C1 != []:
        cmd += ' -C1 %s'%C1[0]

    if C2 != []:
        cmd += ' -C2 %s'%C2[0]

    if snpLog != '':
        cmd += ' --snpLog %s'%snpLog

    if snpSum != '':
        cmd += ' --snpSum %s'%snpSum

    if snpSum2 != '':
        cmd += ' --snpSum2 %s'%snpSum2

    pdb.set_trace()
    run_cmd(cmd)

    #organize res
    if '--copyRes' in args:
        print('copy res to tmp/%s/ ...'%note)
        #pdb.set_trace()

        cmd = 'mkdir -p tmp/%s'%note
        run_cmd(cmd)

        cmd = 'cp %s tmp/%s/%s'%(snpLog, note, snpLog_fn)
        run_cmd(cmd)

        cmd = 'cp %s tmp/%s/%s'%(snpSum, note, snpSum_fn)
        run_cmd(cmd)

        cmd = 'cp %s tmp/%s/%s'%(snpSum2, note, snpSum2_fn)
        run_cmd(cmd)

        #keep more results
        
        for file in files:
            cmd = 'cp %s tmp/%s/%s'%(file[1],note,file[2])
            run_cmd(cmd)

        #cmd = 'cp %s tmp/%s/%s'%(C1[0], note, C1[1])
        #run_cmd(cmd)

        #cmd = 'cp %s tmp/%s/%s'%(C2[0], note, C2[1])
        #run_cmd(cmd)

    return

if __name__ == "__main__":

    args = sys.argv

    test_loadSnpFiles(args)
