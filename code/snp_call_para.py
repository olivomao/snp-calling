import sys, pdb

from old_code.util import run_cmd
'''
by ...1.py, we mainly modified:
- rel path to external absolute path
- a called ...py file will be replaced by code_folder/...py
'''
from old_code.batch_run_parallel_modi import do_sep_sam, do_para_count, \
                                             do_merge_count_x_and_altInfo, \
                                             do_merge_count_y_and_altInfo, \
                                             do_para_caller, \
                                             do_merge_para_caller_res, \
                                             do_filt3_snp_para

'''
modified based on batch_run_parallel.py
- cleaner codes
- esp, count_y altInfo merge and snp filt modified based on existing (or new) count altInfo format

assume:
- sam (including GATK best practices) already generated

usage: 

python snp_call_para.py  -r refAddress 
                         -c covAddress
                         -s inputSam(i.e. need to be grouped/sorted by read ids)
                         -O outDir(not used)
                         -T threshold
                         -p num_p(i.e. parallel, >1 recommended)
                         [--dupRun 1/0]
                         [--code_folder code_folder]
                         [--filt3_sameAb 1/0 (default 1)] [--filt3_rocT v1,v2,... (default 0,1)]
                         [--clear_intFiles 1/0 (default 0)]

param:

dupRun: 1 - is duplicated run, sam sep/ count made/ just need to do caller using different thresholds; default 0

code_folder: old_code/ or else

use filt3: filt similar ab + poisson filtering (w/ rocT)

output: outDir/caller_output_snp_found_<sam fn>_T<thre>_para.txt,
               caller_output_snp_found_<sam fn>_T<thre>_para_rocT_<rocT>.txt

'''

def main(args):

    #pdb.set_trace()

    ref_address = args[args.index('-r')+1]
    covAddress = args[args.index('-c')+1]
    sam_address = args[args.index('-s')+1]
    sam_fn = sam_address.split('/')[-1]
    outDir = args[args.index('-O')+1]
    #pdb.set_trace()
    run_cmd('mkdir -p %s'%outDir)
    Threshold_num_reads = int(args[args.index('-T')+1])
    num_p = int(args[args.index('-p')+1])


    if '--dupRun' in args:
        dupRun = int(args[args.index('--dupRun')+1])
    else:
        dupRun = 0

    if dupRun == 1:
        dupRun = True
    else:
        dupRun = False

    if '--code_folder' in args:
        code_folder = args[args.index('--code_folder')+1]
    else:
        code_folder = 'old_code/'

    #snp filt3
    if '--filt3_sameAb' in args:
        filtSameAb = int(args[args.index('--filt3_sameAb')+1])
        if filtSameAb==1:
            filtSameAb = True
        else:
            filtSameAb = False
    else:
        filtSameAb = True #default enable

    if '--filt3_rocT' in args:
        rocT = args[args.index('--filt3_rocT')+1].split(',')
        rocT = [float(v) for v in rocT if v != '']
    else:
        rocT = [0.0, 1.0] #(min fp, def)

    if '--clear_intFiles' in args:
        clear_intFiles = int(args[args.index('--clear_intFiles')+1])
    else:
        clear_intFiles = 0 

    #pdb.set_trace()
    [split_sam_dir, split_sam_pre_fn, num_split_sam_files] = do_sep_sam(sam_address, num_p, flag=(True and not dupRun))

    #pdb.set_trace()
    [count_split_dir0, count_split_dir1] = do_para_count( split_sam_dir,
                                                          split_sam_pre_fn,
                                                          num_split_sam_files,
                                                          ref_address,
                                                          covAddress,
                                                          num_p,
                                                          flag=(True and not dupRun),
                                                          code_folder=code_folder,
                                                          clear=clear_intFiles)

    #pdb.set_trace()
    do_merge_count_x_and_altInfo(split_sam_dir,
                                 count_split_dir0,
                                 count_split_dir1,
                                 num_split_sam_files,
                                 num_p,
                                 flag=(True and not dupRun),
                                 code_folder=code_folder,
                                 clear=clear_intFiles)    

    #pdb.set_trace()
    [count_abs_address, count_altInfo_address] = do_merge_count_y_and_altInfo(split_sam_dir,
                                                                              sam_fn,
                                                                              num_p,
                                                                              flag=(True and not dupRun),
                                                                              clear=0) #no clear, to be used by para_caller

    #pdb.set_trace()
    [para_caller_op_pre1, para_caller_op_pre2, para_caller_op_pre3] = do_para_caller( split_sam_dir,
                                                                                       sam_fn,
                                                                                       Threshold_num_reads,
                                                                                       num_p,
                                                                                       flag=True,
                                                                                       code_folder=code_folder,
                                                                                       clear=clear_intFiles)

    #pdb.set_trace()
    [caller_op_file, caller_op_exception_file, caller_op_snp_found_file] = do_merge_para_caller_res( sam_fn,
                                                                                                     para_caller_op_pre1, #caller_op_dir+'caller_op_y',
                                                                                                     para_caller_op_pre2, #caller_op_dir+'caller_op_exception_y',
                                                                                                     para_caller_op_pre3, #caller_op_dir+'caller_op_snp_y',
                                                                                                     Threshold_num_reads,
                                                                                                     num_p,
                                                                                                     flag=True,
                                                                                                     outDir=outDir,
                                                                                                     clear=clear_intFiles)

    do_filt3_snp_para(caller_op_snp_found_file,
                      count_abs_address,
                      count_altInfo_address,
                      filtSameAb,
                      rocT,
                      flag=True)

    return

'''
usage:

# count_generation + snp_call + snp_filtering, compute in parallel; see main() for details
python snp_call_para.py  -r refAddress 
                         -c covAddress
                         -s inputSam(i.e. need to be grouped/sorted by read ids)
                         -O outDir(not used)
                         -T threshold
                         -p num_p(i.e. parallel, >1 recommended)
                         [--dupRun 1/0]
                         [--code_folder code_folder]
                         [--filt3_sameAb 1/0 (default 1)] [--filt3_rocT v1,v2,... (default 0,1)]

'''
if __name__ == "__main__":
    
    args = sys.argv

    main(args)