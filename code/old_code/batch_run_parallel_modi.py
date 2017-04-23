'''
modify batch_run_parallel.py
- cleaner codes
- esp, count_y altInfo merge and snp filt modified based on existing (or new) count altInfo format

assume:
- sam (including GATK best practices) already generated

usage: 

python batch_run_paralle_modi.py -r refAddress -c covAddress -s inputSam(i.e. need to be grouped/sorted by read ids) 
                                 -O outDir(not used) -T threshold -p num_p(i.e. parallel, >1 recommended) [--dupRun 1/0]

param:

dupRun: 1 - is duplicated run, sam sep/ count made/ just need to do caller using different thresholds; default 0

output: outDir/caller_output_snp_found_dedupped_T<thre>_para.txt and caller_output_snp_found_dedupped_T<thre>_para_filt.txt

'''

import sys
from para_operations import *
from util import *
from Address import *
from final_caller19thjuly_m import *
from snp_res_statistics import * #filter_snp_lam_half_filt3

def do_sep_sam(sam_address, num_p, flag=True):    

    sam_fn = sam_address.split('/')[-1]
    sam_dir = sam_address[0:len(sam_address)-len(sam_fn)]
    split_sam_dir = sam_dir + '/split_sorted_sam_'+sam_fn[:-4]+'/' #e.g. /split_sorted_sam_dedupped/, in same directory of sam file

    if flag==False:
        return [split_sam_dir, 'split_sam_', num_p]
            
    [split_sam_dir, split_sam_pre_fn, num_split_sam_files] = seperate_sam_file(sam_dir, sam_fn, split_sam_dir, num_p)

    return [split_sam_dir, split_sam_pre_fn, num_split_sam_files]

#cov_address: count_rsem.txt or coverage.txt
#
#output: count_x/count_x{?? * corresponds to sub sam file}_region_{?? * corresponds to genome region}.txt
#        count_x_altInfo/count_x_altInfo{?? * corresponds to sub sam file}_region_{?? * corresponds to genome region}.txt
def do_para_count(split_sam_dir,
                  split_sam_pre_fn,
                  num_split_sam_files,
                  ref_address,
                  cov_address,
                  num_p,
                  flag=True,
                  code_folder='',
                  clear=0,
                  num_q=-1):

    count_split_dir0 = '/count_x/'
    count_split_dir1 = '/count_y/'

    if flag==False:
        return [count_split_dir0, count_split_dir1]    
    
    if code_folder != '': code_folder += '/'

    if num_q==-1:
        num_q = num_p

    num_batches = num_p / num_q
    for ith_batch in range(num_batches):

        cmd = 'parallel python %scount_read_lambda.py '%code_folder+\
              split_sam_dir+' '+\
              split_sam_pre_fn + '{} '+\
              ref_address + ' ' +\
              cov_address + ' ' +\
              count_split_dir0 + ' '+\
              'count_x{}.txt '+\
              repr(num_p)+' '+\
              ':::'

        for i in xrange(ith_batch*num_q, (ith_batch+1)*num_q):
            cmd = cmd + ' %02d'%i

        #pdb.set_trace()
        run_cmd(cmd)

    if clear == 1:        
        for i in range(num_split_sam_files):
            split_sam_file = '%s/%s%02d'%(split_sam_dir, split_sam_pre_fn, i)
            cmd = 'rm %s'%split_sam_file
            #pdb.set_trace()
            run_cmd(cmd)

    return [count_split_dir0, count_split_dir1]

#merge count_x/count_x{?? * corresponds to sub sam file}_region_{?? * corresponds to genome region}.txt
#into  count_y/count_y{?? * corresponds to genome region}.txt
#(see merge_count_x2 and merge_line_list2 for details)
#
#also
#merge count_x_altInfo/count_x_altInfo{?? * corresponds to sub sam file}_region_{?? * corresponds to genome region}.txt
#into  count_y_altInfo/count_y_altInfo{?? * corresponds to genome region}.txt
#(see merge_count_x_altInfo and merge_line_list_altInfo2 for details)

def do_merge_count_x_and_altInfo(split_sam_dir,
                                 count_split_dir0,
                                 count_split_dir1,
                                 num_split_sam_files,
                                 num_p,
                                 flag=True,
                                 code_folder='',
                                 clear=0):

    if flag==False:
        return

    if code_folder != '': code_folder+='/'
    cmd = 'parallel python %spara_operations.py para_count_merge '%code_folder+\
          split_sam_dir+count_split_dir0+' '\
          'count_x'+' '+\
          'region_{}.txt'+' '+\
          repr(num_p)+' '+\
          split_sam_dir+count_split_dir1 + ' '+\
          'count_y{}.txt '+\
          ':::'
    for i in range(num_split_sam_files):
        cmd = cmd + ' %02d'%i
    #pdb.set_trace()
    run_cmd(cmd)

    if clear==1:
        for i in range(num_split_sam_files):
            for j in range(num_split_sam_files):

                xf = '%s/%s/count_x%02d_region_%02d.txt'%(split_sam_dir, count_split_dir0, i, j)
                cmd = 'rm %s'%xf
                #pdb.set_trace()
                run_cmd(cmd)

                #xf_alt = '%s/%s_altInfo/count_x_altInfo%02d_region_%02d.txt'%(split_sam_dir, count_split_dir0, i, j)
                xf_alt = '%s/count_x_altInfo/count_x_altInfo%02d_region_%02d.txt'%(split_sam_dir, i, j)
                cmd = 'rm %s'%xf_alt
                run_cmd(cmd)

    return

#merge count_y/count_y{?? * corresponds to genome region}.txt
#into  count_y/count_<sam fn>.txt
#
#also
#merge count_y_altInfo/count_y_altInfo{?? * corresponds to genome region}.txt
#into  count_y_altInfo/count_<sam fn>_altInfo.txt
#
def do_merge_count_y_and_altInfo(split_sam_dir,
                                 sam_fn,
                                 num_p,
                                 flag=True,
                                 clear=0):

    count_abs_address     = split_sam_dir + '/count_y/'+          '/count_'+sam_fn[:-4]+'.txt'
    count_altInfo_address = split_sam_dir + '/count_y_altInfo/' + '/count_'+sam_fn[:-4]+'_altInfo.txt'

    if flag==False:
        return [count_abs_address, count_altInfo_address]
    
    for i in range(num_p):
        src = split_sam_dir + '/count_y/' + 'count_y%02d.txt'%i
        cmd = 'cat %s >> %s'%(src, count_abs_address)
        run_cmd(cmd)

    merge_count_y_altInfo(split_sam_dir + '/count_y_altInfo/count_y_altInfo', num_p, count_altInfo_address)

    if clear==1:
        for i in range(num_p):
            src = split_sam_dir + '/count_y/' + 'count_y%02d.txt'%i
            cmd = 'rm %s'%src
            #pdb.set_trace()
            run_cmd(cmd)

            src_alt = split_sam_dir + '/count_y_altInfo/' + 'count_y_altInfo%02d.txt'%i
            cmd = 'rm %s'%src_alt
            #pdb.set_trace()
            run_cmd(cmd)

    return [count_abs_address, count_altInfo_address]

#note: final_caller19thjuly_m.py uses absolute input paths now
#call count_y/count_y{?? * corresponds to genome region}.txt
#res into: caller_op_T<thre>/caller_op_y{?? * corresponds to genome region}.txt
#                           /caller_op_exception_y{?? * corresponds to genome region}.txt
#                           /caller_op_snp_y{?? * corresponds to genome region}.txt
def do_para_caller(split_sam_dir,
                   sam_fn,
                   Threshold_num_reads,
                   num_p,
                   flag=True,
                   code_folder='',
                   clear=0):
        
    count_split_dir = split_sam_dir+  '/count_y/'
    caller_op_dir  =  split_sam_dir + '/caller_op_T' + repr(Threshold_num_reads) + '/'
    run_cmd('mkdir -p '+caller_op_dir)

    para_caller_op_pre1 = caller_op_dir+'caller_op_y'
    para_caller_op_pre2 = caller_op_dir+'caller_op_exception_y'
    para_caller_op_pre3 = caller_op_dir+'caller_op_snp_y'

    if flag==False:
        return [para_caller_op_pre1, para_caller_op_pre2, para_caller_op_pre3]
    
    if code_folder != '': code_folder+='/'
    cmd = 'parallel python %sfinal_caller19thjuly_m.py '%code_folder+\
          count_split_dir + ' '+\
          'count_y{}.txt '+\
          caller_op_dir+'caller_op_y{}.txt '+\
          caller_op_dir+'caller_op_exception_y{}.txt '+\
          caller_op_dir+'caller_op_snp_y{}.txt '+\
          repr(Threshold_num_reads)+' '\
          ':::'
    for i in range(num_p):
        cmd = cmd + ' %02d'%i
    run_cmd(cmd)

    if clear==1:
        for i in range(num_p):
            src = split_sam_dir + '/count_y/' + 'count_y%02d.txt'%i
            cmd = 'rm %s'%src
            #pdb.set_trace()
            run_cmd(cmd)

            src_alt = split_sam_dir + '/count_y_altInfo/' + 'count_y_altInfo%02d.txt'%i
            cmd = 'rm %s'%src_alt
            #pdb.set_trace()
            run_cmd(cmd)

    return [para_caller_op_pre1, para_caller_op_pre2, para_caller_op_pre3]

def do_merge_para_caller_res(sam_fn,
                             para_caller_op_pre1, #caller_op_dir+'caller_op_y',
                             para_caller_op_pre2, #caller_op_dir+'caller_op_exception_y',
                             para_caller_op_pre3, #caller_op_dir+'caller_op_snp_y',
                             Threshold_num_reads,
                             num_p,
                             flag=True,
                             outDir=Default_Ref_Path,
                             clear=0):

    caller_op_file = outDir + '/caller_output_'+sam_fn[:-4]+'_T'+repr(Threshold_num_reads) + '_para.txt'
    caller_op_exception_file = outDir + '/caller_output_exception_'+sam_fn[:-4]+'_T'+repr(Threshold_num_reads) + '_para.txt'
    caller_op_snp_found_file = outDir + '/caller_output_snp_found_'+sam_fn[:-4]+'_T'+repr(Threshold_num_reads) + '_para.txt'

    if flag==False:
        return [caller_op_file, caller_op_exception_file, caller_op_snp_found_file]

    #pdb.set_trace()
        
    merge_caller_res(para_caller_op_pre1,
                     para_caller_op_pre2,
                     para_caller_op_pre3,
                     num_p,
                     caller_op_file,
                     caller_op_exception_file,
                     caller_op_snp_found_file,
                     outDir,
                     clear)

    return [caller_op_file, caller_op_exception_file, caller_op_snp_found_file]

'''
def do_merge_para_caller_res(sam_fn,
                             para_caller_op_pre1, #caller_op_dir+'caller_op_y',
                             para_caller_op_pre2, #caller_op_dir+'caller_op_exception_y',
                             para_caller_op_pre3, #caller_op_dir+'caller_op_snp_y',
                             Threshold_num_reads,
                             num_p,
                             flag=True):

    caller_op_file = Default_Ref_Path + '/caller_output_'+sam_fn[:-4]+'_T'+repr(Threshold_num_reads) + '_para.txt'
    caller_op_exception_file = Default_Ref_Path + '/caller_output_exception_'+sam_fn[:-4]+'_T'+repr(Threshold_num_reads) + '_para.txt'
    caller_op_snp_found_file = Default_Ref_Path + '/caller_output_snp_found_'+sam_fn[:-4]+'_T'+repr(Threshold_num_reads) + '_para.txt'

    if flag==False:
        return [caller_op_file, caller_op_exception_file, caller_op_snp_found_file]
        
    merge_caller_res(Default_Ref_Path,
                     para_caller_op_pre1,
                     para_caller_op_pre2,
                     para_caller_op_pre3,
                     num_p,
                     caller_op_file,
                     caller_op_exception_file,
                     caller_op_snp_found_file)

    return [caller_op_file, caller_op_exception_file, caller_op_snp_found_file]

#use external output files
def do_merge_para_caller_res1(sam_fn,
                              para_caller_op_pre1, #caller_op_dir+'caller_op_y',
                              para_caller_op_pre2, #caller_op_dir+'caller_op_exception_y',
                              para_caller_op_pre3, #caller_op_dir+'caller_op_snp_y',
                              Threshold_num_reads,
                              num_p,
                              outDir,
                              flag=True):

    caller_op_file = outDir + '/caller_output_'+sam_fn[:-4]+'_T'+repr(Threshold_num_reads) + '_para.txt'
    caller_op_exception_file = outDir + '/caller_output_exception_'+sam_fn[:-4]+'_T'+repr(Threshold_num_reads) + '_para.txt'
    caller_op_snp_found_file = outDir + '/caller_output_snp_found_'+sam_fn[:-4]+'_T'+repr(Threshold_num_reads) + '_para.txt'

    if flag==False:
        return [caller_op_file, caller_op_exception_file, caller_op_snp_found_file]

    pdb.set_trace()
        
    merge_caller_res1(para_caller_op_pre1,
                      para_caller_op_pre2,
                      para_caller_op_pre3,
                      num_p,
                      caller_op_file,
                      caller_op_exception_file,
                      caller_op_snp_found_file)

    return [caller_op_file, caller_op_exception_file, caller_op_snp_found_file]
'''



#modified based on do_filt_snp from batch_run_case1plus6
#
def do_filt_snp_para( sam_address,
                      caller_op_snp_found_addr,
                      count_abs_address,
                      count_altInfo_address,
                      flag=True):
    
    caller_op_snp_found_addr = caller_op_snp_found_addr.split('/')[-1]
    
    if flag == True:
        #prepare input/output file names
        snp_res_address = Default_Ref_Path + caller_op_snp_found_addr
        #lam half filtering
        filt_snp_res_address2 = Default_Ref_Path + caller_op_snp_found_addr[:-4] + '_filt.txt'
        
        sam_fn = sam_address.split('/')[-1]
        sam_dir = sam_address[0:len(sam_address)-len(sam_fn)]
        #count_altInfo_address = sam_dir + '/count_'+sam_fn[:-4]+'_altInfo.txt'
        
        #altCount modified
        filter_snp_lam_half_filt2(snp_res_address,
                                  filt_snp_res_address2,
                                  count_altInfo_address,
                                  count_abs_address)  
    else:
        filt_snp_res_address2 = Default_Ref_Path + caller_op_snp_found_addr[:-4] + '_filt.txt'
        
    return [filt_snp_res_address2]

#use external address
#depending on rocT, several filt snps (snp_res_address[:-4]+'_rocT_<rocT>.txt') will be generated
def do_filt3_snp_para(caller_op_snp_found_addr,
                      count_abs_address,
                      count_altInfo_address,
                      filtSameAb=True,
                      rocT=[0.0, 1.0],
                      flag=True):
    
    if flag == True:
        #prepare input/output file names
        snp_res_address = caller_op_snp_found_addr
        
        filter_snp_lam_half_filt3(snp_res_address,
                                  #filt_snp_res_address,
                                  count_altInfo_address,
                                  count_abs_address,
                                  filtSameAb=filtSameAb,
                                  rocThre=rocT)
    #pdb.set_trace()
    return

def main(args): #args: list of arguments

    ref_address = args[args.index('-r')+1]
    covAddress = args[args.index('-c')+1]
    sam_address = args[args.index('-s')+1]
    sam_fn = sam_address.split('/')[-1]
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

    #pdb.set_trace()
    [split_sam_dir, split_sam_pre_fn, num_split_sam_files] = do_sep_sam(sam_address, num_p, flag=(True and not dupRun))

    #pdb.set_trace()
    [count_split_dir0, count_split_dir1] = do_para_count( split_sam_dir,
                                                          split_sam_pre_fn,
                                                          num_split_sam_files,
                                                          ref_address,
                                                          covAddress,
                                                          num_p,
                                                          flag=(True and not dupRun))

    #pdb.set_trace()
    do_merge_count_x_and_altInfo(split_sam_dir,
                                 count_split_dir0,
                                 count_split_dir1,
                                 num_split_sam_files,
                                 num_p,
                                 flag=(True and not dupRun))    

    #pdb.set_trace()
    [count_abs_address, count_altInfo_address] = do_merge_count_y_and_altInfo(split_sam_dir,
                                                                              sam_fn,
                                                                              num_p,
                                                                              flag=(True and not dupRun))

    #pdb.set_trace()
    [para_caller_op_pre1, para_caller_op_pre2, para_caller_op_pre3] = do_para_caller(  split_sam_dir,
                                                                                         sam_fn,
                                                                                         Threshold_num_reads,
                                                                                         num_p,
                                                                                         flag=True)

    #pdb.set_trace()
    [caller_op_file, caller_op_exception_file, caller_op_snp_found_file] = do_merge_para_caller_res( sam_fn,
                                                                                                     para_caller_op_pre1, #caller_op_dir+'caller_op_y',
                                                                                                     para_caller_op_pre2, #caller_op_dir+'caller_op_exception_y',
                                                                                                     para_caller_op_pre3, #caller_op_dir+'caller_op_snp_y',
                                                                                                     Threshold_num_reads,
                                                                                                     num_p,
                                                                                                     flag=True,
                                                                                                     outDir = Default_Ref_Path)
    
    print('do_filt_snp_para uses filt2 in batch_run_parallel_modi')                                                                                                 
    pdb.set_trace()
    [caller_op_snp_found_file2] = do_filt_snp_para(   sam_address,
                                                      caller_op_snp_found_file,
                                                      count_abs_address,
                                                      count_altInfo_address,
                                                      flag=False)


    #pdb.set_trace()
    return

if __name__ == "__main__":

    #pdb.set_trace()

    if len(sys.argv)<2: #test mode

        ref_address = '/data1/shunfu1/SNPCalling/data/Chr15.fa'
        covAddress = '/data1/shunfu1/SNPCalling/data/count_rsem.txt'
        sam_address = '/data1/shunfu1/SNPCalling/data/data_GATK/2pass/dedupped.sam'
        output_address = '/data1/shunfu1/SNPCalling/data/' #not used yet


        testArg = '-r %s -c %s -s %s -O %s -T 1 -p 20'%(ref_address, covAddress, sam_address, output_address)
        main(testArg.split())

    else:

        main(sys.argv)
