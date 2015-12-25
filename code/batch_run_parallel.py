"""
- 9/28: try to convert complex blocks to be run parallelly (tag: para_comp)
- 9/30: try to generate count file also parallelly
- 10/3: 0) no-cross-check & cross-check: depends on sam file used 

        1) non-parallel naming structure:
           data/sam_dir/<z>.sam & count_<z>.txt & count_<z>_altInfo.txt
           data/caller_output_(|snp|exception)_<z>_T[thre val].txt & snp_res_<z>_T[thre val].txt
           
        2) parallel naming structure:
           below are files split to be processed parallelly:
           data/sam_dir/split_sorted_sam_<z>/split_sam_[00~19]
                                            /count_x/count_x[00~19 *split_sam]_region_[00~19].txt
                                            /count_y/count_y[00~19 *region].txt
                                            /count_x_altInfo/count_x_altInfo[00~19 *split_sam]_region_[00~19].txt
                                            /count_y_altInfo/count_y_altInfo[00~19 *region].txt
                                                            /count_<z>_altInfo.txt
                                            /caller_op_T[thre val]/caller_op_(|snp|exception)_y[00~19].txt
           below are results:
           data/caller_output_(|snp|exception)_<z>_T[thre val]_para.txt & snp_res_<z>_T[thre val]_para.txt
           
        3) final_caller accepts abs & rel file address
- 10/14: try to dump additional info for count (file structure updated)
         analyze_fp
- 10/23: do_rsem
         -- generate expression estimation (rsem)
         -- generate count using rsem res (both non-parallel & parallel mode)
         -- do_rsem & not do_rsem currently have same name structure for output files
- 11/22: lambda half filtering
         -- filter_snp_lam_half_filt
"""

from Address import *
from Synthesis import *
from ReadProcess import *
from count_read_lambda import *
from final_caller19thjuly_m import *

from debug_MACRO import *
from sim_stat import *
from snp_res_statistics import *
from find_pos import *
from para_operations import *
import pdb
import subprocess

def batch_run_parallel(do_cross_check=False,
                       Threshold_num_reads=1):
                           
    #pdb.set_trace()
    
    para_comp = False #True
    num_p = 1 #20 #20 #number of threads available
    print('Enable parallel comp (tophat, counter, final_caller etc) is %s with num of threads=%d'%(para_comp, num_p))
    
    print('Default data folder is %s'%Default_Ref_Path)
    
    Stat_m = sim_stat('/sim_stat_dmp_m.txt')
    Stat_p = sim_stat('/sim_stat_dmp_p.txt')
    
    ## ---------- synthesis data
    ref_address = Default_Ref_Path + '/Chr15.fa'
    BED_address = Default_Ref_Path + '/hg19_chr15-UCSC.bed'
    
    EXP_fn_m = '/exp_m.txt' #for true exp
    EXP_fn_p = '/exp_p.txt'
    
    COV_fn_m = '/coverage_m.txt'
    COV_fn_p = '/coverage_p.txt'
    
    Num_SNP = 20 #1000 #20 #for true snp
    
    tar_address_m = Default_Ref_Path + '/Tar_m.txt'
    tar_address_p = Default_Ref_Path + '/Tar_p.txt'

    SNP_address_m = Default_Ref_Path + '/SNP_m.txt' 
    SNP_address_p = Default_Ref_Path + '/SNP_p.txt'
    
    N=100000 #10000000 #100000 #Number of Reads
    L=100       #Read Length
    error_rate = 0
    
    do_gen_exp = False
    if do_gen_exp==True:
        [BED_sorted_address, exp_address_m] = BED2ExpressionLevel(BED_address, EXP_fn_m) # generate a random expression level file
        [BED_sorted_address, exp_address_p] = BED2ExpressionLevel(BED_address, EXP_fn_p)
    else:
        exp_address_m = Default_Ref_Path + EXP_fn_m #for test purpose
        exp_address_p = Default_Ref_Path + EXP_fn_p
        BED_sorted_address = Default_Ref_Path + '/hg19_chr15-UCSC-sorted.bed'
    
    do_gen_cov = False
    if do_gen_cov==True:
        #old cov
        coverage_address_m = ExpressionLevel2Coverage(BED_sorted_address, exp_address_m, COV_fn_m, Stat_m) # Find the exonic parts with coverage not less than a ceryain value
        coverage_address_p = ExpressionLevel2Coverage(BED_sorted_address, exp_address_p, COV_fn_p, Stat_p)
    else:
        coverage_address_m = Default_Ref_Path + COV_fn_m #for test purpose
        coverage_address_p = Default_Ref_Path + COV_fn_p
    
    #pdb.set_trace()
    
    do_gen_fil_cov = False
    qt = [90]
    if do_gen_fil_cov==True:
        
        Stat_m.set_qt(qt) 
        Stat_m.get_acc_cov_hd(COV_fn_m)
        Stat_m.set_acc_cov_qt()
        line_cover_threshold = Stat_m.acc_cover_qt[0]
        COV_fn_filtered_m = COV_fn_m[:-4] + '_qt'+repr(qt[0])+'.txt'
        coverage_address_filtered_m = Default_Ref_Path + COV_fn_filtered_m
        #use filtered coverage.txt to make snps generated in high exp-level region
        FilterCoverage(coverage_address_m, coverage_address_filtered_m, line_cover_threshold) 
        
        Stat_p.set_qt(qt) 
        Stat_p.get_acc_cov_hd(COV_fn_p)
        Stat_p.set_acc_cov_qt()
        line_cover_threshold = Stat_p.acc_cover_qt[0]
        COV_fn_filtered_p = COV_fn_p[:-4] + '_qt'+repr(qt[0])+'.txt'
        coverage_address_filtered_p = Default_Ref_Path + COV_fn_filtered_p
        #use filtered coverage.txt to make snps generated in high exp-level region
        FilterCoverage(coverage_address_p, coverage_address_filtered_p, line_cover_threshold)
        
    else:
        
        COV_fn_filtered_m = COV_fn_m[:-4] + '_qt'+repr(qt[0])+'.txt'
        coverage_address_filtered_m = Default_Ref_Path + COV_fn_filtered_m
        
        COV_fn_filtered_p = COV_fn_p[:-4] + '_qt'+repr(qt[0])+'.txt'
        coverage_address_filtered_p = Default_Ref_Path + COV_fn_filtered_p
        
    do_gen_tar = False
    if do_gen_tar==True:   
        GenTarget(ref_address, coverage_address_filtered_m, Num_SNP, tar_address_m, SNP_address_m ) # Generate 2 random target2 (for paternal and maternal) sequence from ref_address by insering Num_SNP SNPs in the exonic positions
        GenTarget(ref_address, coverage_address_filtered_p, Num_SNP, tar_address_p, SNP_address_p )
    
    do_gen_reads = False
    if do_gen_reads==True:
        [readBED_address_m, readFA_address_m, readFQ_address_m] = ReadGeneration(tar_address_m, BED_address, exp_address_m,  N, L, error_rate)
        [readBED_address_p, readFA_address_p, readFQ_address_p] = ReadGeneration(tar_address_p, BED_address, exp_address_p,  N, L, error_rate)
        readFQ_address = Default_Ref_Path + '/Tar_read_l'+ repr(L) +'.fastq'  #for test purpose
        subprocess.call('cat '  + readFQ_address_m + ' ' + readFQ_address_p + ' > ' + readFQ_address, shell=True )
    else:
        readBED_address_m = Default_Ref_Path + '/Tar_m_read_l'+ repr(L) +'.bed'
        readFA_address_m = Default_Ref_Path + '/Tar_m_read_l'+ repr(L) +'.fasta'
        readFQ_address_m = Default_Ref_Path + '/Tar_m_read_l'+ repr(L) +'.fastq' #for test purpose
        
        readBED_address_p = Default_Ref_Path + '/Tar_p_read_l'+ repr(L) +'.bed'
        readFA_address_p = Default_Ref_Path + '/Tar_p_read_l'+ repr(L) +'.fasta'
        readFQ_address_p = Default_Ref_Path + '/Tar_p_read_l'+ repr(L) +'.fastq' #for test purpose
        
        readFQ_address = Default_Ref_Path + '/Tar_read_l'+ repr(L) +'.fastq'  #for test purpose
            
    #pdb.set_trace()
    ## ---------- read alignment & abundance estimation
    aligner ='tophat'
    gtf_address = Default_Ref_Path + '/hg19_chr15-UCSC.gtf'
    tophat_dir = '/tophat_out/'
    EXON_fn='/exon.txt'
    
    #en_debug_0811 sam_address = Align(ref_address, readFQ_address_m + ' ' + readFQ_address_p, aligner, Default_Ref_Path+tophat_dir)
    do_align = False
    if do_align==True:
        if para_comp==False:
            sam_address = Align(ref_address, readFQ_address, aligner, Default_Ref_Path+tophat_dir)
        else:
            sam_address = Align(ref_address, readFQ_address, aligner, Default_Ref_Path+tophat_dir, num_p)
    else:
        sam_address = Default_Ref_Path + tophat_dir + '/accepted_hits.sam' #for test purpose
    
    do_rsem = True
    if do_rsem==True:
        #pdb.set_trace()
        do_RSEM = False
        if do_RSEM == True:
            RSEM_result_address = RSEM(ref_address, readFQ_address, BED_address, gtf_address)
        else:
            RSEM_result_address = Default_Ref_Path + '/rsem/Chr15.isoforms.results'
        
        do_BED2Exon = False 
        if do_BED2Exon == True:
            exon_address = BED2Exon(BED_address, EXON_fn)
        else:
            exon_address = Default_Ref_Path + EXON_fn

        do_RSEM2Coverage = False
        if do_RSEM2Coverage == True:
            #pdb.set_trace()
            RSEM2Coverage(RSEM_result_address, exon_address,  L)
            count_rsem_address = Default_Ref_Path + '/count_rsem.txt'
        else:
            count_rsem_address = Default_Ref_Path + '/count_rsem.txt' #for test purpose
    else: #do_rsem==False
        count_rsem_address = Default_Ref_Path + '/count_rsem.txt'
    #cross check
    #do_cross_check = True
    if do_cross_check==True:
        
        flder = Default_Ref_Path +'/data_GATK/2pass/'
        name = 'dedupped.bam' #Aligned.out.sam
        
        #prepare sam file    
        if 'bam' in name:
            bam_address = flder + name
            sam_address = flder + name[:-4] + '.sam'
            
            if os.path.exists(sam_address)==False: #need to convert bam to sam
                SamtoolsProgram = SamPath + '/samtools'
                subprocess.call(SamtoolsProgram +  ' view -h ' + bam_address + ' > ' + sam_address, shell=True )        
        else:
            sam_address = flder + name
    
    ## ---------- generate count file
    #pdb.set_trace()
    do_count = False #True
    if do_count==True:
        if para_comp==True:
            
            #seperate sorted sam file into num_p smaller files, with reads of same read group in same file
            do_sep_sam = True
            
            sam_fn = sam_address.split('/')[-1]
            sam_dir = sam_address[0:len(sam_address)-len(sam_fn)]
            split_sam_dir = sam_dir + '/split_sorted_sam_'+sam_fn[:-4]+'/' #e.g. /split_sorted_sam_dedupped/
            
            if do_sep_sam == True:
                [split_sam_dir, split_sam_pre_fn, num_split_sam_files] = seperate_sam_file(sam_dir, sam_fn, split_sam_dir, num_p)
            else:
                split_sam_pre_fn = 'split_sam_'
                num_split_sam_files = num_p #should be
                
            #generate count files seperatedly based on split_sam (to clear split_sam)
            do_para_count = True
            count_split_dir0 = '/count_x/'
            count_split_dir1 = '/count_y/'
            
            if do_para_count == True:
                if True:
                    if do_rsem==False:
                        cmd = 'parallel python count_read_lambda.py '+\
                              split_sam_dir+' '+\
                              split_sam_pre_fn + '{} '+\
                              ref_address + ' ' +\
                              coverage_address + ' ' +\
                              count_split_dir0 + ' '+\
                              'count_x{}.txt '+\
                              repr(num_p)+' '+\
                              ':::'
                    else:
                        cmd = 'parallel python count_read_lambda.py '+\
                              split_sam_dir+' '+\
                              split_sam_pre_fn + '{} '+\
                              ref_address + ' ' +\
                              count_rsem_address + ' ' +\
                              count_split_dir0 + ' '+\
                              'count_x{}.txt '+\
                              repr(num_p)+' '+\
                              ':::'
                    for i in range(num_split_sam_files):
                        cmd = cmd + ' %02d'%i
                        
                    subprocess.call(cmd, shell = True)
                else: #for debug purpose
                    #pdb.set_trace()
                    for i in range(20):
                        tag = '%02d'%i
                        if do_rsem==False:
                            cmd = 'python count_read_lambda.py '+\
                                  split_sam_dir+' '+\
                                  split_sam_pre_fn + tag + ' '+\
                                  ref_address + ' ' +\
                                  coverage_address + ' ' +\
                                  count_split_dir0 + ' '+\
                                  'count_x'+tag+'.txt '+\
                                  repr(num_p)+' '
                        else:
                            cmd = 'python count_read_lambda.py '+\
                                  split_sam_dir+' '+\
                                  split_sam_pre_fn + tag + ' '+\
                                  ref_address + ' ' +\
                                  count_rsem_address + ' ' +\
                                  count_split_dir0 + ' '+\
                                  'count_x'+tag+'.txt '+\
                                  repr(num_p)+' '
                        subprocess.call(cmd, shell = True)
                #pdb.set_trace()  
                #subprocess.call('rm '+split_sam_dir+'/'+split_sam_pre_fn+'*', shell = True)
            
            #merge generated count files (to clear split_count)
            do_merge_count = True
            if do_merge_count == True:
                #pdb.set_trace()
                if True:
                    cmd = 'parallel python para_operations.py para_count_merge '+\
                          split_sam_dir+count_split_dir0+' '\
                          'count_x'+' '+\
                          'region_{}.txt'+' '+\
                          repr(num_p)+' '+\
                          split_sam_dir+count_split_dir1 + ' '+\
                          'count_y{}.txt '+\
                          ':::'
                    for i in range(num_split_sam_files):
                        cmd = cmd + ' %02d'%i
                    subprocess.call(cmd, shell = True)
                else: #for debug purpose
                    for i in range(20):
                        tag = '%02d'%i
                        cmd = 'python para_operations.py para_count_merge '+\
                              split_sam_dir+count_split_dir0+' '\
                              'count_x'+' '+\
                              'region_'+tag+'.txt'+' '+\
                              repr(num_p)+' '+\
                              split_sam_dir+count_split_dir1 + ' '+\
                              'count_y'+tag+'.txt'
                        subprocess.call(cmd, shell = True)
                        #pdb.set_trace()                          
                #pdb.set_trace()
                
            #pdb.set_trace()
            do_merge_count_altInfo = True
            if do_merge_count_altInfo==True:
                #pdb.set_trace()
                count_altInfo_address = split_sam_dir + '/count_y_altInfo/' + '/count_'+sam_fn[:-4]+'_altInfo.txt'
                merge_count_y_altInfo(split_sam_dir + '/count_y_altInfo/count_y_altInfo', num_p, count_altInfo_address)
                #pdb.set_trace()
                #merge_count_x(split_sam_dir + count_split_dir0,
                #              'count_x',
                #              num_split_sam_files,
                #              num_p,
                #              Default_Ref_Path + tophat_dir + count_split_dir1,
                #              'count_y')        
        else:
            #pdb.set_trace()
            
            sam_fn = sam_address.split('/')[-1]
            sam_dir = sam_address[0:len(sam_address)-len(sam_fn)]
            
            count_fn = '/count_'+sam_fn[:-4]+'.txt'
            if do_rsem==False:
                generate_count_file(sam_address, ref_address, coverage_address, count_fn) #en_debug_0817
            else:
                #pdb.set_trace()
                generate_count_file(sam_address, ref_address, count_rsem_address, count_fn)
            
                    
            count_abs_address = sam_dir + count_fn
            
            #pdb.set_trace()
    else:
        #print('count_rsem?')
        #pdb.set_trace()
        if para_comp==True:
            sam_fn = sam_address.split('/')[-1]
            sam_dir = sam_address[0:len(sam_address)-len(sam_fn)]
            split_sam_dir = sam_dir + '/split_sorted_sam_'+sam_fn[:-4]+'/' #e.g. /split_sorted_sam_dedupped/
            split_sam_pre_fn = 'split_sam_'
            num_split_sam_files = num_p #should be
                
            count_split_dir0 = '/count_x/'
            count_split_dir1 = '/count_y/'
            
            
        else:
            pdb.set_trace()
            
            sam_fn = sam_address.split('/')[-1]
            sam_dir = sam_address[0:len(sam_address)-len(sam_fn)]
            count_fn = '/count_'+sam_fn[:-4]+'.txt'
            count_abs_address = sam_dir + count_fn
            
    ## ---------- pre-processing for parallel operation
    do_caller = False
    #caller_op_fn = '/caller_output_para'
    #caller_op_exception_fn = '/caller_output_exception_para'
    #caller_op_snp_found_fn = '/caller_output_snp_found_para'
    #Threshold_num_reads = 10 #do SNP only if # of reads >= Threshold
    caller_op_fn = '/caller_output_'+sam_fn[:-4]+'_T'+repr(Threshold_num_reads)
    caller_op_exception_fn = '/caller_output_exception_'+sam_fn[:-4]+'_T'+repr(Threshold_num_reads)
    caller_op_snp_found_fn = '/caller_output_snp_found_'+sam_fn[:-4]+'_T'+repr(Threshold_num_reads)
    ft = '.txt' #file type
    
    if do_caller==True:
        if para_comp==True:
            
            count_split_dir = split_sam_dir+count_split_dir1
            caller_op_dir = split_sam_dir + '/caller_op_T' + repr(Threshold_num_reads) + '/'
            subprocess.call('mkdir -p '+caller_op_dir, shell=True)
            
            #pdb.set_trace()
            
            if True:
                cmd = 'parallel python final_caller19thjuly_m.py '+\
                      count_split_dir + ' '+\
                      'count_y{}.txt '+\
                      caller_op_dir+'caller_op_y{}.txt '+\
                      caller_op_dir+'caller_op_exception_y{}.txt '+\
                      caller_op_dir+'caller_op_snp_y{}.txt '+\
                      repr(Threshold_num_reads)+' '\
                      ':::'
                for i in range(num_p):
                    cmd = cmd + ' %02d'%i
            else:
                cmd = 'python final_caller19thjuly_m.py '+\
                      count_split_dir + ' '+\
                      'count_y00.txt '+\
                      caller_op_dir+'caller_op_y00.txt '+\
                      caller_op_dir+'caller_op_exception_y00.txt '+\
                      caller_op_dir+'caller_op_snp_y00.txt '+\
                      repr(Threshold_num_reads)
            #pdb.set_trace()              
            subprocess.call(cmd, shell = True)
            #parallel python final_caller19thjuly_m.py x{} caller_op_x{}.txt caller_op_exception_x{}.txt caller_op_snp_x{}.txt ::: 08 10 18 #00 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19
            #pdb.set_trace()
            #subprocess.call('rm '+count_split_dir+'/x*', shell = True)
            
            #pdb.set_trace()
            
            caller_op_fn += '_para'
            caller_op_exception_fn += '_para'
            caller_op_snp_found_fn += '_para'
                
            #pdb.set_trace()
    
            merge_caller_res(Default_Ref_Path,
                             caller_op_dir+'caller_op_y',
                             caller_op_dir+'caller_op_exception_y',
                             caller_op_dir+'caller_op_snp_y',
                             num_p,
                             caller_op_fn+ft,
                             caller_op_exception_fn+ft,
                             caller_op_snp_found_fn+ft)
            
        else: #final caller, non-parallel
    ## ---------- caller   
            #pdb.set_trace()
            final_caller(count_abs_address,
                         caller_op_fn+ft,
                         caller_op_exception_fn+ft,
                         caller_op_snp_found_fn+ft,
                         Threshold_num_reads)
    else:
        caller_op_snp_found_fn = '/caller_output_snp_found_'+sam_fn[:-4]+'_T'+repr(Threshold_num_reads)
        if para_comp==True:
            caller_op_snp_found_fn += '_para'
    
    ##---------- snp res
    #pdb.set_trace()
    do_snp_res_stat = False #True
    snp_res_stat = []
    if do_snp_res_stat==True:
        #snp_res_address = Default_Ref_Path + '/data_GATK/GATK_out/raw_variants.vcf'
        snp_res_address = Default_Ref_Path + caller_op_snp_found_fn + ft
        snp_res_stat_fn = Default_Ref_Path + 'snp_res_' + sam_fn[:-4]+'_T'+repr(Threshold_num_reads)
        if para_comp==True:
            snp_res_stat_fn = snp_res_stat_fn + '_para'+ft
        else:
            snp_res_stat_fn = snp_res_stat_fn + ft
        snp_res_stat = do_snp_res_statistics(snp_res_address, SNP_address_m, SNP_address_p, snp_res_stat_fn)
        
    do_fp_analysis = False
    if do_fp_analysis == True and snp_res_stat != []:
        
        m_snps_cd = snp_res_stat[0]
        m_snps_md = snp_res_stat[1]
        p_snps_cd = snp_res_stat[2]
        p_snps_md = snp_res_stat[3]
        r_snps_fp = snp_res_stat[4]
        
        if para_comp==True:
            #pdb.set_trace()
            
            sam_fn = sam_address.split('/')[-1]
            sam_dir = sam_address[0:len(sam_address)-len(sam_fn)]
            split_sam_dir = sam_dir + '/split_sorted_sam_'+sam_fn[:-4]+'/'
            count_altInfo_address = split_sam_dir + '/count_y_altInfo/' + '/count_'+sam_fn[:-4]+'_altInfo.txt'
            
            analyze_fp(m_snps_cd,
                       m_snps_md,
                       p_snps_cd,
                       p_snps_md,
                       r_snps_fp,
                       count_altInfo_address)
        
            pdb.set_trace()                        
            
        else:
            pdb.set_trace()
            
            sam_fn = sam_address.split('/')[-1]
            sam_dir = sam_address[0:len(sam_address)-len(sam_fn)]
            count_altInfo_address = sam_dir + '/count_'+sam_fn[:-4]+'_altInfo.txt'
            
            analyze_fp(m_snps_cd,
                       m_snps_md,
                       p_snps_cd,
                       p_snps_md,
                       r_snps_fp,
                       count_altInfo_address)
        
            pdb.set_trace()
    #check snp statistics based counts_alt etc
    
    pdb.set_trace()    
    
    do_filt_snp = True
    if do_filt_snp == True:
        #prepare input/output file names
        snp_res_address = Default_Ref_Path + caller_op_snp_found_fn + ft
        #filt_snp_res_address = Default_Ref_Path + caller_op_snp_found_fn + '_filt' + ft
        
        #lam half filtering
        filt_snp_res_address2 = Default_Ref_Path + caller_op_snp_found_fn + '_filt' + ft
        if para_comp==True:
            sam_fn = sam_address.split('/')[-1]
            sam_dir = sam_address[0:len(sam_address)-len(sam_fn)]
            split_sam_dir = sam_dir + '/split_sorted_sam_'+sam_fn[:-4]+'/'
            count_altInfo_address = split_sam_dir + '/count_y_altInfo/' + '/count_'+sam_fn[:-4]+'_altInfo.txt'
            
            print('need count file also')
            pdb.set_trace()
        else:
            sam_fn = sam_address.split('/')[-1]
            sam_dir = sam_address[0:len(sam_address)-len(sam_fn)]
            count_altInfo_address = sam_dir + '/count_'+sam_fn[:-4]+'_altInfo.txt'
            #count_abs_address
        
        pdb.set_trace()
        
        #filter_snp(snp_res_address,
        #           filt_snp_res_address,
        #           count_altInfo_address,
        #           count_abs_address)
        
        filter_snp_lam_half_filt(snp_res_address,
                                 filt_snp_res_address2,
                                 count_altInfo_address,
                                 count_abs_address)
    else:
        pdb.set_trace()
    #end if do_filt_snp
    
    do_snp_res_stat2 = True
    snp_res_stat = []
    if do_snp_res_stat2==True:
        #pdb.set_trace()
        #snp_res_address = Default_Ref_Path + '/data_GATK/GATK_out/raw_variants.vcf'
        snp_res_address2 = Default_Ref_Path + caller_op_snp_found_fn + '_filt' + ft
        #snp_res_address2 = Default_Ref_Path + caller_op_snp_found_fn + '_lamHalfFilt' + ft
        snp_res_stat_fn2 = Default_Ref_Path + 'snp_res_' + sam_fn[:-4]+'_T'+repr(Threshold_num_reads)
        if para_comp==True:
            snp_res_stat_fn2 = snp_res_stat_fn2 + '_para_filt'+ft
        else:
            snp_res_stat_fn2 = snp_res_stat_fn2 + '_filt' + ft
        snp_res_stat = do_snp_res_statistics(snp_res_address2, SNP_address_m, SNP_address_p, snp_res_stat_fn2)
    
    pdb.set_trace()    
    
    return
        
if __name__ == "__main__":
    
    if len(sys.argv)>=2:
        do_cross_check= False
        if int(sys.argv[1])==1:
            do_cross_check= True
    
        Threshold_num_reads = int(sys.argv[2])
    else:
        do_cross_check= False #False
        Threshold_num_reads = 1
    
    #pdb.set_trace()
    batch_run_parallel(do_cross_check,
                       Threshold_num_reads)