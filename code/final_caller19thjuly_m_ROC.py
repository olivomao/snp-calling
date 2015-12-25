# -*- coding: utf-8 -*-
"""
Created on Sat Jul 18 12:12:03 2015

@author: Jatindeep
"""


import math
import csv
import operator
import pickle
from Address import *
import pdb
from debug_MACRO import *
from progress.bar import Bar
from count_read_lambda import *
import sys
import numpy as np
from snp_res_statistics import *

#pdb.set_trace()

def JP(Qi,Na,Nc,Ng,Nt,R,Lamda_Knot, Ly, Lamda,X,Y):
    
    #pdb.set_trace()
    #print('JP, X=%s'%X)
    
    N=Na+Nc+Ng+Nt
    Prob_log=0   
    
    for q in range(0,N):
        
        #if X=='A' and Y[q]=='G':
        #    pdb.set_trace()
        Prob_logq =Prob_funct(Qi[q],R,Lamda_Knot, Ly[q], Lamda[q],X,Y[q])
        
        Prob_log = Prob_log + Prob_logq
        
        #print('...%d/%d: X=%s, y=%s, Prob_logq:%f, Prob_log(acc):%f'%(q, N, X, Y[q], Prob_logq, Prob_log))
    #print('')
#        if(i==N-1):
    return Prob_log
            
        
def Prob_funct(Qi,R,Lamda_Knot, Ly, Lamda,X,Y):
    
    if en_debug==1: #handel exception condition #en_debug_0811
        
        tot_Lamda = Lamda #already including Lamda_Knot and Ly
        
        
        if tot_Lamda==0:
            Lamda_Knot = 1
            Ly = 0
            Lamda = 1
            tot_Lamda = 1
        
        if (X==Y)and(Y==R):
            
            Q = ((Qi*Lamda_Knot) + Ly)/tot_Lamda
            
        elif (Y!=X)and(Y!=R):
            
            Q = (((1-Qi)*Lamda_Knot/3) + Ly)/tot_Lamda  
        
        else:
            
            Q = ((Qi*Lamda_Knot/2) + Ly)/tot_Lamda
        
    else: #old code
    
        if (X==Y)and(Y==R):
            
            Q = ((Qi*Lamda_Knot) + j)/(Lamda+Lamda_Knot)
            
        elif (Y!=X)and(Y!=R):
            
            Q = (((1-Qi)*Lamda_Knot/3) + j)/(Lamda+Lamda_Knot)    
        
        else:
            
            Q = ((Qi*Lamda_Knot/2) + j)/(Lamda+Lamda_Knot)
 
        
    return math.log(Q)
        
def process_count_line(row):

    #pdb.set_trace() 

    #j.append(row[1]) # ref position
    ref_pos = row[1]
    #R = row[2].upper() # ref base
    #pdb.set_trace()
    R = row[2] #.upper() hope to keep snp res as original (lower case: intron?)
    Ru = R.upper()
    
    Lamda_Knot = float(row[3]) # expression level
    [Na,Nc,Ng,Nt] = map(int, row[4:8]) # counts
    a = row[8:]
    
    if 0:
        b=[]                    
        for i in range(0,len(a)):
            if i==0:
                b=a[i].split(",")
            else:
                c=a[i].split(",")
                b=b+c
                   
        Lamda = []
        mi=[]
        Y=[]    
        Qi = []  
       
        for m in range(0,len(b),4):
            Y.append(b[m])
            
        for m in range(1,len(b),4):
            if b[m]=='I':
                Qi.append(0.999)
                
        for m in range(2,len(b),4):
            mi.append(int(b[m]))
            
        for m in range(3,len(b),4):
            Lamda.append(int(b[m]))
        
    else:
        #speed-up trial
        #pdb.set_trace()
        Lamda = []
        Ly = []
        Y = []
        Qi = []
        
        for i in range(0, len(a)):
            b=a[i].split(",")
            Y.append(b[0])
            Qi.append(0.999) #convert b[1] to score in future
            #mi.append(int(b[2]))
            Ly.append(float(b[2])) #L_y
            #Lamda.append(int(b[3]))
            Lamda.append(float(b[3]))#L_sum
        
    prob_A = 0
    prob_C = 0
    prob_G = 0
    prob_T = 0
    
    #pdb.set_trace()
    
    prob_A = (JP(Qi,Na,Nc,Ng,Nt,Ru,Lamda_Knot, Ly, Lamda,'A',Y))

    prob_C = (JP(Qi,Na,Nc,Ng,Nt,Ru,Lamda_Knot, Ly, Lamda,'C',Y))

    prob_G = (JP(Qi,Na,Nc,Ng,Nt,Ru,Lamda_Knot, Ly, Lamda,'G',Y))

    prob_T = (JP(Qi,Na,Nc,Ng,Nt,Ru,Lamda_Knot, Ly, Lamda,'T',Y))

    return [ref_pos, R, prob_A, prob_C, prob_G, prob_T] #=process_count_line(row)
#----------------------------------------
#if 1:     

def write_row_res(row_res, out_file): #row_res = [ref_pos, ref_base, prob_A, prob_C, prob_G, prob_T]
    try:
        #caller_output
        out_file.write(str(row_res[0])+'\t')
        out_file.write(str(row_res[1])+'\t') #('ref_base:'+str(ref_base)+'\t')
        out_file.write('%+10.15f\t'%row_res[2]) #('t_base_A:'+str(prob_A)+"\t")
        out_file.write('%+10.15f\t'%row_res[3])
        out_file.write('%+10.15f\t'%row_res[4])
        out_file.write('%+10.15f\t'%row_res[5])
        out_file.write("\n")
    except:
        pdb.set_trace()
        print('caller output error')

def snp_found(row_res, out_file, alpha=0):#row_res = [ref_pos, ref_base, prob_A, prob_C, prob_G, prob_T]
    dProb = None # Prob(target)-Prob(ref), to check statistics for ROC curve purpose
    try:
        #special case - D:
        prob_A=row_res[2]
        prob_C=row_res[3]
        prob_G=row_res[4]
        prob_T=row_res[5]
        
        if prob_A==0 and prob_C==0 and prob_G==0 and prob_T==0:
            return
            
        #pdb.set_trace()
        rB = row_res[1] #.upper()
        rBu = rB.upper()
        ACGT = 'ACGT'
        #debug 9/12
        tBs = []
        vals = []
        ref_Prob = None
        for i in range(4):
            if ACGT[i] != rBu:
                tBs.append(ACGT[i])
                vals.append(row_res[2+i])
            else:
                ref_Prob = row_res[2+i]
        
        #pdb.set_trace()
        #vals = row_res[2::]
        possible_SNP_Prob = max(vals)
        possible_SNP_idx = vals.index(max(vals))
        possible_SNP_tB = tBs[possible_SNP_idx]
        #if ref_Prob != None:
        #    dProb = possible_SNP_Prob - ref_Prob
            
        if ref_Prob != None and possible_SNP_Prob > ref_Prob + alpha:
            dProb = possible_SNP_Prob - ref_Prob
            #if dProb > 3000:
            #    pdb.set_trace()
            out_file.write(str(row_res[0])+'\t')
            out_file.write(rB + '\t-->\t')
            out_file.write(possible_SNP_tB + '\t\n')        
    except:
        pdb.set_trace()
        print('snp found error')
    
    return [dProb, row_res[0], rB, possible_SNP_tB] #may add other results for future usage

def final_caller(input_fn_count_file = '/count_sample.txt', 
                 output_fn = '/caller_output.txt', 
                 exception_output_fn='/caller_output_exception.txt',
                 found_snp_fn = '/caller_output_snp_found.txt',
                 alpha = 0):
    
    if en_debug == 0:
        Count_file = '/home/kmazooji/jatindeep/Bioinformatics/Simulation/data/tophat_out_m/output/count.txt'
        ofile = open("/home/kmazooji/jatindeep/Bioinformatics/Simulation/data/tophat_out_m/output/caller_output.txt", "w")
    else:
        #Count_file = '/home/olivo/Desktop/SNP-Calling-Summer15 (1)/data/count_sample.txt' #local linux
        Count_file = Default_Ref_Path + input_fn_count_file #remote linux
        ofile = open(Default_Ref_Path+output_fn, 'w')       
        ofile_exception = open(Default_Ref_Path + exception_output_fn, 'w') 
        ofile_snp_found = open(Default_Ref_Path + found_snp_fn, 'w')
        ofile_dProb_dmp = open(Default_Ref_Path + output_fn + '.dProb_dmp', 'w')               

    snps_res = {}
    
    #pdb.set_trace()
    num_lines = sum(1 for line in open(Count_file))
    bar = Bar('Processing', max=num_lines)
    #pdb.set_trace()
    
    with open(Count_file) as counts:
        
        reader = csv.reader(counts, delimiter='\t')
        
        for row in reader:
            bar.next()
            if len(row) >= 8:
                row_res = []
                try: 
                    #pdb.set_trace() 
                    #if int(row[0])==940:
                    #    pdb.set_trace()
                    row_res =process_count_line(row) #[ref_pos, ref_base, prob_A, prob_C, prob_G, prob_T]
                except:
                    pdb.set_trace()
                    row_str = ' '.join(row)
                    print('\nprocess_count_line exception with row=%s' % row_str)
                    ofile_exception.write('%s\n'%row_str)
                    #pdb.set_trace()
                    #row_res =process_count_line(row)
                    #further debug
                    #[ref_pos, ref_base, prob_A, prob_C, prob_G, prob_T] =process_count_line(row)
                                    
                if len(row_res)>0:
                    #pdb.set_trace()
                    write_row_res(row_res, ofile)
                    snp_found_itm = snp_found(row_res, ofile_snp_found, alpha)
                    #pdb.set_trace()
                    dProb = snp_found_itm[0]
                    if dProb != None: #snp found condition
                        ofile_dProb_dmp.write("%f\n"%dProb)
                        res_itm = snps_res.setdefault(int(snp_found_itm[1]), []) #append([snp_found_itm[1], snp_found_itm[2], snp_found_itm[3], dProb])
                        res_itm.append([snp_found_itm[2], snp_found_itm[3], dProb])
                        #pos, rB, tB, dProb
                        #pdb.set_trace()
                    
        bar.finish()
                
    #Count_file.close()
    ofile.close()
    ofile_exception.close()
    ofile_snp_found.close()
    ofile_dProb_dmp.close()
    
    return snps_res
    
def gen_snps_res_from_caller_output(caller_output_addr): #snps_res = 
    snps_res = {}
    ACGT = 'ACGT'
    
    with open(caller_output_addr) as co:
        
        reader = csv.reader(co, delimiter='\t')
        
        for row in reader:
            pos = int(row[0])
            rBu = row[1].upper()
            
            ref_Prob = None
            tBs = []
            vals = []
            for i in range(4):
                if ACGT[i] != rBu:
                    tBs.append(ACGT[i])
                    vals.append(float(row[2+i]))
                else:
                    ref_Prob = float(row[2+i])
                    
            possible_SNP_Prob = max(vals)
            possible_SNP_idx = vals.index(max(vals))
            possible_SNP_tB = tBs[possible_SNP_idx]
            
            if ref_Prob != None and possible_SNP_Prob > ref_Prob:
                #pdb.set_trace()
                dProb = possible_SNP_Prob - ref_Prob
                res_itm = snps_res.setdefault(pos, []) #append([snp_found_itm[1], snp_found_itm[2], snp_found_itm[3], dProb])
                res_itm.append([rBu, possible_SNP_tB, dProb])
            #else:
            #    print('gen_snps_res_from_caller_output: unexpected')
            #    pdb.set_trace()                
    
    return snps_res

def generate_ROC_statistics():
    
    #pdb.set_trace()
    
    Default_Ref_Path = '/home/sreeramkannan/singleCell/SNP_Calling_Summer15/data_0814/'
    
    #count_dir = '/tophat_out/'
    #count_dir = '/data_GATK/1pass/'
    count_dir = '/data_GATK/2pass/'
        
    #count_fn = "/count_alt_mapping_debug_0911.txt" 
    #count_fn = '/count_altMap_cross_check_star1pass.txt'
    count_fn = '/count_altMap_cross_check_split.txt'
    
    count_rel_address = count_dir + count_fn    
    #caller_output_res_fn = '/caller_output_roc.txt'
    #caller_output_exception_fn = '/caller_output_exception_roc.txt'
    #caller_output_snp_found_fn = '/caller_output_snp_found_roc.txt'
    #caller_output_res_fn = '/caller_output_altMap_cross_check_star1pass.txt'
    #caller_output_exception_fn = '/caller_output_exception_altMap_cross_check_star1pass.txt'
    #caller_output_snp_found_fn = '/caller_output_snp_found_altMap_cross_check_star1pass.txt'
    caller_output_res_fn = '/caller_output_altMap_cross_check_split.txt'
    caller_output_exception_fn = '/caller_output_exception_altMap_cross_check_split.txt'
    caller_output_snp_found_fn = '/caller_output_snp_found_altMap_cross_check_split.txt'
    
    use_caller_output = True
    
    print('generate_ROC_statistics --> get snps_res')       
    pdb.set_trace()
    snps_res = []
    if use_caller_output == False:
        alpha = 0
        snps_res = final_caller(count_rel_address, 
                                caller_output_res_fn, 
                                caller_output_exception_fn,
                                caller_output_snp_found_fn,
                                alpha)
    else:
        snps_res = gen_snps_res_from_caller_output(Default_Ref_Path + caller_output_res_fn)
    
    print('generate_ROC_statistics --> get md0 and fp0 etc statistics')             
    pdb.set_trace()
    
    snp_res_address = Default_Ref_Path + caller_output_snp_found_fn    
    snp_m_address = Default_Ref_Path + '/SNP_m.txt'
    snp_p_address = Default_Ref_Path + '/SNP_p.txt'
    
    snps_list_sorted = get_snp_res_statistics(snp_res_address,
                                              snp_m_address,
                                              snp_p_address)
                                              
    [m_snps_cd, m_snps_md, p_snps_cd, p_snps_md, r_snps_fp] = group_snps(snps_list_sorted)
    
    print('generate_ROC_statistics --> show md0 and fp0 etc statistics')
    pdb.set_trace()
    print('# of mis-detection (m & p):%d'%(len(m_snps_md)+len(p_snps_md)))
    print('...# of mis-detection (m):%d'%len(m_snps_md))
    print('...# of mis-detection (p):%d'%len(p_snps_md))
    
    print('# of false-positive:%d'%len(r_snps_fp))
    
    print('# of correct-detection (m & p):%d'%(len(m_snps_cd)+len(p_snps_cd)))
    print('...# of correct-detection (m):%d'%len(m_snps_cd))
    print('...# of correct-detection (p):%d'%len(p_snps_cd))
    
    print('generate_ROC_statistics --> generate md and fp statistics')
    pdb.set_trace()
    
    alpha0 = 0
    alpha_list = np.arange(2, 200, 2) #1st 0
    md_list = [] # num of mis-detections, wrt alpha
    md0 = (len(m_snps_md)+len(p_snps_md))
    fp_list = [] # num of false-positives, wrt alpha
    fp0 = (len(r_snps_fp))
    
    for i in range(len(alpha_list)):
        md = md0
        fp = fp0
        for itm in m_snps_cd:
            dProb = snps_res.get(itm[0])[0][2]
            if dProb < alpha_list[i]: #previously detected item will not be detected
                #pdb.set_trace()
                md = md + 1
        for itm in p_snps_cd:
            dProb = snps_res.get(itm[0])[0][2]
            if dProb < alpha_list[i]:
                #pdb.set_trace()
                md = md + 1
        for itm in r_snps_fp:
            dProb = snps_res.get(itm[0])[0][2]
            if dProb < alpha_list[i]: #previously falsely detected item will not be detected
                #pdb.set_trace()
                fp = fp - 1
        md_list.append(md)
        fp_list.append(fp)
    
    alpha_list = [alpha0] + [a for a in alpha_list]
    md_list = [md0] + md_list
    fp_list = [fp0] + fp_list
    
    print('generate_ROC_statistics --> dump md and fp statistics')
    pdb.set_trace()
    
    dmp_md_fp_addr = Default_Ref_Path + caller_output_res_fn + '.dmp_md_fp'   
    dmp_md_fp = open(dmp_md_fp_addr, 'w+')
    for i in range(len(alpha_list)):
        dmp_md_fp.write('%.2f\t'%float(alpha_list[i]))
        dmp_md_fp.write('%d\t'%md_list[i])
        dmp_md_fp.write('%d\n'%fp_list[i])
    dmp_md_fp.close()
    
    print('generate_ROC_statistics --> exit')
    pdb.set_trace()
    #ready to draw ROC curve
    
    
    
    

# main
#final_caller()
if __name__ == "__main__":
    
    pdb.set_trace()
    
    generate_ROC_statistics()
    
    """ ==========
    # check count (alt mapping)
    Default_Ref_Path = '/home/sreeramkannan/singleCell/SNP_Calling_Summer15/data_0814/'
    #COV_fn = '/coverage_qt90.txt'
    #COV_fn = '/coverage.txt'
    #coverage_address = Default_Ref_Path + COV_fn
    
    tophat_dir = '/tophat_out/'
    
    #sam_address = Default_Ref_Path + tophat_dir + '/accepted_hits.sam' #for test purpose
    #ref_address = Default_Ref_Path + '/Chr15.fa'
    
    #count_fn = '/count_filtCovQt90.txt'
    #count_fn = '/count_allCov.txt'
    #count_fn = "/count_alt_mapping_debug.txt"
    #count_fn = "/count_alt_mapping_debug_90346907.txt"
    count_fn = "/count_alt_mapping_debug_0911.txt"
    
    #generate_count_file(sam_address, ref_address, coverage_address, count_fn) #en_debug_0817
    
    count_rel_address = tophat_dir + count_fn #for test purpose
    ## ---------- caller
    
    caller_output_res_fn = '/caller_output_roc.txt'
    caller_output_exception_fn = '/caller_output_exception_roc.txt'
    caller_output_snp_found_fn = '/caller_output_snp_found_roc.txt'
    
    final_caller(tophat_dir + count_fn, 
                 caller_output_res_fn, 
                 caller_output_exception_fn,
                 caller_output_snp_found_fn)    
    """
    
    """ ==========
    # caller using parallel computing
    
    print('final caller: starts')

    count_fn_rel_dir = '/tophat_out/count_split/'
    if len(sys.argv)>=2:
        count_fn_rel_dir = '/tophat_out/count_split/'
        count_fn = count_fn_rel_dir + sys.argv[1]
    else:
        count_fn = '/tophat_out/count.txt'

    if len(sys.argv)>=3:
        caller_op_fn = sys.argv[2]
    else:
        caller_op_fn = '/caller_output.txt'
    
    if len(sys.argv)>=4:
        caller_op_exception_fn = sys.argv[3]
    else:
        caller_op_exception_fn = '/caller_output_exception.txt'

    if len(sys.argv)>=5:
        caller_op_snp_fn = sys.argv[4]
    else:
        caller_op_snp_fn = '/caller_output_snp_found.txt'

    #if len(sys.argv)>=6:
    #    caller_log = 
 
    print('final caller: starts: %s\nRef data folder: %s'%(count_fn,Default_Ref_Path))
    final_caller(count_fn,
                 caller_op_fn,
                 caller_op_exception_fn,
                 caller_op_snp_fn)

    print('final caller: ready to exit: %s'%count_fn)
    #pdb.set_trace()
    """

