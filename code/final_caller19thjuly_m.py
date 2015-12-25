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
        
def process_count_line(row, Threshold_num_reads=1):

    #pdb.set_trace()
    
    if sum(map(int, row[4:8]))<Threshold_num_reads:
        return [] # too small # of reads, skip snp calling

    #j.append(row[1]) # ref position
    #pdb.set_trace()
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

def snp_found(row_res, out_file):#row_res = [ref_pos, ref_base, prob_A, prob_C, prob_G, prob_T]
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
        vals = row_res[2::]
        tB_idx = vals.index(max(vals))
        tB = ACGT[tB_idx]
        if rBu != tB:
            out_file.write(str(row_res[0])+'\t')
            out_file.write(rB + '\t-->\t')
            out_file.write(tB + '\t\n')        
    except:
        pdb.set_trace()
        print('snp found error')

def final_caller(input_fn_count_file = '/count_sample.txt', 
                 output_fn = '/caller_output.txt', 
                 exception_output_fn='/caller_output_exception.txt',
                 found_snp_fn = '/caller_output_snp_found.txt',
                 Threshold_num_reads = 1):
    
    if en_debug == 0:
        Count_file = '/home/kmazooji/jatindeep/Bioinformatics/Simulation/data/tophat_out_m/output/count.txt'
        ofile = open("/home/kmazooji/jatindeep/Bioinformatics/Simulation/data/tophat_out_m/output/caller_output.txt", "w")
    else:
        #Count_file = '/home/olivo/Desktop/SNP-Calling-Summer15 (1)/data/count_sample.txt' #local linux
        if Default_Ref_Path not in input_fn_count_file:        
            Count_file = Default_Ref_Path + input_fn_count_file #rel addr ; remote linux
        else:
            Count_file = input_fn_count_file #abs addr
        
        if Default_Ref_Path not in output_fn:
            ofile = open(Default_Ref_Path+output_fn, 'w')
        else:
            ofile = open(output_fn, 'w')
        
        if Default_Ref_Path not in exception_output_fn:
            ofile_exception = open(Default_Ref_Path + exception_output_fn, 'w')
        else:
            ofile_exception = open(exception_output_fn, 'w') 
        
        if Default_Ref_Path not in found_snp_fn:
            ofile_snp_found = open(Default_Ref_Path + found_snp_fn, 'w')                 
        else:               
            ofile_snp_found = open(found_snp_fn, 'w')
        
    #pdb.set_trace()
    num_lines = sum(1 for line in open(Count_file))
    bar = Bar('Processing', max=num_lines)
    #pdb.set_trace()
    
    with open(Count_file) as counts:
        
        reader = csv.reader(counts, delimiter='\t')
        #reader = csv.reader(counts) #, delimiter='\t')
        
        for row in reader:
            if len(row)==1:
                row = row[0].split()
            bar.next()
            if float(row[3])==0:
                continue
            elif len(row) >= 8:
                row_res = []
                try: 
                    #pdb.set_trace() 
                    #if int(row[0])==940:
                    #    pdb.set_trace()
                    row_res =process_count_line(row, Threshold_num_reads) #[ref_pos, ref_base, prob_A, prob_C, prob_G, prob_T]
                except:
                    #pdb.set_trace()
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
                    snp_found(row_res, ofile_snp_found)
        bar.finish()
                
    #Count_file.close()
    ofile.close()
    ofile_exception.close()
    ofile_snp_found.close()

def merge_caller_res(Default_Ref_Path,
              caller_op_x,  caller_op_exception_x,  caller_op_snp_x, num_p, #inputs
              caller_op_fn, caller_op_exception_fn, caller_op_snp_fn): #outputs
    
    if Default_Ref_Path in caller_op_fn: #abs address
        caller_op_file = open(caller_op_fn, 'w')
    else:
        caller_op_file = open(Default_Ref_Path+caller_op_fn, 'w')
    for i in range(num_p):
        caller_op_x_fn = '%s%02d.txt'%(caller_op_x, i)
        if Default_Ref_Path not in caller_op_x_fn:
            caller_op_x_fn = Default_Ref_Path + caller_op_x_fn
        #print('merge: %s'%caller_op_x_fn)
        with open(caller_op_x_fn) as op:            
            reader = csv.reader(op)            
            for row in reader:
                caller_op_file.write('%s\n'%row[0])
        #subprocess.call('rm '+Default_Ref_Path + caller_op_x_fn, shell = True)    
    caller_op_file.close()
    
    if Default_Ref_Path in caller_op_exception_fn:
        caller_op_exception_file = open(caller_op_exception_fn, 'w')
    else:
        caller_op_exception_file = open(Default_Ref_Path+caller_op_exception_fn, 'w')
    for i in range(num_p):
        caller_op_exception_x_fn = '%s%02d.txt'%(caller_op_exception_x, i)
        if Default_Ref_Path not in caller_op_exception_x_fn:
            caller_op_exception_x_fn = Default_Ref_Path + caller_op_exception_x_fn    
        #print('merge: %s'%caller_op_exception_x_fn)
        with open(caller_op_exception_x_fn) as op:            
            reader = csv.reader(op)            
            for row in reader:
                caller_op_exception_file.write('%s\n'%row[0])
        #subprocess.call('rm '+Default_Ref_Path + caller_op_exception_x_fn, shell = True)    
    caller_op_exception_file.close()
    
    if Default_Ref_Path in caller_op_snp_fn:
        caller_op_snp_file = open(caller_op_snp_fn, 'w')
    else:
        caller_op_snp_file = open(Default_Ref_Path+caller_op_snp_fn, 'w')
    for i in range(num_p):
        caller_op_snp_x_fn = '%s%02d.txt'%(caller_op_snp_x, i)
        if Default_Ref_Path not in caller_op_snp_x_fn:
            caller_op_snp_x_fn = Default_Ref_Path + caller_op_snp_x_fn
        print('merge: %s'%caller_op_snp_x_fn)
        with open(caller_op_snp_x_fn) as op:            
            reader = csv.reader(op)            
            for row in reader:
                #pdb.set_trace()
                caller_op_snp_file.write('%s\n'%row[0])
        #subprocess.call('rm '+Default_Ref_Path + caller_op_snp_x_fn, shell = True)    
    caller_op_snp_file.close()

# main
#final_caller()
if __name__ == "__main__":
    
    if len(sys.argv)<2: #non-parallel
        pdb.set_trace()
        
        # check count (alt mapping)
        #Default_Ref_Path = '/home/sreeramkannan/singleCell/SNP_Calling_Summer15/data_0814/'
        #Default_Ref_Path = '/home/sreeramkannan/singleCell/SNP_Calling_Summer15/data_SNP1k_Reads10M/'
        Default_Ref_Path = '/home/sreeramkannan/singleCell/SNP_Calling_Summer15/data_SNP20_Reads100K/'
        #Default_Ref_Path = '/home/olivo/Downloads/0907_count_alt_mapping/bkp/data_0814/'
        #COV_fn = '/coverage_qt90.txt'
        #COV_fn = '/coverage.txt'
        #coverage_address = Default_Ref_Path + COV_fn
        
        tophat_dir = '' #'tophat_out/count_y/' #'/tophat_out/' #"data_GATK/2pass/" 
        
        #sam_address = Default_Ref_Path + tophat_dir + '/accepted_hits.sam' #for test purpose
        #ref_address = Default_Ref_Path + '/Chr15.fa'
        
        #count_fn = '/count_filtCovQt90.txt'
        #count_fn = '/count_allCov.txt'
        #count_fn = "/count_alt_mapping_debug.txt"
        #count_fn = "/count_alt_mapping_debug_90346907.txt"
        #count_fn = "/count_alt_mapping_debug_0911.txt"
        #count_fn = "/count_altMap_cross_check_dedupped_debug_case2.txt"
        #count_fn = "/count_dirAltMap_cross_check_dedupped.txt"
        #count_fn = '/count_dirAltMap.txt'
        #count_fn = '/count_y03.txt'
        count_fn = 'caller_output_exception_accepted_hits_T1.txt' #exception analysis
        #generate_count_file(sam_address, ref_address, coverage_address, count_fn) #en_debug_0817
        
        #count_rel_address = tophat_dir + count_fn #for test purpose
        ## ---------- caller
        #final_caller(count_rel_address, '/caller_output.txt')
        
        
        
        final_caller(tophat_dir + count_fn,
                     '/caller_op_tmp.txt', #_cross_check_dedupped.txt',
                     '/caller_op_exception_tmp.txt', #_cross_check_dedupped.txt',
                     '/caller_op_snp_tmp.txt') #_cross_check_dedupped.txt')
        
    
    else:
        
        # caller using parallel computing
        
        print('final caller: starts')
    
        count_fn_rel_dir = sys.argv[1] #'/tophat_out/count_split/'; can also be abs dir now
        count_fn = count_fn_rel_dir + sys.argv[2]
        
        #pdb.set_trace()
        
        caller_op_fn = sys.argv[3]
        
        caller_op_exception_fn = sys.argv[4]
        
        caller_op_snp_fn = sys.argv[5]
        
        Threshold_num_reads = int(sys.argv[6])
        
        print('starts: %s'%(count_fn))
        
        final_caller(count_fn,
                     caller_op_fn,
                     caller_op_exception_fn,
                     caller_op_snp_fn,
                     Threshold_num_reads)
    
        print('exit: %s'%count_fn)
        #pdb.set_trace()
        

