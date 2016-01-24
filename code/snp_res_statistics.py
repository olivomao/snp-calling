import math
import csv
import operator
import pickle
from Address import *
import pdb
from debug_MACRO import *
from progress.bar import Bar
from operator import itemgetter, attrgetter
from itertools import izip
from scipy.stats import poisson

"""
16/1/24: 
- snp filtering based on modified expression levels of alt mapping
  (the alt mappings of the reads at the loci are grouped based on the direction and base of the reads)
"""

def get_snps(snp_dic, snp_addr, tag=''):
    
    """
    description
    """
    print('... build snps dic from: %s'%snp_addr)    
    
    with open(snp_addr) as snp_file:        
        reader = csv.reader(snp_file, delimiter='\t')

        i = 0        
        for row in reader:
            i = i+1
            
            rPos = int(row[0])
            rB = row[1].split()[0].upper()
            nB = row[3].split()[0].upper()
            #print('...%dth snp row: %s'%(i, str(row)))
            
            if rPos in snp_dic.keys():
                #print('...\tadd to existing key = %d'%rPos)
                
                #pdb.set_trace()
                snp_dic[rPos].append([rB, nB, tag])
                #print('...\texisting key\'s val: %s'%str(snp_dic[rPos]))
                #print(snp_dic[rPos])
                #print('')
            else:
                #print('...\tadd new key = %d'%rPos)                
                #pdb.set_trace()
                snp_dic[rPos] = [[rB, nB, tag]]
                
                #print('')

def get_snp_res_statistics(snp_res_address,
                           snp_m_address,
                           snp_p_address):
                               
    """
    description
    """
    print('get_snp_res_statistics:')
    print('snp_found_addr:%s'%snp_res_address)
    print('snp_m_addr:%s'%snp_m_address)
    print('snp_p_addr:%s'%snp_p_address) 
    print('')
    #pdb.set_trace()
    
    snps = {}
    get_snps(snps, snp_m_address, 'm')
    get_snps(snps, snp_p_address, 'p')
    get_snps(snps, snp_res_address, 'r')   
    
    snps_list = snps.items()
    snps_list_sorted = sorted(snps_list, key=itemgetter(0) )    
    
    return snps_list_sorted

def group_snps(snps_list_sorted):
    
    """
    description
    """
    print('\ngroups snps into different types for further analysis')
    #pdb.set_trace()
    
    m_snps_cd = [] # snps correctly detected (m)
    m_snps_md = [] # snps miss-detected
    p_snps_cd = [] # snps correctly detected (p)
    p_snps_md = [] # snps miss-detected
    r_snps_fp = [] # snps in result but are false positive (not in true m or p SNPs)

    for i in range(len(snps_list_sorted)):
        #print('process %d/%d-th snp node'%(i+1, len(snps_list_sorted)))
        
        snp_node = snps_list_sorted[i]        
        snp_Pos = snp_node[0]
        snp_Info = snp_node[1]
        
        #print('...pos=%d'%snp_Pos)
        #print('...info=%s'%str(snp_Info))
        
        if(len(snp_Info)==1): #false positive or mis-detection
            if(snp_Info[0][2]=='m'):
                #print('...m snp mis-detected')
                m_snps_md.append(snp_node)
                #pdb.set_trace()
            elif(snp_Info[0][2]=='p'):
                #print('...p snp mis-detected')
                p_snps_md.append(snp_node)
                #pdb.set_trace()
            elif(snp_Info[0][2]=='r'):
                #print('...res snp false positive')
                r_snps_fp.append(snp_node)  
                #pdb.set_trace()
            else:
                #pdb.set_trace()
                print('exception1')
                
        else:
            snp_node_tags = [itm[2] for itm in snp_Info]
            if 'r' in snp_node_tags and 'm' in snp_node_tags and len(snp_node_tags)==2:
                #print('...m snp correctly detected')
                m_snps_cd.append(snp_node)
                #pdb.set_trace()
            elif 'r' in snp_node_tags and 'p' in snp_node_tags and len(snp_node_tags)==2:
                #print('...p snp correctly detected')
                p_snps_cd.append(snp_node)
                #pdb.set_trace()    
            else:
                #pdb.set_trace()
                print('exception2')
                
        
        #pdb.set_trace()
    
    return [m_snps_cd, m_snps_md, p_snps_cd, p_snps_md, r_snps_fp]
    
def dmp_snp_res(snps_list, res_file):
    for i in range(len(snps_list)):
        itm = snps_list[i]
        pos = itm[0]
        info = itm[1]
        res_file.write('%d\n'%pos)
        for j in range(len(info)):
            res_file.write('%s\n'%str(info[j]))
        res_file.write('\n')
    res_file.write('\n')
    return

def isVCF(snp_res_address):
    #pdb.set_trace()
    
    res = 0
    if '.vcf' in snp_res_address:
        res = 1    
    return res
    
def extractVCF(snp_res_address): # x.vcf --> x.vcf.txt
    #pdb.set_trace()
    
    res_address = snp_res_address + '.txt'
    res_file = open(res_address, 'w+')
    res_address2 = snp_res_address + '.2.txt' #indel and deletion
    res_file2 = open(res_address2, 'w+')
    
    """
    description
    """
    #pdb.set_trace()
    print('extractVCF from: %s to %s'%(snp_res_address, res_address))    
    
    with open(snp_res_address) as snp_file:        
        reader = csv.reader(snp_file, delimiter='\t')

        for row in reader:
            #pdb.set_trace()
            if row[0][0] != '#':
                #pdb.set_trace()
                rPos = int(row[1])-1 # convert 1-base pos to 0-base pos
                rB = row[3].upper()
                nB = row[4].upper()
                if len(rB)==1 and len(nB)==1:
                    res_file.write('%d\t%s\t-->\t%s\n'%(rPos, rB, nB))
                else:
                    #pdb.set_trace()
                    res_file2.write('%d\t%s\t-->\t%s\n'%(rPos, rB, nB)) 
    
    res_file.close()
    res_file2.close()
    
    return res_address

def do_snp_res_statistics(snp_res_address, snp_m_address, snp_p_address, snp_res_stat_fn):

    if isVCF(snp_res_address):
        snp_res_address = extractVCF(snp_res_address) # x.vcf --> x.vcf.txt
    
    #pdb.set_trace()
    snps_list_sorted = get_snp_res_statistics(snp_res_address,
                                              snp_m_address,
                                              snp_p_address)
    
    #pdb.set_trace()
    [m_snps_cd, m_snps_md, p_snps_cd, p_snps_md, r_snps_fp] = group_snps(snps_list_sorted)
    
    print('# of mis-detection (m & p):%d'%(len(m_snps_md)+len(p_snps_md)))
    print('...# of mis-detection (m):%d'%len(m_snps_md))
    print('...# of mis-detection (p):%d'%len(p_snps_md))
    
    print('# of false-positive:%d'%len(r_snps_fp))
    
    print('# of correct-detection (m & p):%d'%(len(m_snps_cd)+len(p_snps_cd)))
    print('...# of correct-detection (m):%d'%len(m_snps_cd))
    print('...# of correct-detection (p):%d'%len(p_snps_cd))
    
    
    #pdb.set_trace()
    if Default_Ref_Path in snp_res_stat_fn:
        res_file = open(snp_res_stat_fn, 'w')
    else:
        res_file = open(Default_Ref_Path + snp_res_stat_fn, 'w')
    
    res_file.write('#========== Input Files:\n')
    res_file.write('' + snp_res_address+'\n')
    res_file.write('' + snp_m_address+'\n')
    res_file.write('' + snp_p_address+'\n\n')
    
    res_file.write('#========== Results Summary:\n')
    res_file.write('num of mis-detection (m & p):%d\n'%(len(m_snps_md)+len(p_snps_md)))
    res_file.write('...num of mis-detection (m):%d\n'%len(m_snps_md))
    res_file.write('...num of mis-detection (p):%d\n'%len(p_snps_md))    
    res_file.write('num of false-positive:%d\n'%len(r_snps_fp))    
    res_file.write('num of correct-detection (m & p):%d\n'%(len(m_snps_cd)+len(p_snps_cd)))
    res_file.write('...num of correct-detection (m):%d\n'%len(m_snps_cd))
    res_file.write('...num of correct-detection (p):%d\n\n'%len(p_snps_cd))
    
    res_file.write('#========== Results (m SNP - correct detection):\n')
    dmp_snp_res(m_snps_cd, res_file)
    
    res_file.write('#========== Results (m SNP - miss detection):\n')
    dmp_snp_res(m_snps_md, res_file)
    
    res_file.write('#========== Results (p SNP - correct detection):\n')
    dmp_snp_res(p_snps_cd, res_file)
    
    res_file.write('#========== Results (p SNP - miss detection):\n')
    dmp_snp_res(p_snps_md, res_file)
    
    res_file.write('#========== Results (res SNP - false positive):\n')
    dmp_snp_res(r_snps_fp, res_file)        
    
    res_file.close()
    
    return [m_snps_cd, m_snps_md, p_snps_cd, p_snps_md, r_snps_fp]

def filter_snp_lam_half_filt2(snp_res_address,
                              filt_snp_res_address,
                              count_altInfo_address,
                              count_abs_address): # check shadow snp from counts_altInfo.txt
    
    bases = {"A": 0, "C": 1, "G": 2, "T": 3}
    
    print('lambda half filtering 2 (alt count modified)')
    #pdb.set_trace()
    
    found_snps = {}
    get_snps(found_snps, snp_res_address)
    
    num_removed_snps = 0
    print('== filter_snp == ')
    
    #pdb.set_trace()
    #with open(count_altInfo_address) as count_altInfo_file:
    with open(count_altInfo_address) as count_altInfo_file, open(count_abs_address) as count_file:
        #for line in count_altInfo_file:
        for xline, yline in izip(count_altInfo_file, count_file):
            
            #pdb.set_trace()
            
            x = xline.split()
            y = yline.split()

            gp = int(y[1])
            
            if gp not in found_snps:
                continue
            
            #gp considered as a snp
            found_snp_base = found_snps[gp][0][1].upper()           
            N_snp_base = y[4+bases[found_snp_base]]
            N_snp_base = int(N_snp_base)
            
            #detected snp has alt mappings
            if len(x)>=3:
                #pdb.set_trace()
                
                positions = []
                lambdas = []
                poisson_pmf = []
                
                #add first itm
                itm = x[1].split(',')
                q = int(itm[0])
                q_e = float(itm[1])
                pmf = poisson(q_e/2.0).pmf(N_snp_base)
                
                positions.append(q)
                lambdas.append(q_e)
                poisson_pmf.append(pmf)
                
                #add alt mappings - find (pos,l) pairs
                pos_el_dic = {}
                itms = x[2].split('[')
                for itm in itms:
                    if len(itm)>0 and itm[1]==found_snp_base:
                        #pdb.set_trace()
                        pos_el = itm[3:len(itm)-2].split(',')
                        for i in range(len(pos_el)/2):
                            q = int(pos_el[i*2])
                            q_e = float(pos_el[i*2+1])
                            if q not in pos_el_dic: # +A and -A may have same alt mappings
                                pos_el_dic[q]=q_e
                                
                #add alt mappings - calculate (pos, l) pairs
                for q, q_e in pos_el_dic.iteritems():
                    pmf = poisson(q_e/2.0).pmf(N_snp_base)
                    positions.append(q)
                    lambdas.append(q_e)
                    poisson_pmf.append(pmf)
                
                #filter snp
                if poisson_pmf.index(max(poisson_pmf)) != 0:
                        #pdb.set_trace()
                        del found_snps[gp]
            
    #pdb.set_trace()
    filt_snp_res_file = open(filt_snp_res_address, 'w')
    
    for itm in found_snps.items():
        #pdb.set_trace()
        gp = itm[0]
        rB = itm[1][0][0]
        tB = itm[1][0][1]
        filt_snp_res_file.write(str(gp)+'\t')
        filt_snp_res_file.write(rB + '\t-->\t')
        filt_snp_res_file.write(tB + '\t\n')
    
    filt_snp_res_file.close()
    
    return

def filter_snp_lam_half_filt(snp_res_address,
                             filt_snp_res_address,
                             count_altInfo_address,
                             count_abs_address): # check shadow snp from counts_altInfo.txt
    
    bases = {"A": 0, "C": 1, "G": 2, "T": 3}
    
    print('lambda half filtering')
    #pdb.set_trace()
    
    found_snps = {}
    get_snps(found_snps, snp_res_address)
    
    num_removed_snps = 0
    print('== filter_snp == ')
    
    #pdb.set_trace()
    #with open(count_altInfo_address) as count_altInfo_file:
    with open(count_altInfo_address) as count_altInfo_file, open(count_abs_address) as count_file:
        #for line in count_altInfo_file:
        for xline, yline in izip(count_altInfo_file, count_file):
            
            #pdb.set_trace()
            
            x = xline.split()
            y = yline.split()            
            
            if len(x)>=3:
                #pdb.set_trace()
                
                """
                to modify
                """
                p = int(x[1].split(',')[0])
                #p_e = float(x[1].split(',')[1])
                
                if p in found_snps:
                    #pdb.set_trace()
                    
                    found_snp_base = found_snps[p][0][1]
                    N_snp_base = y[4+bases[found_snp_base]]
                    N_snp_base = int(N_snp_base)
                    
                    positions = []
                    lambdas = []
                    poisson_pmf = []
                    for itm in x[1:]:
                        q = int(itm.split(',')[0])
                        q_e = float(itm.split(',')[1])
                        positions.append(q)
                        lambdas.append(q_e)
                        pmf = poisson(q_e/2.0).pmf(N_snp_base)
                        poisson_pmf.append(pmf)
                    
                    if poisson_pmf.index(max(poisson_pmf)) != 0:
                        #pdb.set_trace()
                        del found_snps[p]
                
                """
                flag_remove = False
                if p in found_snps:
                    
                    for itm in x[2:]:
                        q = int(itm.split(',')[0])
                        q_e = float(itm.split(',')[1])
                        if q in found_snps and q_e>p_e: #a potential shadow snp
                            flag_remove = True
                            num_removed_snps += 1
                            #print('p=%d, p_e=%f, q=%d, q_e=%f (num_removed_snps=%d)'%(p, p_e, q, q_e, num_removed_snps))
                            break
                    #end for itm in x[2:]:
                    
                    if flag_remove==True:
                        del found_snps[p]                                
                        #pdb.set_trace() 
                """
    #end with
    #pdb.set_trace()
    filt_snp_res_file = open(filt_snp_res_address, 'w')
    
    for itm in found_snps.items():
        #pdb.set_trace()
        gp = itm[0]
        rB = itm[1][0][0]
        tB = itm[1][0][1]
        filt_snp_res_file.write(str(gp)+'\t')
        filt_snp_res_file.write(rB + '\t-->\t')
        filt_snp_res_file.write(tB + '\t\n')
    
    filt_snp_res_file.close()
    
    return

def filter_snp(snp_res_address,
               filt_snp_res_address,
               count_altInfo_address,
               count_abs_address): # check shadow snp from counts_altInfo.txt
    
    pdb.set_trace()
    
    found_snps = {}
    get_snps(found_snps, snp_res_address)        
    #pdb.set_trace()
    
    #filt_snps = []
    num_removed_snps = 0
    print('== filter_snp == ')
    
    pdb.set_trace()
    #with open(count_altInfo_address) as count_altInfo_file:
    with open(count_altInfo_address) as count_altInfo_file, open(count_abs_address) as count_file:
        #for line in count_altInfo_file:
        for x, y in izip(count_altInfo_file, count_file):
            
            pdb.set_trace()
            
            x = line.split()            
            
            if len(x)>=3:
                #pdb.set_trace()
                p = int(x[1].split(',')[0])
                p_e = float(x[1].split(',')[1])
                
                flag_remove = False
                if p in found_snps:
                    
                    for itm in x[2:]:
                        q = int(itm.split(',')[0])
                        q_e = float(itm.split(',')[1])
                        if q in found_snps and q_e>p_e: #a potential shadow snp
                            flag_remove = True
                            num_removed_snps += 1
                            #print('p=%d, p_e=%f, q=%d, q_e=%f (num_removed_snps=%d)'%(p, p_e, q, q_e, num_removed_snps))
                            break
                    #end for itm in x[2:]:
                    
                    if flag_remove==True:
                        del found_snps[p]                                
                        #pdb.set_trace() 
    #end with
    
    filt_snp_res_file = open(filt_snp_res_address, 'w')
    
    for itm in found_snps.items():
        #pdb.set_trace()
        gp = itm[0]
        rB = itm[1][0][0]
        tB = itm[1][0][1]
        filt_snp_res_file.write(str(gp)+'\t')
        filt_snp_res_file.write(rB + '\t-->\t')
        filt_snp_res_file.write(tB + '\t\n')
    
    filt_snp_res_file.close()
    
    return

def analyze_fp(m_snps_cd,
               m_snps_md,
               p_snps_cd,
               p_snps_md,
               r_snps_fp,
               count_altInfo_address): # check shadow snp from counts_altInfo.txt
    
    snps_dic = {}
    for itm in m_snps_cd:
        snps_dic[itm[0]]=itm[1]
    for itm in m_snps_md:
        snps_dic[itm[0]]=itm[1]
    for itm in p_snps_cd:
        snps_dic[itm[0]]=itm[1]
    for itm in p_snps_md:
        snps_dic[itm[0]]=itm[1]
    
    fp_dic = {}
    for itm in r_snps_fp:
        fp_dic[itm[0]]=itm[1]
        
    #pdb.set_trace()
    
    print('== analyze_fp == shadow snp')
    num_fp_from_shadow_snp = 0
    with open(count_altInfo_address) as count_altInfo_file:
        for line in count_altInfo_file:
            x = line.split()            
            
            if len(x)>=3:
                #pdb.set_trace()
                p = int(x[1].split(',')[0])
                
                if p in fp_dic:
                
                    for itm in x[2:]:
                        q = int(itm.split(',')[0])
                        if p in fp_dic and q in snps_dic: #shadow snp
                            num_fp_from_shadow_snp += 1
                            print('[' + '; '.join(x[1:])+'] [fp=%d true_snp=%d] (# fp from shadow snp=%d)'%(p,q, num_fp_from_shadow_snp))
                            break                            
                            #pdb.set_trace() 
    
    print('\n\n')
    print('== analyze_fp == 0 lambda')
    
    num_fp_from_0_lambda0 = 0
    with open(count_altInfo_address) as count_altInfo_file:
        for line in count_altInfo_file:
            x = line.split()
            
            gp = int(x[1].split(',')[0])
            el = float(x[1].split(',')[1])
            
            if gp in fp_dic and el==0:
                num_fp_from_0_lambda0 += 1
                print('[' + '; '.join(x[1:])+'] (# fp from 0 lambda =%d)'%(num_fp_from_0_lambda0))
            
    print('== analyze_fp == end')
    
    return
                   

if __name__ == "__main__":
    
    
    #pdb.set_trace()
    
    #Default_Ref_Path = '/home/olivo/Downloads/0827/server_data_0827_SNP1k_Reads10M/'
    #Default_Ref_Path = '/home/sreeramkannan/singleCell/SNP_Calling_Summer15/data_0814/' #data_0827_SNP1k_Reads10M/'
    Default_Ref_Path = '/home/sreeramkannan/singleCell/SNP_Calling_Summer15/data_SNP1k_Reads10M/'
    
    #snp_res_address = Default_Ref_Path + '/caller_output_snp_found_cross_check_split_trimS3.txt' #'/data_GATK/GATK_out/raw_variants.vcf'
    #snp_res_address = Default_Ref_Path + '/caller_output_snp_found_filtCovQt90.txt'
    #snp_res_address = Default_Ref_Path + '/caller_output_snp_found_cross_check_accepted_hits.txt'
    #snp_res_address = Default_Ref_Path + '/caller_output_snp_found_cross_check_split_trimS3.txt'
    #snp_res_address = Default_Ref_Path + '/caller_output_snp_found_cross_check_star1pass_modCase7.txt'
    #snp_res_address = Default_Ref_Path + '/caller_output_snp_found_alt_mapping_debug.txt'
    #snp_res_address = Default_Ref_Path + '/caller_output_snp_found_alt_mapping_debug_2.txt'
    #snp_res_address = Default_Ref_Path + '/caller_output_snp_found_alt_mapping_debug_0911.txt'
    #snp_res_address = Default_Ref_Path + '/caller_output_snp_found_allCovQt.txt'
    #snp_res_address = Default_Ref_Path + '/caller_output_snp_found_altMap_cross_check_star1pass.txt'
    #snp_res_address = Default_Ref_Path + '/caller_output_snp_found_altMap_cross_check_split.txt'
    #snp_res_address = Default_Ref_Path + '/caller_output_snp_found_altMap_cross_check_dedupped.txt'    
    #snp_res_address = Default_Ref_Path + '/caller_output_snp_found_altMap_cross_check_rg_added_sorted.txt'
    #snp_res_address = Default_Ref_Path + '/caller_output_snp_found_altMap_cross_check_star2pass.txt'
    #snp_res_address = Default_Ref_Path + '/caller_output_snp_found_dirAltMap_cross_check_dedupped.txt'
    #snp_res_address = Default_Ref_Path + '/caller_output_snp_found_dirAltMap'
    #for case2    
    snp_res_address = Default_Ref_Path + '/data_GATK/GATK_out/raw_variants.vcf'    
    
    if isVCF(snp_res_address):
        snp_res_address = extractVCF(snp_res_address) # x.vcf --> x.vcf.txt
    
    pdb.set_trace()
    
    snp_m_address = Default_Ref_Path + '/SNP_m.txt'
    snp_p_address = Default_Ref_Path + '/SNP_p.txt'
    
    snps_list_sorted = get_snp_res_statistics(snp_res_address,
                                              snp_m_address,
                                              snp_p_address)
                                              
    [m_snps_cd, m_snps_md, p_snps_cd, p_snps_md, r_snps_fp] = group_snps(snps_list_sorted)
    
    pdb.set_trace()
    print('# of mis-detection (m & p):%d'%(len(m_snps_md)+len(p_snps_md)))
    print('...# of mis-detection (m):%d'%len(m_snps_md))
    print('...# of mis-detection (p):%d'%len(p_snps_md))
    
    print('# of false-positive:%d'%len(r_snps_fp))
    
    print('# of correct-detection (m & p):%d'%(len(m_snps_cd)+len(p_snps_cd)))
    print('...# of correct-detection (m):%d'%len(m_snps_cd))
    print('...# of correct-detection (p):%d'%len(p_snps_cd))
    
    pdb.set_trace()
    #dmp results
    #res_file = open(Default_Ref_Path + '/snp_res_cross_check_star1pass_modCase7.txt', 'w')
    #res_file = open(Default_Ref_Path + '/snp_res_split_trimS3.txt', 'w')
    #res_file = open(Default_Ref_Path + '/snp_res_alt_mapping_debug_2.txt', 'w')
    #res_file = open(Default_Ref_Path + '/snp_res_allCovQt.txt', 'w')
    #res_file = open(Default_Ref_Path + '/snp_res_alt_mapping_debug_0911.txt', 'w')
    #res_file = open(Default_Ref_Path + '/snp_res_altMap_cross_check_star1pass.txt', 'w')
    #res_file = open(Default_Ref_Path + '/snp_res_altMap_cross_check_split.txt', 'w')
    #res_file = open(Default_Ref_Path + '/snp_res_altMap_cross_check_dedupped.txt', 'w')
    #res_file = open(Default_Ref_Path + '/snp_res_altMap_cross_check_rg_added_sorted.txt', 'w')
    #res_file = open(Default_Ref_Path + '/snp_res_altMap_cross_check_star2pass.txt', 'w')
    #res_file = open(Default_Ref_Path + '/snp_res_dirAltMap_cross_check_dedupped.txt', 'w')
    res_file = open(Default_Ref_Path + '/snp_res_dirAltMap.txt', 'w')
    
    res_file.write('#========== Input Files:\n')
    res_file.write('' + snp_res_address+'\n')
    res_file.write('' + snp_m_address+'\n')
    res_file.write('' + snp_p_address+'\n\n')
    
    res_file.write('#========== Results Summary:\n')
    res_file.write('num of mis-detection (m & p):%d\n'%(len(m_snps_md)+len(p_snps_md)))
    res_file.write('...num of mis-detection (m):%d\n'%len(m_snps_md))
    res_file.write('...num of mis-detection (p):%d\n'%len(p_snps_md))    
    res_file.write('num of false-positive:%d\n'%len(r_snps_fp))    
    res_file.write('num of correct-detection (m & p):%d\n'%(len(m_snps_cd)+len(p_snps_cd)))
    res_file.write('...num of correct-detection (m):%d\n'%len(m_snps_cd))
    res_file.write('...num of correct-detection (p):%d\n\n'%len(p_snps_cd))
    
    res_file.write('#========== Results (m SNP - correct detection):\n')
    dmp_snp_res(m_snps_cd, res_file)
    
    res_file.write('#========== Results (m SNP - miss detection):\n')
    dmp_snp_res(m_snps_md, res_file)
    
    res_file.write('#========== Results (p SNP - correct detection):\n')
    dmp_snp_res(p_snps_cd, res_file)
    
    res_file.write('#========== Results (p SNP - miss detection):\n')
    dmp_snp_res(p_snps_md, res_file)
    
    res_file.write('#========== Results (res SNP - false positive):\n')
    dmp_snp_res(r_snps_fp, res_file)        
    
    res_file.close()
    
    """
    snp_res_address = Default_Ref_Path + caller_op_snp_found_fn
    snp_res_stat_fn = '/snp_res_stat.txt'    
    do_snp_res_statistics(snp_res_address, SNP_address_m, SNP_address_p, snp_res_stat_fn)
    """