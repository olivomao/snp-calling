# -*- coding: utf-8 -*-
"""
Created on Wed Sep 30 12:40:46 2015

@author: olivo
- 10/14: merge counts_altInfo
         
         1) merge_count_x_altInfo: /count_x_altInfo/count_x_altInfo[00~19 *split_sam]_region_??.txt ==>
                                   /count_y_altInfo/count_y_altInfo??.txt
            -- merge_line_list_altInfo
         2) merge_count_y_altInfo: /count_y_altInfo/count_y_altInfo??[00~19 *region].txt ==> /count_<z>_altInfo.txt

"""

import os.path
import subprocess
import pdb
from read import *
from progress.bar import Bar
import sys

def merge_line_list(line_list, op_file):
    
    #pdb.set_trace()
    idx = 0
    line = line_list[0].split()
    cnt = int(line[0])
    ref_pos = int(line[1])
    R = line[2]
    L0 = float(line[3])
    [Na, Nc, Ng, Nt] = [0, 0, 0, 0]
    read_info = []
    
    if len(line)>=8:
        [a, c, g, t] = [int(line[4]), int(line[5]), int(line[6]), int(line[7])]
        n = a + c + g + t
        Na += a
        Nc += c
        Ng += g
        Nt += t
        if n > 0:
            read_info = read_info + line[8:8+n]
    
    for idx in range(1, len(line_list)):
        line = line_list[idx].split()
        if cnt != int(line[0]) or ref_pos != int(line[1]) or R != line[2] or L0 != float(line[3]):
            print('exception at merge_line_list: inconsistent line_list')
            pdb.set_trace()
        if len(line)>=8:
            [a, c, g, t] = [int(line[4]), int(line[5]), int(line[6]), int(line[7])]
            n = a + c + g + t
            Na += a
            Nc += c
            Ng += g
            Nt += t
            if n > 0:
                read_info = read_info + line[8:8+n]
    
    res_str = '%d\t%d\t%s\t%02f'%(cnt, ref_pos, R, L0)
    if len(read_info)>0:
        res_str = '%s\t%d\t%d\t%d\t%d'%(res_str, Na, Nc, Ng, Nt)
        for r in read_info:
            res_str = '%s\t%s'%(res_str, r)
    op_file.write(res_str+'\n')
    
    #if len(read_info)>0:
    #    pdb.set_trace()
    
    return

"""
 /count_y_altInfo/count_y_altInfo[00~19 *region].txt ==> /count_<z>_altInfo.txt
 for fp analysis
 
"""
def merge_count_y_altInfo(in_file_pre, num_p, count_altInfo_address):
    
    count_y_address_list = []
    for i in range(num_p):
        addr = in_file_pre + '%02d'%i + '.txt'
        count_y_address_list.append(addr)
    #pdb.set_trace()
    
    ofile = open(count_altInfo_address, 'w+')
    for i in range(num_p):
        with open(count_y_address_list[i]) as infile:
            for line in infile:
                ofile.write(line) #+'\n')
        print('%s merged'%count_y_address_list[i])
    ofile.close()
    #pdb.set_trace()    
    return
    
def merge_count_x(count_x_dir, #split_sam_dir + count_split_dir0,
                  count_x_fn_pre, #'count_x',
                  num_count_x, #num_split_sam
                  num_p, #num of threads
                  count_y_dir, #Default_Ref_Path + tophat_dir + count_split_dir1,
                  count_y_fn_pre): #'count_y'
    #pdb.set_trace()
    
    cmd = 'mkdir -p '+\
          count_y_dir
    subprocess.call(cmd, shell=True)
    
    count_x_address_list = []
    for i in range(num_count_x):
        addr = count_x_dir + count_x_fn_pre + '%02d'%i + '.txt'
        count_x_address_list.append(addr)
    
    #pdb.set_trace()
    count_x_file_list = []
    for i in range(num_count_x):
        f = open(count_x_address_list[i], 'r')
        count_x_file_list.append(f)
        
    num_tot_pos = sum(1 for line in open(count_x_address_list[0]))
    num_pos_per_count_y = int(num_tot_pos/num_p)+1
    
    count_y_idx = 0
    curr_count_y_addr = count_y_dir+count_y_fn_pre+'%02d'%count_y_idx+'.txt'
    curr_count_y_file = open(curr_count_y_addr, 'w')
    num_lines_processed = 0
    
    #pdb.set_trace()
    
    #bar = Bar('Processing', max=num_tot_pos)
    for i in range(num_tot_pos):
        #bar.next()
        
        #read lines
        line_list = []
        for j in range(num_count_x):
            line_list.append(count_x_file_list[j].readline())
        num_lines_processed += 1
        
        #pdb.set_trace()
        merge_line_list(line_list, curr_count_y_file)
        
        if num_lines_processed > num_pos_per_count_y:
            curr_count_y_file.close()
            print('%s generated (%d lines)'%(curr_count_y_addr, num_lines_processed))
            count_y_idx += 1
            curr_count_y_addr = count_y_dir+count_y_fn_pre+'%02d'%count_y_idx+'.txt'
            curr_count_y_file = open(curr_count_y_addr, 'w')
            num_lines_processed = 0
    
    #bar.finish()
    
    #pdb.set_trace()
    
    if num_lines_processed <= num_pos_per_count_y:
        curr_count_y_file.close()
        print('%s generated (%d lines)'%(curr_count_y_addr, num_lines_processed))    
    
    for i in range(num_count_x):
        count_x_file_list[i].close()
    
    #pdb.set_trace()
    
    return

def merge_line_list2(line_list, op_file):
    
    pt1 = line_list[0].split()[0:4]
    
    pt2 = [0, 0, 0, 0] #a c g t
    pt3 = []
    for i in range(len(line_list)):
        line = line_list[i].split()
        if(len(line))>8:
            pt2[0] += int(line[4])
            pt2[1] += int(line[5])
            pt2[2] += int(line[6])
            pt2[3] += int(line[7])
            pt3 += line[8:]
    
    if sum(pt2)>0:
        #pdb.set_trace()
        res = '\t'.join(pt1) + '\t%d\t%d\t%d\t%d\t'%(pt2[0],pt2[1],pt2[2],pt2[3]) + '\t'.join(pt3)
    else:
        res = '\t'.join(pt1)
    op_file.write(res+'\n')
    
    """
    #pdb.set_trace()
    idx = 0
    line = line_list[0].split()
    cnt = int(line[0])
    ref_pos = int(line[1])
    R = line[2]
    L0 = float(line[3])
    [Na, Nc, Ng, Nt] = [0, 0, 0, 0]
    read_info = []
    
    if len(line)>=8:
        #pdb.set_trace()
        [a, c, g, t] = [int(line[4]), int(line[5]), int(line[6]), int(line[7])]
        n = a + c + g + t
        Na += a
        Nc += c
        Ng += g
        Nt += t
        if n > 0:
            read_info = read_info + line[8:8+n]
        if n>=2:
            pdb.set_trace()
    
    for idx in range(1, len(line_list)):
        line = line_list[idx].split()
        if cnt != int(line[0]) or ref_pos != int(line[1]) or R != line[2] or L0 != float(line[3]):
            print('exception at merge_line_list: inconsistent line_list')
            pdb.set_trace()
        if len(line)>=8:
            #pdb.set_trace()
            [a, c, g, t] = [int(line[4]), int(line[5]), int(line[6]), int(line[7])]
            n = a + c + g + t
            if n>=2:
                pdb.set_trace()
            Na += a
            Nc += c
            Ng += g
            Nt += t
            if n > 0:
                read_info = read_info + line[8:8+n]
    
    res_str = '%d\t%d\t%s\t%02f'%(cnt, ref_pos, R, L0)
    if len(read_info)>0:
        res_str = '%s\t%d\t%d\t%d\t%d'%(res_str, Na, Nc, Ng, Nt)
        for r in read_info:
            res_str = '%s\t%s'%(res_str, r)
    op_file.write(res_str+'\n')
    
    #if len(read_info)>0:
    #    pdb.set_trace()
    """
    
    return
    
def merge_count_x2(count_x_dir, #split_sam_dir + count_split_dir0,
                   count_x_fn1, #pre
                   count_x_fn2, #suffix
                   num_f, #num_count_x to process
                   count_y_dir, #Default_Ref_Path + tophat_dir + count_split_dir1,
                   count_y_fn): #'count_y'
    
    #pdb.set_trace()        
    
    count_x_address_list = []
    for i in range(num_f):
        addr = count_x_dir + count_x_fn1 + '%02d_'%i + count_x_fn2
        #pdb.set_trace()
        count_x_address_list.append(addr)
    
    #pdb.set_trace()
    count_x_file_list = []
    for i in range(num_f):
        f = open(count_x_address_list[i], 'r')
        count_x_file_list.append(f)
        
    num_tot_pos = sum(1 for line in open(count_x_address_list[0]))
    
    curr_count_y_addr = count_y_dir+count_y_fn
    #pdb.set_trace()
    curr_count_y_file = open(curr_count_y_addr, 'w')
    
    for i in range(num_tot_pos):
        line_list = []
        for j in range(num_f):
            line_list.append(count_x_file_list[j].readline())
        merge_line_list2(line_list, curr_count_y_file)
    
    curr_count_y_file.close()
    print('%s generated (%d lines)'%(curr_count_y_addr, num_tot_pos))    
    
    for i in range(num_f):
        count_x_file_list[i].close()
    
    #pdb.set_trace()
    
    return
    
def merge_line_list_altInfo(line_list, curr_count_y_file):
    #pdb.set_trace()

    shift = 0
    if line_list[0].split()[0][0]=='[':
        shift = 1
    
    rel_pos = int(line_list[0].split()[shift])
    itm = line_list[0].split()[1+shift]
    res = ['%d'%rel_pos, itm]
    
    alt_pos_list = []
    for line in line_list:
        itms = line.split()
        if len(itms)>2+shift:
            
            #pdb.set_trace()
            #res = res + itms[2+shift:]
        
            #to avoid duplicates
            candidates = itms[2+shift:]
            for c in candidates:
                gp = int(c.split(',')[0])
                #pdb.set_trace()
                if gp not in alt_pos_list:
                    #pdb.set_trace()
                    alt_pos_list.append(gp)
                    res = res + [c]
                else:
                    #pdb.set_trace()
                    continue
                    
    #pdb.set_trace()
    
    res_str = '\t'.join(res)
    curr_count_y_file.write(res_str+'\n')   
    
    return
    
def merge_count_x_altInfo(count_x_dir, #split_sam_dir + count_split_dir0,
                          count_x_fn1, #pre
                          count_x_fn2, #suffix
                          num_f, #num_count_x to process
                          count_y_dir, #Default_Ref_Path + tophat_dir + count_split_dir1,
                          count_y_fn): #'count_y'
    #pdb.set_trace()
    
    count_x_address_list = []
    for i in range(num_f):
        addr = count_x_dir + count_x_fn1 + '%02d_'%i + count_x_fn2
        #pdb.set_trace()
        count_x_address_list.append(addr)
    
    #pdb.set_trace()
    
    count_x_file_list = []
    for i in range(num_f):
        f = open(count_x_address_list[i], 'r')
        count_x_file_list.append(f)
    
    #pdb.set_trace()
    
    num_tot_pos = sum(1 for line in open(count_x_address_list[0]))
    
    curr_count_y_addr = count_y_dir+count_y_fn
    #pdb.set_trace()
    curr_count_y_file = open(curr_count_y_addr, 'w')
    
    for i in range(num_tot_pos):
        line_list = []
        for j in range(num_f):
            line_list.append(count_x_file_list[j].readline())
        merge_line_list_altInfo(line_list, curr_count_y_file)
    
    curr_count_y_file.close()
    print('%s generated (%d lines)'%(curr_count_y_addr, num_tot_pos))    
    
    for i in range(num_f):
        count_x_file_list[i].close()
    
    #pdb.set_trace()
        
    
    return

def seperate_sam_file(in_dir, in_sam, out_dir, num_p):
    
    #pdb.set_trace()
    
    #sort sam file
    orig_sam = in_dir + in_sam
    
    sorted_sam = ''
    if not 'sorted' in in_sam:
        sorted_sam = in_dir + in_sam[:-4] + '_sorted.sam'
        if os.path.exists(sorted_sam)==False:
            sort_command = "sort " + orig_sam + " > " + sorted_sam
            subprocess.call(sort_command, shell=True)
    else:
        sorted_sam = in_dir + in_sam #same as old one
    
    #sorted_sam statistics
    num_reads = 0
    with open(sorted_sam) as sam_file:
        for line in sam_file:
            if line[0] != '@': # not a comment
                num_reads += 1
    
    #pdb.set_trace()
    
    num_max_reads_per_split_sam = int(num_reads/num_p)+1    
    
    subprocess.call('mkdir -p '+out_dir, shell=True)
    
    curr_split_sam_idx = 0
    curr_split_sam_file_addr = out_dir+'split_sam_%02d'%curr_split_sam_idx
    curr_split_sam_file = open(curr_split_sam_file_addr, 'w+')
    
    #counter_num_split_sam_file = 1
    counter_processed_reads = 0
    
    with open(sorted_sam) as sam_file:
        
        read_group = []
        for line in sam_file:
            if line[0] == '@':
                continue
            
            read_id = line.split()[0]
            if not read_group or read_group[-1][0] == read_id:
                read_group.append([read_id, line])
            elif read_group[-1][0] != read_id:
                #process the current read_group
                for itm in read_group:
                    curr_split_sam_file.write('%s'%itm[1])
                    #if len(read_group)>1:
                    #    curr_split_sam_file.write('[g] %s\n'%itm[1])
                    #else:
                    #    curr_split_sam_file.write('[ ] %s\n'%itm[1])
                    counter_processed_reads += 1
                
                if counter_processed_reads > num_max_reads_per_split_sam:
                    print('processed: %s (%d reads)'%(curr_split_sam_file_addr, counter_processed_reads))
                    curr_split_sam_file.close()
                    
                    curr_split_sam_idx += 1
                    curr_split_sam_file_addr = out_dir+'split_sam_%02d'%curr_split_sam_idx
                    curr_split_sam_file = open(curr_split_sam_file_addr, 'w+')
                    
                    counter_processed_reads = 0
                
                read_group = [[read_id, line]]
            
        #process the final read_group
        for itm in read_group:
            curr_split_sam_file.write('%s'%itm[1])
            #if len(read_group)>1:
            #    curr_split_sam_file.write('[g] %s\n'%itm[1])
            #else:
            #    curr_split_sam_file.write('[ ] %s\n'%itm[1])
            counter_processed_reads += 1
        
        print('processed: %s (%d reads)'%(curr_split_sam_file_addr, counter_processed_reads))
        curr_split_sam_file.close()
    
    #pdb.set_trace()
    
    return [out_dir, 'split_sam_', curr_split_sam_idx+1] #used for further processing

if __name__ == "__main__":
    
    if len(sys.argv)<2: #simple debug and test
        """
        in_dir = '/home/sreeramkannan/singleCell/SNP_Calling_Summer15/data_SNP20_Reads100K_para/tophat_out/' #'/home/sreeramkannan/singleCell/SNP_Calling_Summer15/data_SNP1k_Reads10M/tophat_out/'
        in_sam = 'accepted_hits.sam'
        out_dir = in_dir + '/split_sorted_sam/'
        num_p = 20
        seperate_sam_file(in_dir, in_sam, out_dir, num_p)
        """
        
        count_x_dir = '/home/sreeramkannan/singleCell/SNP_Calling_Summer15/data_SNP20_Reads100K_para/tophat_out/split_sorted_sam/count_split_from_sam_split/'
        count_x_fn_pre = 'count_x'
        num_count_x = 20
        num_p = 20
        count_y_dir = '/home/sreeramkannan/singleCell/SNP_Calling_Summer15/data_SNP20_Reads100K_para/tophat_out/count_y/'
        count_y_fn_pre = 'count_y'
        
        merge_count_x(count_x_dir, #split_sam_dir + count_split_dir0,
                      count_x_fn_pre, #'count_x',
                      num_count_x, #num_split_sam
                      num_p, #num of threads
                      count_y_dir, #Default_Ref_Path + tophat_dir + count_split_dir1,
                      count_y_fn_pre) #'count_y'
    else:
        
        #pdb.set_trace()
        """
        'python para_operations.py para_count_merge '+\
                  split_sam_dir+count_split_dir0+' '\
                  'count_x'+' '+\
                  'region_00.txt'+' '+\
                  repr(num_p)+' '+\
                  Default_Ref_Path + tophat_dir + count_split_dir1 + ' '+\
                  'count_y00.txt'
        """
        mode = sys.argv[1]
        if mode=='para_count_merge':
            count_x_dir = sys.argv[2]
            count_x_fn1 = sys.argv[3] #count_x_??_region_00.txt
            count_x_fn2 = sys.argv[4]
            num_f = int(sys.argv[5])
            count_y_dir = sys.argv[6]
            subprocess.call('mkdir -p '+count_y_dir, shell=True)
            count_y_fn = sys.argv[7]
            #pdb.set_trace()
            merge_count_x2(count_x_dir, #split_sam_dir + count_split_dir0,
                           count_x_fn1, #pre
                           count_x_fn2, #suffix
                           num_f, #num_count_x to process
                           count_y_dir, #Default_Ref_Path + tophat_dir + count_split_dir1,
                           count_y_fn) #'count_y'
            #pdb.set_trace()
            
            """
            merge counts_altInfo
            """
            dir_info = [itm for itm in count_x_dir.split('/') if itm != '']
            count_x_dir = '/' + '/'.join(dir_info[:-1]) + '/' + dir_info[-1] + '_altInfo' + '/'
            count_x_fn1 = count_x_fn1 + '_altInfo'
            count_x_fn2 = count_x_fn2
            num_f = num_f
            
            #pdb.set_trace()
            dir_info = [itm for itm in count_y_dir.split('/') if itm != ''] #count_y_dir.split('/')            
            count_y_dir = '/' + '/'.join(dir_info[:-1]) + '/' + dir_info[-1] + '_altInfo' + '/'
            subprocess.call('mkdir -p '+count_y_dir, shell=True)
            count_y_fn = count_y_fn[:-4]
            count_y_fn = count_y_fn[:-2] + '_altInfo' + count_y_fn[len(count_y_fn)-2:] + '.txt'
            
            #pdb.set_trace()
            
            merge_count_x_altInfo(count_x_dir, #split_sam_dir + count_split_dir0,
                                  count_x_fn1, #pre
                                  count_x_fn2, #suffix
                                  num_f, #num_count_x to process
                                  count_y_dir, #Default_Ref_Path + tophat_dir + count_split_dir1,
                                  count_y_fn) #'count_y'
            #pdb.set_trace()
                        
            
        else:
            print('other modes in para_operations.py')
            pdb.set_trace()
        
                