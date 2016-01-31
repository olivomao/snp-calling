"""
a debug/dmp version for count_read_lambda
- unnecessary codes not executed for efficient run
- for certain pos, dmp aligned reads over it, as well as each read's alt mapping info
- 
- 16/1/31 to avoid null item (e.g. read_group[0]==None) in read_group by try & except
"""

from itertools import *
from test import *
from read import *
import pickle
import time
import sys
import subprocess
import pdb
from debug_MACRO import *
from Address import *
import os.path

ACGT = ("A", "C", "G", "T")
BASES = {"A": 0, "C": 1, "G": 2, "T": 3}
NUM_TO_BASE = {0: "A", 1: "C", 2: "G", 3: "T"}

def update_counts_no_alt_mapping(ref, D, positions, counts, read_group, counts_alt):
    #pdb.set_trace() #to debug
    #lp = len(positions[0])
    for position in positions[0]:
        if position != None:
            genome_pos = position[0]
            if do_debug==True and genome_pos == md_fp_loc:
                #pdb.set_trace()
                
                print('process new read alignment at %d\n'%md_fp_loc)
                md_fp_dmp_file.write('== new read alignment ==\n')
                md_fp_dmp_file.write('read_line:\n%s\n'%read_group[0].line)
                stt_pos = read_group[0].pos
                curr_pos = genome_pos
                try:
                    curr_base = read_group[0].read[curr_pos-stt_pos]
                    #can be misleading
                    #if 'S' in Cigar string, stt_pos and curr_pos corresp to read
                    #with 'S' part skipped
                except:
                    curr_base = '?'
                md_fp_dmp_file.write('stt pos:%d\n'%stt_pos)
                md_fp_dmp_file.write('cur pos:%d\n'%curr_pos)
                md_fp_dmp_file.write('cur read base:%s\n'%curr_base)
                md_fp_dmp_file.write('cur ref base:%s\n'%ref[genome_pos])
                md_fp_dmp_file.write('cur ref L:%.2f\n\n'%float(D.get(genome_pos,0)))
                
                #pdb.set_trace()
                
            read_res = position[1]
            #pdb.set_trace()
            Lsum = "%.2f"%float(D.get(genome_pos, 0))
            read_res[3] = Lsum
            count = counts.setdefault(genome_pos, [])
            count.append(read_res)
            
            #pdb.set_trace()            
            counts_alt.setdefault(genome_pos, {})
    return

def update_counts_with_alt_mapping(ref, D, positions, counts, read_group):
    #pdb.set_trace() #to debug

    for t in zip(*positions): # locii all correlated
    
        if do_debug == True:
            found_loci = 0
        
        exps = []
        
        # filter out None and insertions/deletions
        #f = [read_base for read_base in t if (read_base and read_base[0] in BASES)]
        f = [read_base for read_base in t if (read_base and read_base[1][0] in BASES)]
        if not f:
            continue

        ##merge same pos
        g = []
        pos_list = []
        for read_base in f:
            
            if do_debug==True and read_base[0] == md_fp_loc:
                found_loci = 1
                #pdb.set_trace()
                print('process new read alignment at %d\n'%md_fp_loc)
                md_fp_dmp_file.write('== new read alignment ==\n')
                N_rg = len(read_group)
                for i_rg in range(N_rg):
                    md_fp_dmp_file.write('read_line:\n%s\n'%read_group[i_rg].line)
                    try:
                        stt_pos = read_group[i_rg].pos
                        curr_pos = f[i_rg][0]
                    except:
                        pdb.set_trace()
                    try:
                        curr_base = read_group[i_rg].read[curr_pos-stt_pos]  
                    except:
                        curr_base = '?'
                    md_fp_dmp_file.write('stt pos:%d\n'%stt_pos)
                    md_fp_dmp_file.write('cur pos:%d\n'%curr_pos)
                    md_fp_dmp_file.write('cur read base:%s\n'%curr_base)
                    md_fp_dmp_file.write('cur ref base:%s\n'%ref[curr_pos])
                    md_fp_dmp_file.write('cur ref L:%.2f\n\n'%float(D.get(curr_pos,0)))
                
                #pdb.set_trace()
                
            if read_base[0] not in pos_list:
                g.append(read_base)
                pos_list.append(read_base[0])
            #else:
            #    print('debug here')
            #    pdb.set_trace()
        f = g

        #if len(f)>1:
        #    pdb.set_trace()

        #pdb.set_trace()
        s = zip(*f)
        locii = set(s[0])

        for i in locii:
            exps.append((ref[i], float(D.get(i, 0))))

        L = [sum(e[1] for e in exps if e[0] == x) for x in ACGT]
        Lsum = sum(L)        
        
        if do_debug==True and found_loci==1:
            #pdb.set_trace()
            md_fp_dmp_file.write('L=%s\n'%repr(L))
            md_fp_dmp_file.write('Lsum=%s\n\n'%repr(Lsum))
            
            print('manual debug here (L info for inv/non-inv read alignments)')
            pdb.set_trace()
            
            found_loci = 0

        for gp, r in f:
            #pdb.set_trace()
            #L0 = float(D.get(gp, 0))
            if ref[gp]==r[0]:
                #pdb.set_trace()
                r[2] = "%.2f" % (L[BASES[r[0]]]-float(D.get(gp,0))) #excluding lambda zero
            else:
                #if len(f)>1:
                #    pdb.set_trace()
                #check read_group
                r[2] = "%.2f" % L[BASES[r[0]]] #including lambda zero
            r[3] = "%.2f" % Lsum
            count = counts.setdefault(gp, [])
            count.append(r)    
    return

# function needed for dir_alt_map
def get_directions(read_group):
    directions = [] #record direction info +/- of aligned reads from read_group
    num_pos = 0
    num_neg = 0

    for r in read_group:
        if r.is_reversed():
            directions.append('-')
            num_neg += 1
        else:
            directions.append('+')
            num_pos += 1
    
    return [directions, num_neg, num_pos]

def get_comp_base(base):
    res = ''
    if   base == 'A':
        res = 'T'
    elif base == 'C':
        res = 'G'
    elif base == 'T':
        res = 'A'
    elif base == 'G':
        res = 'C'
    else:
        print('exception at get_comp_base, base=%s'%base)
        pdb.set_trace()
    return res
# function needed for dir_alt_map
# read_positions = [(42841071, ['T', 'I', '0', '0']), (42841072, ['G', 'I', '0', '0']), ...]
# => 100 items, may contain 'None' object generated in increament_read
def re_order_positions(read_positions, direction):
    res = []
    if direction == '+':
        L = len(read_positions)
        for i in range(L):
            rp = read_positions[i]
            if rp:
                pos = rp[0]
                read_res = rp[1]
                res.append([pos, read_res, '+'])
            else: #None
                #pdb.set_trace()
                res.append(rp)
    else:
        L = len(read_positions)
        for i in range(L):
            rp = read_positions[L-1-i]
            if rp:
                pos = rp[0]
                read_res = rp[1]
                res.append([pos, read_res, '-'])
            else:
                #pdb.set_trace()
                res.append(rp)
    return res
    
# function needed for dir_alt_map
def update_counts_with_dir_alt_mapping(ref, D, positions, counts, read_group, directions=[], counts_alt=None):
    
    #pdb.set_trace() #to debug
    
    for i in range(len(directions)):
        positions[i] = re_order_positions(positions[i], directions[i])
    
    #pdb.set_trace()

    for t in zip(*positions): # locii all correlated
    
        if do_debug == True:
            found_loci = 0
        
        exps = []
        
        # filter out None and insertions/deletions
        f = [read_base for read_base in t if (read_base and read_base[1][0] in BASES)]
        if not f:
            continue

        ##merge same pos
        #pdb.set_trace()
        
        g = []
        pos_list = []
        for read_base in f:
            
            if do_debug==True and read_base[0] == md_fp_loc:
                found_loci = 1
                #pdb.set_trace()
                print('process new read alignment at %d\n'%md_fp_loc)
                md_fp_dmp_file.write('== new read alignment ==\n')
                N_rg = len(read_group)
                for i_rg in range(N_rg):
                    md_fp_dmp_file.write('read_line:\n%s\n'%read_group[i_rg].line)
                    try:
                        stt_pos = read_group[i_rg].pos
                        curr_pos = f[i_rg][0]
                    except:
                        md_fp_dmp_file.write('[exception]\n')
                        continue
                    
                    try:
                        curr_base = read_group[i_rg].read[curr_pos-stt_pos]  
                    except:
                        curr_base = '?'
                    md_fp_dmp_file.write('stt pos:%d\n'%stt_pos)
                    md_fp_dmp_file.write('cur pos:%d\n'%curr_pos)
                    md_fp_dmp_file.write('cur read base:%s\n'%curr_base)
                    md_fp_dmp_file.write('cur ref base:%s\n'%ref[curr_pos])
                    md_fp_dmp_file.write('cur ref L:%.2f\n\n'%float(D.get(curr_pos,0)))
                
                #pdb.set_trace()
                
            if read_base[0] not in pos_list:
                g.append(read_base)
                pos_list.append(read_base[0])
            #else:
            #    print('debug here')
            #    pdb.set_trace()
        f = g

        #if len(f)>1:
        #    pdb.set_trace()

        #pdb.set_trace()
        """
        to handle dir alt map here carefully
        """
        
        for gp, r, d in f:
            exps.append((ref[gp].upper(), float(D.get(gp, 0)), d)) #e.g. (A, 120, '+')
        
        L_pos = [sum(e[1] for e in exps if e[0]==x and e[2]=='+') for x in ACGT]
        L_neg = [sum(e[1] for e in exps if e[0]==x and e[2]=='-') for x in ACGT]
        Lsum = sum(L_pos) + sum(L_neg)
        
        #if do_debug==True and found_loci==1:
            #pdb.set_trace()
        #    md_fp_dmp_file.write('L=%s\n'%repr(L))
        #    md_fp_dmp_file.write('Lsum=%s\n\n'%repr(Lsum))
            
        #    print('manual debug here (L info for inv/non-inv read alignments)')
        #    pdb.set_trace()
            
        #    found_loci = 0

        
        for gp, r, d in f:
            
            ref_base = ref[gp].upper()
            aligned_base = r[0].upper()
            
            if d=='+':             
                Lsum_y = L_pos[BASES[aligned_base]] + L_neg[BASES[get_comp_base(aligned_base)]]                
            else:# d=='-'
                Lsum_y = L_pos[BASES[get_comp_base(aligned_base)]] + L_neg[BASES[aligned_base]]
            L0 = float(D.get(gp,0))
            
            if aligned_base == ref_base:
                r[2] = "%.2f" % (Lsum_y-L0) #excluding lambda zero
            else:
                r[2] = "%.2f" % Lsum_y #including lambda zero
            r[3] = "%.2f" % Lsum
            
            count = counts.setdefault(gp, [])
            count.append(r)
        
        """
        add alt mapping info in counts_alt
        """
        
        for gp1, r1, d1 in f:
            
            altCount = counts_alt.setdefault(gp1, {})
            
            for gp2, r2, d2 in f:                
                if gp2 != gp1 and gp2 not in altCount:
                    altCount[gp2] =  float(D.get(gp2, 0))
            
            #pdb.set_trace()
    return
    
def generate_count_file(reads, ref_address, cov_address, count_fn='count.txt', count_dir='', num_p=1, seperate_regions=False):
    """Generate convenient count files

    Inputs:
    example sam_address: "../Data/G_30000000-1/read_l100.sam"
    example ref_address: "../Data/Chr15.fa"
    example cov_address: "../Data/G_30000000-1/coverage.txt")

    Outputs (all of these in the /output directory):
    - counts.txt
    exon_pos  reference_pos  ref_base  lambda_x  N_A  N_C  N_G  N_T  [reads]
    [reads_j] = b_j,q_j,lambda_j,lambda_sumj
    """
    #pdb.set_trace() #debug
    start_time = time.clock() # timer

    # initialize reference
    ref = ''
    with open(ref_address,'rU') as ref_file:
        next(ref_file) # skip header lines
        for line in ref_file:
            segment = dna_pattern.search(line)
            if segment:
                segment = segment.group(1)
            if segment == line[:-1]:
                ref = ref + segment
    G = len(ref)
    
    print 'init ref done'
    
    #pdb.set_trace()
    use_RsemRes = False
    if 'rsem' in cov_address:
        use_RsemRes = True

    # initialize global position to expression level map
    exon_pos = [] # exon_pos[i] = j means global position i on the genome
    #exon_to_global = {}
    j = 0 # corresponds to position j on the transcriptome
    D = {}
    #exon_start = []
    #exon_end = []
    #exon_acc_len = [0]

    if use_RsemRes == False:
        with open(cov_address) as cov_file:
            for line in cov_file:
                x = line.split()
                #exon_start.append(int(x[1]))
                #exon_end.append(int(x[2]) )
                #exon_acc_len.append(int(x[3]))
                expression = x[4]
                if 1: #en_debug == 0:
                    for i in range(int(x[1]), int(x[2])):
                        D[i] = expression
                        exon_pos.append(i)
                        j +=1
                        #exon_to_global[i] = j
                else:
                    #to make data structure more concise
                    continue
    else: #use_RsemRes == True
        with open(cov_address) as cov_file:
            for line in cov_file:
                x = line.split()
                #pdb.set_trace()
                e_pos = int(x[0])
                e_lambda = float(x[1])
                D[e_pos] = e_lambda #expression
                exon_pos.append(e_pos)
                j += 1
                #exon_to_global[e_pos] = j
    
    #pdb.set_trace()
    G_eff = j # the total length of the exons
    print("total reference length", G)
    print("total exons length", G_eff)

    #pdb.set_trace() #debug
    if not "sorted" in reads: # assume not sorted, so sort
        read_name = reads.split("/")[-1]
        if en_debug==0:
            sorted_sam = working_dir + read_dir  + read_name[:-4] + "_sorted.sam"
        else:
            sorted_sam = reads[0:len(reads)-len(read_name)] + read_name[:-4] + "_sorted.sam"
        
        #pdb.set_trace()
        if os.path.exists(sorted_sam)==True:
            print('sorted_sam: %s exists'%sorted_sam)
        else:
            sort_command = "sort " + reads + " > " + sorted_sam
            print('sorted_sam: %s does not exists'%sorted_sam)
            print('sorted_sam: gen cmd: %s'%sort_command)
            subprocess.call(sort_command, shell=True)
    else:
        print("assuming " + reads + " is already sorted")
        sorted_sam = reads

    def count(sam_address):
        
        #pdb.set_trace()
        
        """Takes the given reads at sam_address and produces a count of bases
        D is map from genome location to rna location -
            if D is provided, it is assumed to be RNA reads.
        """
        #pdb.set_trace()
        # Locus not in exon region would be mapped to counts[G_eff]
        # count= [A, C, T, G, I, D, read_map_locations]
        counts = dict()
        counts["insertions"] = dict()
        insertions = counts["insertions"]
        
        counts_alt = dict()
        
        counter = 0 # just to report progress

        def properties(read):
            """returns (is_reversed, is_secondary)"""
            #pdb.set_trace()
            x = read.line
            flag_num = int(x[1])
            # if ((flag_num >> 7) & 1) != ((flag_num >> 6) & 1):
            #     print "0x40 and 0x80 flags are the same, something is wrong"
            #     print read.line
            return ((flag_num >> 4) & 1, (flag_num >> 6) & 1)
            
        def increment_by_read_found_md_fp_pos(read, n, md_fp_pos_ref):
            found_md_fp_pos_res = 0
            
            x = read.line
            flag_num = int(x[1])
            flag = '{0:08b}'.format(flag_num)
            pos = int(x[3]) - 1
            s = num_pattern.findall(x[5])
            splice_n = [int(s[i]) for i in range(len(s))]
            splice_c = chr_pattern.findall(x[5])
            #read = x[9] # may be modified if 'S' happens
            #quality_scores = x[10] # may be modified if 'S' happens
            
            if splice_c[0]=='H':
                splice_n = splice_n[1:]
                splice_c = splice_c[1:]
                #pdb.set_trace()
            if splice_c[-1]=='H':
                splice_n = splice_n[0:-1]
                splice_c = splice_c[0:-1]
                #pdb.set_trace()
                
            if splice_c[0]=='S':
                #pdb.set_trace()
                tmp=splice_n[0]
                splice_n = splice_n[1:]
                splice_c = splice_c[1:]
                #read = read[tmp:]
                #quality_scores = quality_scores[tmp:]                
            if splice_c[-1]=='S':
                #pdb.set_trace()
                tmp=splice_n[-1]
                splice_n = splice_n[0:-1]
                splice_c = splice_c[0:-1]
                #read = read[0:-tmp]
                #quality_scores = quality_scores[0:-tmp]                
            
            #positions = [None] * len(read)
            #pdb.set_trace()
            #if flag[-3] == '1':
                #pdb.set_trace()
            #    continue
            #else:
            if flag[-3] == '0':
                read_pos = -1
                genome_pos = pos - 1                
                
                #if len(splice_n)>1:
                #    pdb.set_trace()
                #if splice_c[0]!="M":
                #    pdb.set_trace()
                while len(splice_n)>0:
                    if splice_c[0]=='M':
                        #pdb.set_trace()
                        for i in range(splice_n[0]):
                            read_pos += 1
                            genome_pos += 1
                            
                            if genome_pos==md_fp_pos_ref:
                                found_md_fp_pos_res = 1
                                break
                            ##count = counts.setdefault(genome_pos, []) #if genome_pos in counts, return its val in counts; or add genome_pos:[] in dict
                            #base = read[read_pos]
                            #quality = quality_scores[read_pos]
                            #pdb.set_trace()
                            #if en_debug == 0:
                            #    read_result = [base.upper(), quality, "1", "1"]
                            #else:
                            #    read_result = [base.upper(), quality, "0", "0"]
                            ##count.append(read_result)
                            #pdb.set_trace()
                            #positions[read_pos] = (genome_pos, read_result)
                    
                    elif splice_c[0]=='N':
                        #print('splice_c[0]==N')
                        #pdb.set_trace()
                        genome_pos += splice_n[0]
                        
                    elif splice_c[0]=='D':
                        #print('splice_c[0]==D')
                        #pdb.set_trace()
                        for _ in range(splice_n[0]):
                            genome_pos += 1
                            
                            if genome_pos==md_fp_pos_ref:
                                found_md_fp_pos_res = 1
                                break
                            ##count = counts.setdefault(genome_pos, [])
                            ##count.append(('D', "I", str(n), "1"))                        
                    elif splice_c[0]=='I':
                        #print('splice_c[0]==I')
                        #pdb.set_trace()
                        for _ in range(splice_n[0]):
                            read_pos += 1
                            ##count = insertions.setdefault(genome_pos, [])
                            ##count.append((read[read_pos].upper(), str(n), quality_scores[read_pos]))
                    elif splice_c[0]=='H':
                        print('splice_c[0]==H (unexpected)')
                        pdb.set_trace()
                    elif splice_c[0]=='S':
                        print('splice_c[0]==S (unexpected)')
                        pdb.set_trace()               
                    else:
                        #misc place here
                        print('splice_c[0]==Misc (unexpected)')
                        pdb.set_trace()                    
                    splice_n = splice_n[1::]
                    splice_c = splice_c[1::]
                    
                    if found_md_fp_pos_res == 1:
                        break
                        
            #return positions
            return found_md_fp_pos_res
            
        def increment_by_read(read, n):
            #pdb.set_trace()
            x = read.line
            flag_num = int(x[1])
            flag = '{0:08b}'.format(flag_num)
            pos = int(x[3]) - 1
            s = num_pattern.findall(x[5])
            splice_n = [int(s[i]) for i in range(len(s))]
            splice_c = chr_pattern.findall(x[5])
            read = x[9] # may be modified if 'S' happens
            quality_scores = x[10] # may be modified if 'S' happens
            
            if splice_c[0]=='H':
                splice_n = splice_n[1:]
                splice_c = splice_c[1:]
                #pdb.set_trace()
            if splice_c[-1]=='H':
                splice_n = splice_n[0:-1]
                splice_c = splice_c[0:-1]
                #pdb.set_trace()
            
            if True: #dir_alt_map == True:
                read_pos = -1
                genome_pos = pos - 1  
                if splice_c[0]=='S':
                    #pdb.set_trace()
                    tmp=splice_n[0]
                    splice_n = splice_n[1:]
                    splice_c = splice_c[1:]
                    #read = read[tmp:]                    
                    #quality_scores = quality_scores[tmp:]
                    read_pos = read_pos + tmp #make sure (1) read len is still 100 so that positions have len 100, 
                                              #can be convenient for later dir_alt_map processing
                                              #(2) read_pos is correct for 'S' condition
                if splice_c[-1]=='S':
                    #pdb.set_trace()
                    #tmp=splice_n[-1]
                    splice_n = splice_n[0:-1]
                    splice_c = splice_c[0:-1]
                    #read = read[0:-tmp]
                    #quality_scores = quality_scores[0:-tmp]             
            else:
                read_pos = -1
                genome_pos = pos - 1  
                if splice_c[0]=='S':
                    #pdb.set_trace()
                    tmp=splice_n[0]
                    splice_n = splice_n[1:]
                    splice_c = splice_c[1:]
                    read = read[tmp:]
                    quality_scores = quality_scores[tmp:]                
                if splice_c[-1]=='S':
                    #pdb.set_trace()
                    tmp=splice_n[-1]
                    splice_n = splice_n[0:-1]
                    splice_c = splice_c[0:-1]
                    read = read[0:-tmp]
                    quality_scores = quality_scores[0:-tmp]   
                         
            
            positions = [None] * len(read)
            #pdb.set_trace()
            #if flag[-3] == '1':
                #pdb.set_trace()
            #    continue
            #else:
            if flag[-3] == '0':
                #read_pos = -1
                #genome_pos = pos - 1                
                
                #if len(splice_n)>1:
                #    pdb.set_trace()
                #if splice_c[0]!="M":
                #    pdb.set_trace()
                while len(splice_n)>0:
                    if splice_c[0]=='M':
                        #pdb.set_trace()
                        for i in range(splice_n[0]):
                            read_pos += 1
                            genome_pos += 1
                            ##count = counts.setdefault(genome_pos, []) #if genome_pos in counts, return its val in counts; or add genome_pos:[] in dict
                            base = read[read_pos]
                            quality = quality_scores[read_pos]
                            #pdb.set_trace()
                            if en_debug == 0:
                                read_result = [base.upper(), quality, "1", "1"]
                            else:
                                read_result = [base.upper(), quality, "0", "0"]
                            ##count.append(read_result)
                            #pdb.set_trace()
                            positions[read_pos] = (genome_pos, read_result)
                    
                    elif splice_c[0]=='N':
                        #print('splice_c[0]==N')
                        #pdb.set_trace()
                        genome_pos += splice_n[0]                    
                    elif splice_c[0]=='D':
                        #print('splice_c[0]==D')
                        #pdb.set_trace()
                        for _ in range(splice_n[0]):
                            genome_pos += 1
                            ##count = counts.setdefault(genome_pos, [])
                            ##count.append(('D', "I", str(n), "1"))                        
                    elif splice_c[0]=='I':
                        #print('splice_c[0]==I')
                        #pdb.set_trace()
                        for _ in range(splice_n[0]):
                            read_pos += 1
                            ##count = insertions.setdefault(genome_pos, [])
                            ##count.append((read[read_pos].upper(), str(n), quality_scores[read_pos]))
                    elif splice_c[0]=='H':
                        print('splice_c[0]==H (unexpected)')
                        pdb.set_trace()
                    elif splice_c[0]=='S':
                        print('splice_c[0]==S (unexpected)')
                        pdb.set_trace()               
                    else:
                        #misc place here
                        print('splice_c[0]==Misc (unexpected)')
                        pdb.set_trace()                    
                    splice_n = splice_n[1::]
                    splice_c = splice_c[1::]
                        
            return positions

        def increment_by_read_orig(read, n): #try to modify this function
            #pdb.set_trace()
            x = read.line
            flag_num = int(x[1])
            flag = '{0:08b}'.format(flag_num)
            pos = int(x[3]) - 1
            s = num_pattern.findall(x[5])
            splice_n = [int(s[i]) for i in range(len(s))]
            splice_c = chr_pattern.findall(x[5])
            read = x[9]
            quality_scores = x[10]
            positions = [None] * len(read)
            #pdb.set_trace()
            #if flag[-3] == '1':
                #pdb.set_trace()
            #    continue
            #else:
            if flag[-3] == '0':
                read_pos = -1
                genome_pos = pos - 1
                #if len(splice_n)>1:
                #    pdb.set_trace()
                #if splice_c[0]!="M":
                #    pdb.set_trace()
                while len(splice_n)>0:
                    for i in range(splice_n[0]):
                        read_pos += 1
                        genome_pos += 1
                        count = counts.setdefault(genome_pos, []) #if genome_pos in counts, return its val in counts; or add genome_pos:[] in dict
                        base = read[read_pos]
                        quality = quality_scores[read_pos]
                        #pdb.set_trace()
                        if en_debug == 0:
                            read_result = [base.upper(), quality, "1", "1"]
                        else:
                            read_result = [base.upper(), quality, "0", "0"]
                        count.append(read_result)
                        #pdb.set_trace()
                        positions[read_pos] = (genome_pos, read_result)
                    splice_n = splice_n[1::]
                    splice_c = splice_c[1::]
                    while splice_n and splice_c[0] in ["N", "D", "I"]:
                        if splice_c[0]=='N':
                            genome_pos += splice_n[0]
                            splice_n = splice_n[1:]
                            splice_c = splice_c[1:]
                        if splice_c[0]=='D':
                            #pdb.set_trace()
                            for _ in range(splice_n[0]):
                                genome_pos += 1
                                count = counts.setdefault(genome_pos, [])
                                count.append(('D', "I", str(n), "1"))
                            splice_n = splice_n[1:]
                            splice_c = splice_c[1:]
                        if splice_c[0]=='I':
                            #pdb.set_trace()
                            for _ in range(splice_n[0]):
                                read_pos += 1
                                count = insertions.setdefault(genome_pos, [])
                                count.append((read[read_pos].upper(), str(n), quality_scores[read_pos]))
                            splice_n = splice_n[1:]
                            splice_c = splice_c[1:]
                        
            return positions

        if en_check_multi_repeats==1:
            num_multi_repeats_cases=0
            tmp_file = open(Default_Ref_Path+'/dmp_multi_repeats_cases.txt', 'w+')
        with open(sam_address) as sam_file:
            if do_debug_for_loc_pc == True:            
                pdb.set_trace()
            if en_debug == 0:
                sam_file.readline() # throw away the first line
            read_group = [] # temporarily stores reads
            for line in sam_file:
                if line[0] != '@' and line[0] != '\n': # not a comment
                    read = Read(line)
                    if not read_group or read_group[-1].id == read.id:
                        read_group.append(read)
                    elif read_group[-1].id != read.id:
                        if do_debug_for_loc_pc == True:
                            #check read group of interest
                            found_md_fp_loc = 0
                            for r in read_group:
                                #pdb.set_trace()
                                found_md_fp_loc = increment_by_read_found_md_fp_pos(r, 0, md_fp_loc)
                                if found_md_fp_loc == 1:
                                    break
                            if found_md_fp_loc == 0:
                                read_group = [read]
                                
                                if counter % 10000 == 0:
                                    print(counter, "lines of the .sam file are processed!",
                                          round(time.clock() - start_time, 2))
                                counter += 1
                                #pdb.set_trace()
                                continue
                        
                        # split into forward and backward cases
                        segment_split = [[r for r in read_group if not r.is_first_segment()], [r for r in read_group if r.is_first_segment()]]
                        if en_check_multi_repeats==1:
                            count_i = 1
                            count_j = 1
                            #count_multi = 0
                            #pdb.set_trace()
                            #print('')
                        for split in segment_split:
                            if en_check_multi_repeats==1:
                                count_i = (count_i + 1)%2
                            M = len(split)
                            groups = [split]
                            #tmp = [[r for r in split if not r.is_reversed()], [r for r in split if r.is_reversed()]]
                            #if len(tmp[0])>0 and len(tmp[1])>0:
                            #    pdb.set_trace()
                            for read_group in groups:
                                if en_check_multi_repeats==1:
                                    count_j = (count_j + 1)%2
                                    if len(read_group)>1:
                                        num_multi_repeats_cases += 1
                                        tmp_id = read_group[0].id
                                        tmp_file.write('case: %d\n'%num_multi_repeats_cases)
                                        tmp_file.write('%d%d[%d]\n'%(count_i, count_j, len(read_group)))
                                        for tmp_r in read_group:
                                            tmp_file.write(' '.join(tmp_r.line)+'\n')
                                    #count_multi = count_multi + len(read_group)
                                if True: #dir_alt_map == True:
                                    directions_info = get_directions(read_group)
                                    directions = directions_info[0]
                                    num_neg = directions_info[1]
                                    num_pos = directions_info[2]
                                    
                                positions = [increment_by_read(r, M) for r in read_group] #counts updated here
                                if len(positions)==0:
                                    #empty read_group
                                    continue
                                    
                                elif len(positions) == 1:
                                    update_counts_no_alt_mapping(ref, D, positions, counts, read_group, counts_alt)
                                    
                                elif len(positions) > 1: # deal with repeat regions
                                    #print('debug: deal with repeat regions')
                                    #pdb.set_trace()
                                    if True: #dir_alt_map == True:
                                        update_counts_with_dir_alt_mapping(ref, D, positions, counts, read_group, directions, counts_alt)
                                    else:
                                        update_counts_with_alt_mapping(ref, D, positions, counts, read_group)
                                    
                                else:
                                    print('update counts exception: unexpected len(positions)')
                                    pdb.set_trace()
                        #if en_check_multi_repeats==1:
                        #    if count_multi > 1:
                        #        try:
                        #            num_multi_repeats_cases += 1
                        #            #pdb.set_trace()
                        #            tmp = segment_split[0]
                        #            if len(tmp)>0:
                        #                tmp_id = tmp[0].id
                        #            else:
                        #                tmp = segment_split[1]
                        #                tmp_id = tmp[0].id
                        #            print('case %d: %s'%(num_multi_repeats_cases,  tmp_id))
                        #            #pdb.set_trace()
                        #        except:
                        #            pdb.set_trace()
                        #            print('check multi-repeats-cases exception')
                        read_group = [read]

                    if counter % 10000 == 0:
                        print(counter, "lines of the .sam file are processed!",
                            round(time.clock() - start_time, 2))
                    counter += 1
                    
            # deal with final read group after the loop
            segment_split = [[r for r in read_group if not r.is_first_segment()], [r for r in read_group if r.is_first_segment()]]

            for split in segment_split:
                M = len(split)
                groups = [split] #[[r for r in split if not r.is_reversed()], [r for r in split if r.is_reversed()]]
                for read_group in groups: 
                    if True: #dir_alt_map == True:
                        directions_info = get_directions(read_group)
                        directions = directions_info[0]
                        num_neg = directions_info[1]
                        num_pos = directions_info[2]
                            
                    positions = [increment_by_read(r, M) for r in read_group] #counts updated here
                    if len(positions)==0:
                        continue
                    
                    elif len(positions) == 1:
                        update_counts_no_alt_mapping(ref, D, positions, counts, read_group, counts_alt)
                        
                    elif len(positions) > 1: # deal with repeat regions
                        #print('debug: deal with repeat regions')
                        #pdb.set_trace()
                        if True: #dir_alt_map == True:
                            update_counts_with_dir_alt_mapping(ref, D, positions, counts, read_group, directions, counts_alt)
                        else:
                            update_counts_with_alt_mapping(ref, D, positions, counts, read_group)
                        
                    else:
                        print('update counts exception: unexpected len(positions)')
                        pdb.set_trace()                
                """
                for read_group in groups:
                    positions = [increment_by_read(r, M) for r in read_group]
                    if len(positions) > 1: # deal with repeat regions
                        for t in zip(*positions): # locii all correlated
                            exps = []

                            # filter out None and insertions/deletions
                            f = [read_base for read_base in t if read_base and (read_base[0] in BASES)]
                            if not f:
                                continue

                            s = zip(*f)
                            locii = set(s[0])

                            for i in locii:
                                exps.append((ref[i], float(D.get(i, 0))))

                            L = [sum(e[1] for e in exps if e[0] == x) for x in ACGT]
                            Lsum = sum(L)

                            for _, r in f:
                                r[2] = "%.2f" % L[BASES[r[0]]]
                                r[3] = "%.2f" % Lsum
                """
        if en_check_multi_repeats==1:
            tmp_file.close()
        

        print("Done processing " + sam_address)
        return [counts, counts_alt]
    
    def dump(counts, file, num_p=1, seperate_regions=False, counts_alt=None):
        
        def println(i, count, e, of):
                reference_base = ref[i]
                of.write(str(e) + "\t")
                of.write(str(i) + "\t") # reference position
                of.write(reference_base + "\t")
                of.write(str(round(float(D.get(i, "0")), 3))+ "\t") # expression level
                if (count):
                    N_A = sum(1 for c in count if c[0] == "A")
                    N_C = sum(1 for c in count if c[0] == "C")
                    N_G = sum(1 for c in count if c[0] == "G")
                    N_T = sum(1 for c in count if c[0] == "T")
                    of.write(str(N_A) + "\t")
                    of.write(str(N_C) + "\t")
                    of.write(str(N_G) + "\t")
                    of.write(str(N_T) + "\t")
                    of.write("\t".join([",".join(c) for c in count]))
                of.write("\n")
                return
        
        def println2(i, count, e, of):

                #if len(count)>1:
                #    of.write('[%d]\t'%len(count))
                #else:
                #    of.write('[]\t')
            
                of.write(str(e) + "\t")
                of.write(str(i) + ",") # reference position
                of.write(str(round(float(D.get(i, "0")), 3))+ "\t") # expression level
                
                for gp, el in count.iteritems():
                    of.write(str(gp) + ",")
                    of.write(str(el) + "\t")
                
                of.write("\n")
                return
                
        if seperate_regions==False: #non-parallel     
        
            #pdb.set_trace()
            out_address_2 = file[:-4]+'_altInfo.txt' #for counts_alt
            
            #pdb.set_trace()
            #with open(file, 'w+') as outfile:
            outfile = open(file, 'w+') 
            outfile2 = open(out_address_2, 'w+')
            
            for i, j in enumerate(exon_pos):
                x = counts.get(j, [])
                println(j, x, i, outfile)
                y = counts_alt.get(j, {}) #additional alt mapping info
                println2(j, y, i, outfile2)
                
            outfile.close()
            outfile2.close()
            print(file + " written")
            print(out_address_2 + "written")
            # with open('insertions' + file, 'w+') as outfile:
            #     insertions = counts["insertions"]
            #     for pos in insertions.iterkeys():
            #         outfile.write(str(pos) + "\t")
            #         outfile.write(str(round(float(D.get(pos, "0")), 3)))
            #         outfile.write("\t")
            #         outfile.write("\t".join([str(i) for i in map(",".join, insertions[pos])]))
            #         outfile.write("\n")
            # print('insertions' + file + " written")
        else: #seperate_regions==True
            #pdb.set_trace()
            
            curr_f_limit = int(len(exon_pos)/num_p)+1            
            f_pre = file[:-4] #addr without '.txt'
            
            path_info = f_pre.split('/')            
            f_pre2_dir = '/'.join(path_info[:-3]) + '/' + path_info[-2] + '_altInfo/'
            subprocess.call('mkdir -p '+f_pre2_dir, shell=True)
            f_pre2 = f_pre2_dir + path_info[-1][:-2] + '_altInfo' + path_info[-1][len(path_info[-1])-2:]
            #to store counts_alt seperatedly
            
            curr_f_idx = 0
            curr_f_addr = f_pre + '_region_%02d.txt'%curr_f_idx
            curr_f_addr2 = f_pre2 + '_region_%02d.txt'%curr_f_idx
            
            curr_f = open(curr_f_addr, 'w+')
            curr_f2 = open(curr_f_addr2, 'w+')
            
            curr_f_counter = 0
            
            for i, j in enumerate(exon_pos):
                x = counts.get(j, [])
                println(j, x, i, curr_f)
                
                y = counts_alt.get(j, {}) #additional alt mapping info
                println2(j, y, i, curr_f2)
                
                #if len(y)>1:
                #    pdb.set_trace()
                
                curr_f_counter += 1
                if curr_f_counter > curr_f_limit:
                    curr_f.close()
                    curr_f2.close()
                    
                    print('%s (and %s) written (%d pos)'%(curr_f_addr, curr_f_addr2, curr_f_counter))
                    
                    curr_f_idx += 1
                    
                    curr_f_addr = f_pre + '_region_%02d.txt'%curr_f_idx
                    curr_f_addr2 = f_pre2 + '_region_%02d.txt'%curr_f_idx
                    
                    curr_f = open(curr_f_addr, 'w+')
                    curr_f2 = open(curr_f_addr2, 'w+')
                    
                    curr_f_counter = 0
            
            if curr_f_counter <= curr_f_limit:
                curr_f.close()
                curr_f2.close()
                print('%s (and %s) written (%d pos)'%(curr_f_addr, curr_f_addr2, curr_f_counter))
                
        #pdb.set_trace()
        return #dump

    [counts, counts_alt] = count(sorted_sam)

    read_name = reads.split("/")[-1]
    out_dir = reads[0:len(reads)-len(read_name)] + count_dir
    subprocess.call('mkdir -p '+out_dir, shell=True)
    out_address = out_dir + count_fn #for counts    
    
    print('dump: count file')
    pdb.set_trace()
    #dump(counts, out_address, num_p, seperate_regions, counts_alt) 
    del counts
    del counts_alt

"""
#reads = "../../Data/G_30000000-1/read_l100_sorted.sam"
reads = sys.argv[1] # "G_30000000-1/read_l100_sorted.sam"
#reads = "/data/soheil/Chr15/13-4_accepted_hits_chr15.bam"
reference = sys.argv[2] # "G_30000000-1/Chr15.fa"
#reference = "/data/soheil/Chr15/cher15.fa"
coverage = sys.argv[3] # "G_30000000-1/coverage.txt"
resolve_target(reads, reference, coverage)
# AAA
"""

if __name__ == "__main__":
    do_debug = True
    do_debug_for_loc_pc = True #comment some codes for function to run faster

    if do_debug:
        #working_dir = '/home/sreeramkannan/singleCell/SNP_Calling_Summer15/data_0827_SNP1k_Reads10M/'
        #working_dir = '/home/sreeramkannan/singleCell/SNP_Calling_Summer15/data_0814/'
        #working_dir = '/home/olivo/Downloads/0907_count_alt_mapping/bkp/data_0814/'
        #working_dir = '/home/sreeramkannan/singleCell/SNP_Calling_Summer15/data_SNP20_Reads100K_diffMPExp/'
        working_dir = '/home/sreeramkannan/singleCell/SNP_Calling_Summer15/data_SNP20_Reads100K_diffMPExp_case1plus6_altCountModi_2ndrun/'
        #read_dir = "/tophat_out/"
        #read_dir = "data_GATK/1pass/"
        read_dir = "data_GATK/2pass/"

        #reads = working_dir + read_dir + "accepted_hits.sam"
        #reads = working_dir + read_dir + "Aligned.out.sam"
        #reads = working_dir + read_dir + "split.sam"
        reads = working_dir + read_dir + "dedupped.sam"
        reference = working_dir + "Chr15.fa"
        coverage = working_dir + "count_rsem.txt" #"coverage.txt"
        #count_fn = "/count_altMap_cross_check_star1pass_debug.txt"
        #count_fn = "/count_accepted_hits.txt"
        #count_fn = "/count_split.txt"
        count_fn = "/count_dedupped.txt"
        
        #md_fp_loc = 30922906 #26361326 #74536339 #45405535 #28765945 #28635981 #['A', 'G', 'r']
        #md_fp_loc = 20935096 #37188833 #34354849 #44776540
        md_fp_loc = 74362860 #83100618 #82937286 #25333635 #70389192 #43069377
        md_fp_dmp_fn = reads + '.md_fp_dmp_' + repr(md_fp_loc) #e.g. /data_GATK/1pass/dedupped.sam.md_fp_dmp
        md_fp_dmp_file = open(md_fp_dmp_fn, 'w+')

        generate_count_file(reads, reference, coverage, count_fn)
        
        #md_fp_dmp_file.close() close before output of count.txt
