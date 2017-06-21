"""
to be executed mainly instead of debug
- debug code should be cleaned later
- 9/25: modify [ref_base, quality, Ly, Lsum] based on direction of read alignment (tag: dir_alt_map)
- 9/28: tag: dir_alt_map is commented (always True)
- 10/14: add counts_alt (additional info to check shadow snp etc)
- 10/24: rsem (check cov file is ideal or from rsem in generate_count...)
         some variables are commented since no where used (e.g. exon_to_global)
- 16/1/24: counts_alt:
  -- the alt mappings of the reads at the loci are grouped based on the direction and base of the reads
  -- relevant dump (println2) also modified
- 17/2/11: re-visit code, add some comments

"""

#from itertools import *
#from test import *
from read import *
#import pickle
import time
import sys
import subprocess
#import pdb
#from debug_MACRO import *
from Address import *
import os.path
import pdb

ACGT = ("A", "C", "G", "T")
BASES = {"A": 0, "C": 1, "G": 2, "T": 3}
NUM_TO_BASE = {0: "A", 1: "C", 2: "G", 3: "T"}

#positions/D/counts: genome_pos 0-based
#sel_snps_loc: genome_pos 1-based
def update_counts_no_alt_mapping(ref, D, positions, counts, read_group, counts_alt, sel_snp_loc=None, cnt_format=0):
    #pdb.set_trace() #to debug
    #lp = len(positions[0])
    for position in positions[0]:
        if position != None and ((sel_snp_loc is None) or (sel_snp_loc is not None and (position[0]) in sel_snp_loc)): #snp loc 0-based?
            #pdb.set_trace()
            genome_pos = position[0]
            read_res = position[1]
            #pdb.set_trace()
            Lsum = "%.2f"%float(D.get(genome_pos, 0))
            read_res[3] = Lsum

            if cnt_format==0:
                count = counts.setdefault(genome_pos, [])
                count.append(read_res)
            elif cnt_format==1:
                #pdb.set_trace()
                count = counts.setdefault(genome_pos, [])
                bj = read_res[0].upper()
                count.append((bj, str(0)))
            else:
                print('undefined cnt_format:%d'%cnt_format)
                pdb.set_trace()
            
            #pdb.set_trace()            
            counts_alt.setdefault(genome_pos, {})
    return

def update_counts_with_alt_mapping(ref, D, positions, counts, read_group):
    #pdb.set_trace() #to debug

    for t in zip(*positions): # locii all correlated
    
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
def update_counts_with_dir_alt_mapping(ref, D, positions, counts, read_group, directions=[], counts_alt=None, sel_snp_loc=None, cnt_format=0):
    
    #pdb.set_trace() #to debug
    
    for i in range(len(directions)):
        positions[i] = re_order_positions(positions[i], directions[i])
    
    #pdb.set_trace()

    for t in zip(*positions): # locii all correlated

        # e.g. t = ([20363038, ['C', 'I', '0', '0'], '+'], [21369645, ['C', 'I', '0', '0'], '+'])
        # same read covers 20363038 (read base C, read base quality I) and 21369645 (read base C, read base quality I)
        # some part of read can be truncated (e.g. S, H, D, I patterns), so relevant entry in t is None

        #if sum([1 for it in t if it[0]==89399684])>0: #debug
        #    pdb.set_trace()
    
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
            if read_base[0] not in pos_list: #filter out overlapping
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
            exps.append((ref[gp].upper(), float(D.get(gp, 0)), d)) #e.g. (gp's ref base=A, gp's abundance = 120, '+')
        
        L_pos = [sum(e[1] for e in exps if e[0]==x and e[2]=='+') for x in ACGT]
        L_neg = [sum(e[1] for e in exps if e[0]==x and e[2]=='-') for x in ACGT]
        Lsum = sum(L_pos) + sum(L_neg)
        
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

            if (sel_snp_loc is None) or (sel_snp_loc is not None and (gp) in sel_snp_loc):#snp loc 0-based?
                #if gp==84911835 or gp==84911836:
                #    pdb.set_trace()
                #if gp==82986260 and flag_alt==1:
                #    pdb.set_trace()
                #if gp==72958147 and aligned_base=='A': # or gp==28772436:
                #    pdb.set_trace()

                if cnt_format==0:
                    count = counts.setdefault(gp, [])
                    count.append(r)
                elif cnt_format==1:
                    #pdb.set_trace()
                    count = counts.setdefault(gp, [])
                    count.append((aligned_base, str(len(f)-1)))
                else:
                    print('unexpected cnt_format: %d'%cnt_format)
                    pdb.set_trace()
            else:
                continue
        
        """
        add alt mapping info in counts_alt
        """

        #if '+' in directions and '-' in directions:
        #    pdb.set_trace()

        # if read_1  aligned onto gp1 (read base A, dir +), and multiply aligned onto gp2 (read base A, dir +)
        #    read_1' also multiply aligned onto gp1 (read base A, dir +) and gp2 (read base A, dir +)
        # at count_alt[gp1][+A]={gp2:lambda(gp2), ...} number of reads not indicated here
        # 

        if len(f)!=1: #it's possible that f has two same entries at first, and collapsed to one entry at g

            for gp1, r1, d1 in f:

                if (sel_snp_loc is None) or (sel_snp_loc is not None and (gp1) in sel_snp_loc):#snp loc 0-based?
                    #pdb.set_trace()

                    #if gp1==72582568: pdb.set_trace()
                
                    altCount = counts_alt.setdefault(gp1, {})
                
                    newKey = d1 + r1[0].upper() #e.g. +A, -A, +T, -T, +C, -C, +G, -G (these are read bases)
                    sub_altCount = altCount.setdefault(newKey, {})
                
                    for gp2, r2, d2 in f:
                        if gp2 != gp1 and gp2 not in sub_altCount:
                            sub_altCount[gp2] = float(D.get(gp2,0))

                else:
                    continue

    return

def generate_count_file(reads, ref_address, cov_address, count_fn='count.txt', count_dir='', num_p=1, seperate_regions=False, sel_snp_loc=None, isRsemSam=1, cnt_format=0):
    """Generate convenient count files

    Inputs:
    reads: a sam file
    ref_address: a fasta file w/ ref genome
    cov_address: specifies abundance for each genome location
    count_fn and count_dir: output
    num_p and seperate_regions used for parallel purpose

    sel_snps_loc: None or set of snp locations (1-based) for final snp analysis purpose
    isRsemSam: 1(default) or 0. If the sam file provided is an intermediate product through RSEM quantification, we start filter read alignments of low ZW-values

    cnt_format: 0 (default counts.txt output format) 1 (modified counts.txt output format)

    Outputs (all of these in the /output directory):
    
    - counts.txt
    
    --cnt_format 0 (default)
    exon_pos  reference_pos (0-based)  ref_base  lambda_0 (abundance at reference_pos)  N_A  N_C  N_G  N_T  [read_l \t read_2 \t ...]
    reads_j = read_j_base (e.g. y),read_j_base_quality,lambda_y (sum of abundance of alt mapping of read_j_base where those ref base is y),lambda_sum (sum of abundance of all mapping of read_j_base)
    
    --cnt_format 1:
    bj,num_alt_loc (exclude current one). -- num_alt_loc: other locations (excluding current read alignment) the read is mapped onto.
    esp for sensitivity analysis (needs to group snp into different multimapping groups)

    - count_altInfo.txt
    exon_pos  ref_pos,lambda_0  [dir_Base:gp1,ab1,[gp2,ab2,...]][dir_Base:gp1,ab1,[gp2,ab2,...]]...
    -- ref_pos: genome location, 0-based
    -- lambda_0: abundance at ref_pos
    -- dir_Base: +A/-A/+C/-C/+T/-T/+G/-G 
    -- for example, [+A:83141033,57.11,82975427,5.48,] means at ref_pos, there're reads aligned in forward direction w/ read base A, 
       and some or all of there reads have been alternatively mapped at 83141033 w/ abundance 57.11 there, and mapped at 82975427 w/ abundance 5.48 there.

    Variables:

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
    
    print('init ref done')
    
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
            #pdb.set_trace()
            if en_debug == 0:
                sam_file.readline() # throw away the first line
            read_group = [] # temporarily stores reads
            num_read_group = 0
            skip_read_group = {} #key-len(read_group), val-count
            non_skip_read_group = {} #key-len(read_group), val-count
            for line in sam_file:
                if line[0] != '@' and line[0] != '\n': # not a comment
                    read = Read(line)
                    if not read_group or read_group[-1].id == read.id:
                        read_group.append(read)
                    elif read_group[-1].id != read.id:

                        # for rsem sam, discard alignments w/o zw or w/ zw<0.1
                        #if sum([1 for r in read_group if r.ZW<0.1])>0:
                        #    pdb.set_trace()
                        if isRsemSam==1:
                            if sel_snp_loc is None:
                                read_group = [r for r in read_group if r.ZW>=0.1]###rsem sam only
                            else: #for analysis purpose
                                read_group = [r for r in read_group if r.ZW>=0.1]###rsem sam only 

                        #if len(read_group)>0 and read_group[0].id=='uc002bna.2_e_32_chr15_89399479': pdb.set_trace()

                        '''
                        sameZW_flag=True
                        for r_idx in xrange(1,len(read_group)):
                            sameZW_flag=sameZW_flag and (read_group[r_idx-1].ZW==read_group[r_idx].ZW)

                        num_read_group += 1
                        if len(read_group)>1 and sameZW_flag==True:
                            #pdb.set_trace()
                            if len(read_group) in skip_read_group:
                                skip_read_group[len(read_group)]+=1
                            else:
                                skip_read_group[len(read_group)]=1

                            read_group = [read]
                            continue # go back to for line in sam_file for processing another read group
                        else:
                            if len(read_group) in non_skip_read_group:
                                non_skip_read_group[len(read_group)]+=1
                            else:
                                non_skip_read_group[len(read_group)]=1
                        '''

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
                                    update_counts_no_alt_mapping(ref, D, positions, counts, read_group, counts_alt, sel_snp_loc=sel_snp_loc, cnt_format=cnt_format)
                                    
                                elif len(positions) > 1: # deal with repeat regions
                                    #print('debug: deal with repeat regions')
                                    #pdb.set_trace()
                                    if True: #dir_alt_map == True:
                                        update_counts_with_dir_alt_mapping(ref, D, positions, counts, read_group, directions, counts_alt, sel_snp_loc=sel_snp_loc, cnt_format=cnt_format)
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
            
            #pdb.set_trace()

            # for rsem sam, discard alignments w/o zw or w/ zw<0.1
            #if sum([1 for r in read_group if r.ZW<0.1])>0:
            #    pdb.set_trace()
            if isRsemSam==1:
                if sel_snp_loc is None:
                    read_group = [r for r in read_group if r.ZW>=0.1]###rsem sam only
                else: #for analysis purpose
                    read_group = [r for r in read_group if r.ZW>=0.1]###rsem sam only 

            '''
            sameZW_flag=True
            for r_idx in xrange(1,len(read_group)):
                sameZW_flag=sameZW_flag and (read_group[r_idx-1].ZW==read_group[r_idx].ZW)

            num_read_group += 1
            if len(read_group)>1 and sameZW_flag==True:
                #pdb.set_trace()
                if len(read_group) in skip_read_group:
                    skip_read_group[len(read_group)]+=1
                else:
                    skip_read_group[len(read_group)]=1

                #read_group = [read]
                #continue # go back to for line in sam_file for processing another read group
            else:
            '''
            if True:
                if len(read_group) in non_skip_read_group:
                    non_skip_read_group[len(read_group)]+=1
                else:
                    non_skip_read_group[len(read_group)]=1

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
                            update_counts_no_alt_mapping(ref, D, positions, counts, read_group, counts_alt, sel_snp_loc=sel_snp_loc, cnt_format=cnt_format)
                            
                        elif len(positions) > 1: # deal with repeat regions
                            #print('debug: deal with repeat regions')
                            #pdb.set_trace()
                            if True: #dir_alt_map == True:
                                update_counts_with_dir_alt_mapping(ref, D, positions, counts, read_group, directions, counts_alt, sel_snp_loc=sel_snp_loc, cnt_format=cnt_format)
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
    
    def dump(counts, file, num_p=1, seperate_regions=False, counts_alt=None, sel_snp_loc=None, cnt_format=0):
        
        def println(i, count, e, of, cnt_format=0):
                #if sel_snp_loc is not None and (i+1) not in sel_snp_loc:
                #    return

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
                    if cnt_format==0:
                        of.write("\t".join([",".join(c) for c in count]))
                    elif cnt_format==1:
                        #pdb.set_trace()
                        st = '\t'.join([','.join(c) for c in count])
                        of.write(st)
                    else:
                        print('unknown cnt format: %d'%cnt_format)
                        pdb.set_trace()
                of.write("\n")
                return
        
        """
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
        """
                
        def println2(i, count, e, of):

                #if sel_snp_loc is not None and (i+1) not in sel_snp_loc:
                #    return

                #if len(count)>1:
                #    of.write('[%d]\t'%len(count))
                #else:
                #    of.write('[]\t')
            
                of.write(str(e) + "\t")
                of.write(str(i) + ",") # reference position
                of.write(str(round(float(D.get(i, "0")), 3))+ "\t") # expression level
                
                for dB, sub_count in count.iteritems():
                    #dB: direction + aligned base  e.g. +A, -A, +T, -T, +C, -C, +G, -G
                    of.write("[%s:"%(dB))
                    for gp, el in sub_count.iteritems():
                        of.write(str(gp) + ",")
                        of.write(str(el) + ",")
                    of.write("]")
                of.write("\n")
                return
                
        if sel_snp_loc is not None: #non-parallel, restricted to selected snps     
        
            #pdb.set_trace()
            out_address_2 = file[:-4]+'_altInfo.txt' #for counts_alt
            
            #pdb.set_trace()
            #with open(file, 'w+') as outfile:
            outfile = open(file, 'w+') 
            outfile2 = open(out_address_2, 'w+')

            sel_snp_loc_list = sorted(list(sel_snp_loc))
            #pdb.set_trace()
            
            for i in sel_snp_loc_list:
                j = i #j = i-1 #1-based to 0-based #already 0-based?
                #j-genomic location                
                x = counts.get(j, [])
                println(j, x, j, outfile, cnt_format=cnt_format)
                y = counts_alt.get(j, {}) #additional alt mapping info
                println2(j, y, j, outfile2)
                
            outfile.close()
            outfile2.close()
            print(file + " written")
            print(out_address_2 + " written\n")

        elif seperate_regions==False: #non-parallel     
        
            #pdb.set_trace()
            out_address_2 = file[:-4]+'_altInfo.txt' #for counts_alt
            
            #pdb.set_trace()
            #with open(file, 'w+') as outfile:
            outfile = open(file, 'w+') 
            outfile2 = open(out_address_2, 'w+')
            
            for i, j in enumerate(exon_pos):
                #j-genomic location                
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
            
            curr_f_limit = int(len(exon_pos)/num_p)+1 #len(exon_pos): total number of exon bases
            f_pre = file[:-4] #addr without '.txt' file is the out address (e.g. .../count.txt)
            
            path_info = f_pre.split('/')            
            f_pre2_dir = '/'.join(path_info[:-3]) + '/' + path_info[-2] + '_altInfo/' #count and count_alt stored in diff folders
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
    
    if sel_snp_loc is not None:
        out_dir = count_dir
    else:
        out_dir = reads[0:len(reads)-len(read_name)] + count_dir
    subprocess.call('mkdir -p '+out_dir, shell=True)
    out_address = out_dir + count_fn #for counts    
    
    print('dump: count file')
    #pdb.set_trace()
    dump(counts, out_address, num_p, seperate_regions, counts_alt, sel_snp_loc=sel_snp_loc, cnt_format=cnt_format) 
    del counts
    del counts_alt

'''
usage:

#generate count and count_alt selectively in order to analyze snps

python count_read_lambda.py --gen_count_selectively [--format 0 (def),1 (for sens analysis)]
                            -s0 snpFile0 [-s1 snpFile1 etc]
                            --sam samFile [--isRsemSam 1/0 (default 1), ZW=0 filtered]
                            --ref refGenome
                            --cov covFile
                            --countDir countDir
                            --countFn countFn


output is count file and count_alt file

count file format is determined by --format:
- 0: default bj,qj,lambda_bj,lambda_sum
- 1: esp for sensitivity analysis (needs to group snp into different multimapping groups): bj,num_alt_loc (exclude current one)
     -- num_alt_loc: other locations (excluding current read alignment) the read is mapped onto
'''
def gen_count_selectively(args):

    dir_alt_map = True

    i=0
    snp_locations = set() #set of locations (1-based), regardless of snp files, pooled
    while '-s%d'%i in args:
        snpFile = args[args.index('-s%d'%i)+1]
        with open(snpFile, 'r') as sF:
            for line in sF:
                if line[0] != '#' and len(line.split())==4:
                    snp_locations.add(int(line.split()[0]))
        print('%s processed'%snpFile)
        i += 1
    #pdb.set_trace()

    samFile = args[args.index('--sam')+1]
    refGenome = args[args.index('--ref')+1]
    covFile = args[args.index('--cov')+1]
    countDir = args[args.index('--countDir')+1]
    countFn = args[args.index('--countFn')+1]

    if '--isRsemSam' in args:
        isRsemSam = int(args[args.index('--isRsemSam')+1])
    else:
        isRsemSam = 1 #s.t. we filter low ZW values

    if '--format' in args:
        cnt_format = int(args[args.index('--format')+1])
    else:
        cnt_format = 0

    print('gen_count_selectively starts')

    generate_count_file(samFile, refGenome, covFile, countFn, countDir, num_p=1, seperate_regions=False, sel_snp_loc=snp_locations, isRsemSam=isRsemSam, cnt_format=cnt_format)

    print('gen_count_selectively ends')

    return

#parse a count_alt_line, return gPos (0-based), lambda, and number of multimappings (excluding itself)
#
#ref:
#println2(i, count, e, of)
def check_lambda_multimapping_per_line(count_alt_line):
    #pdb.set_trace()
    tokens = count_alt_line.split()
    gPos = int(tokens[1].split(',')[0])
    lam = float(tokens[1].split(',')[1])
    num_multimapping = 0
    if len(tokens)>=3:
        alt_info = tokens[2]
        alt_info = [s[3:len(s)-1] for s in alt_info.split('[') if s != '']
        multi_locations = set()
        for a1 in alt_info:
            a1_split = a1.split(',')
            for a2 in range(len(a1_split)/2):
                gp = int(a1_split[a2*2])
                multi_locations.add(gp)
        num_multimapping = len(multi_locations)
    return gPos, lam, num_multimapping

#return res={key=gPos 0-based; val=[lambda, # mappings (excluding itself)]}
def check_lambda_multimapping(count_alt_file):
    res = {}
    with open(count_alt_file, 'r') as f:
        for line in f:
            if line[0]=='#' or line.strip()=='': continue
            [gPos, lam, num_multimapping] = check_lambda_multimapping_per_line(line)
            res[gPos]=[lam, num_multimapping]
    print('check_lambda_multimapping: %s done'%count_alt_file)
    return res

#usage: python count_read_lambda.py --test_check_count_fmt1_lam_altMap --count_file_fmt1 file
def test_check_count_fmt1_lam_altMap(args):
    count_file_fmt1 = args[args.index('--count_file_fmt1')+1]
    res = check_count_fmt1_lam_altMap(count_file_fmt1)
    pdb.set_trace()
    return

#input is a format-1 count file
#return res={key=gPos0; val=[lam, # of altMap (excluded)]}, more specifically for # of altMap (excluded),
#    -1: no snp reads
#    0:  part or all of snp reads uniq mapped
#    1,2,etc (float val): all snp reads alt mapped, # of avg alt mapping (excluding self)
def check_count_fmt1_lam_altMap(count_file_fmt1):

    #pdb.set_trace()
    res = {}
    with open(count_file_fmt1, 'r') as f:
        for line in f:
            if line[0]=='#' or line.strip()=='': continue
            [gPos, lam, num_altMap] = check_count_fmt1_lam_altMap_per_line(line)
            res[gPos]=[lam, num_altMap]
    print('check_count_fmt1_lam_altMap: %s done'%count_file_fmt1)
    #pdb.set_trace()
    return res

#parse a line of count file in format 1, 
#return gPos (0-based), lambda, and # of altMap (excluded), especially
#    -1: no snp reads
#    0:  part or all of snp reads uniq mapped
#    1,2,etc (float val): all snp reads alt mapped, # of avg alt mapping (excluding self)
def check_count_fmt1_lam_altMap_per_line(cnt_fmt1_line):

    #pdb.set_trace()
    tokens = cnt_fmt1_line.split()
    gPos = int(tokens[0])
    rB = tokens[2]
    lam = float(tokens[3])
    readTokens = tokens[8:]
    
    num_tot_reads = 0
    num_snp_reads = 0
    num_snp_reads_altMapped = 0
    acc_num_altMap = 0 #acc sum of num_altMap per snp read

    for rt in readTokens:
        rt_tB = rt.split(',')[0]
        rt_num_altMap = int(rt.split(',')[1])

        num_tot_reads += 1
        if rt_tB != rB:
            num_snp_reads += 1
            if rt_num_altMap > 0:
                num_snp_reads_altMapped += 1
                acc_num_altMap += rt_num_altMap

    if num_snp_reads == 0:
        code = -1
    elif num_snp_reads_altMapped < num_snp_reads:
        code = 0
    else:
        #pdb.set_trace()
        code = float(acc_num_altMap) / num_snp_reads_altMapped
    return gPos, lam, code

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
'''
usage:

#generate count and count_alt selectively in order to analyze snps

python count_read_lambda.py --gen_count_selectively
                            -s0 snpFile0 [-s1 snpFile1 etc]
                            --sam samFile [--isRsemSam 1/0 (default 1), ZW=0 filtered]
                            --ref refGenome
                            --cov covFile
                            --countDir countDir
                            --countFn countFn


'''

if __name__ == "__main__":

    args = sys.argv
    
    if '--gen_count_selectively' in args:

        gen_count_selectively(args)

    #usage: python count_read_lambda.py --test_check_count_fmt1_lam_altMap --count_file_fmt1 file
    elif '--test_check_count_fmt1_lam_altMap' in args:

        test_check_count_fmt1_lam_altMap(args)
        
    elif len(sys.argv)<2: #non-parallel
        
        dir_alt_map = True
    
        #working_dir = '/home/sreeramkannan/singleCell/SNP_Calling_Summer15/data_0827_SNP1k_Reads10M/'
        working_dir = '/home/sreeramkannan/singleCell/SNP_Calling_Summer15/data_0814/'
        #working_dir = '/home/olivo/Downloads/0907_count_alt_mapping/bkp/data_0814/'
        
        read_dir = "/tophat_out/"
        #read_dir = "data_GATK/1pass/"
        #read_dir = "data_GATK/2pass/"
    
        reads = working_dir + read_dir + "accepted_hits.sam"
        #reads = working_dir + read_dir + "Aligned.out.sam"
        #reads = working_dir + read_dir + "split.sam"
        #reads = working_dir + read_dir + "dedupped.sam"
        
        reference = working_dir + "Chr15.fa"
        coverage = working_dir + "coverage.txt"
        
        #count_fn = "/count_altMap_cross_check_star1pass_debug.txt"
        #count_fn = "/count_tmp_debug_0926.txt"
        count_fn = "/count_dirAltMap.txt"
        
        generate_count_file(reads, reference, coverage, count_fn)
    
        print('ready to exit count_read_lambda')
        pdb.set_trace()
    
    else:
        
        #pdb.set_trace()
        #usage:
        #python count_read_lambda.py [1]split_sam_dir [2]split_sam_addr [3]ref [4]coverage [5]count_dir [6]count_fn [7]num_p
        
        # gen count file using parallel computing
        
        dir_alt_map = True
        
        split_sam_dir = sys.argv[1]
        split_sam_addr = split_sam_dir + sys.argv[2]
        
        reference = sys.argv[3]
        coverage = sys.argv[4]
        
        count_dir = sys.argv[5]
        count_fn = sys.argv[6]
        
        num_p = int(sys.argv[7])
        seperate_regions = False
        if num_p > 1:
            seperate_regions = True
        
        print('para count_read_lambda [starts] using: %s'%sys.argv[2])
        
        #pdb.set_trace()
        generate_count_file(split_sam_addr, reference, coverage, count_fn, count_dir, num_p, seperate_regions)
    
        print('para count_read_lambda [ends  ] using: %s'%sys.argv[2])
        #pdb.set_trace()