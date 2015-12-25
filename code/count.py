from itertools import *
from test import *
from read import *
import pickle
import time
import sys
import subprocess

def generate_count_file(reads, ref_address, cov_address):
    """Generate convenient count files

    Inputs:
    example sam_address: "../Data/G_30000000-1/read_l100.sam"
    example ref_address: "../Data/Chr15.fa"
    example cov_address: "../Data/G_30000000-1/coverage.txt")

    Outputs (all of these in the /output directory):
    - counts.txt
    exon_pos, reference_pos, ref_base, N_A, N_C, N_G, N_T, lambda_x, [reads]
    """
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
    
    # initialize global position to expression level map
    exon_pos = [] # exon_pos[i] = j means global position i on the genome
    exon_to_global = {}
    j = 0 # corresponds to position j on the transcriptome
    D = {}
    exon_start = []
    exon_end = []
    exon_acc_len = [0]
    with open(cov_address) as cov_file:
        for line in cov_file:
            x = line.split()
            exon_start.append(int(x[1]))
            exon_end.append(int(x[2]) )
            exon_acc_len.append(int(x[3]))
            expression = x[4]
            for i in range(int(x[1]), int(x[2])):
                D[i] = expression
                exon_pos.append(i)
                j +=1
                exon_to_global[i] = j
    G_eff = j # the total length of the exons
    print("total reference length", G)
    print("total exons length", G_eff)

    if not "sorted" in reads: # assume not sorted, so sort
        read_name = reads.split("/")[-1]
        sorted_sam = working_dir + read_dir + read_name[:-4] + "_sorted.sam"
        sort_command = "sort " + reads + " > " + sorted_sam
        print(sort_command)
        subprocess.call(sort_command, shell=True)
    else:
        print("assuming " + reads + " is already sorted")
        sorted_sam = reads


    correlations = {}

    def count(sam_address):
        """Takes the given reads at sam_address and produces a count of bases
        D is map from genome location to rna location -
            if D is provided, it is assumed to be RNA reads.
        """
        # Locus not in exon region would be mapped to counts[G_eff]
        # count= [A, C, T, G, I, D, read_map_locations]
        counts = dict()
        counts["insertions"] = dict()
        insertions = counts["insertions"]
        counter = 0 # just to report progress


        def properties(read):
            """returns (is_reversed, is_secondary)"""
            x = read.line
            flag_num = int(x[1])
            # if ((flag_num >> 7) & 1) != ((flag_num >> 6) & 1):
            #     print "0x40 and 0x80 flags are the same, something is wrong"
            #     print read.line
            return ((flag_num >> 4) & 1, (flag_num >> 6) & 1)

        def increment_by_read(read, n):
            x = read.line
            flag_num = int(x[1])
            flag = '{0:08b}'.format(flag_num)
            pos = int(x[3]) - 1
            s = num_pattern.findall(x[5])
            splice_n = [int(s[i]) for i in range(len(s))]
            splice_c = chr_pattern.findall(x[5])
            read = x[9]
            quality_scores = x[10]
            positions = [-1] * len(read)
            if flag[-3] == '0':
                read_pos = -1
                genome_pos = pos - 1
                while len(splice_n)>0:
                    for i in range(splice_n[0]):
                        read_pos += 1
                        genome_pos += 1
                        count = counts.setdefault(genome_pos, [])
                        base = read[read_pos]
                        quality = quality_scores[read_pos]
                        read_result = (base, str(n), quality)
                        count.append(read_result)
                        positions[read_pos] = genome_pos
                    splice_n = splice_n[1::]
                    splice_c = splice_c[1::]
                    while splice_n and splice_c[0] in ["N", "D", "I"]:
                        if splice_c[0]=='N':
                            genome_pos += splice_n[0]
                            splice_n = splice_n[1:]
                            splice_c = splice_c[1:]
                        if splice_c[0]=='D':
                            for _ in range(splice_n[0]):
                                genome_pos += 1
                                count = counts.setdefault(genome_pos, [])
                                count.append(('D', str(n), "I"))
                            splice_n = splice_n[1:]
                            splice_c = splice_c[1:]
                        if splice_c[0]=='I':
                            for _ in range(splice_n[0]):
                                read_pos += 1
                                count = insertions.setdefault(genome_pos, [])
                                count.append((read[read_pos].upper(), str(n), quality_scores[read_pos]))
                            splice_n = splice_n[1:]
                            splice_c = splice_c[1:]
            return positions

        with open(sam_address) as sam_file:
            sam_file.readline() # throw away the first line
            read_group = [] # temporarily stores reads
            for line in sam_file:
                if line[0] != '@': # not a comment
                    read = Read(line)
                    if not read_group or read_group[-1].id == read.id:
                        read_group.append(read)
                    elif read_group[-1].id != read.id:
                        # split into forward and backward cases
                        segment_split = [[r for r in read_group if not r.is_first_segment()], [r for r in read_group if r.is_first_segment()]]

                        for split in segment_split:
                            L = len(split)
                            groups = [[r for r in split if not r.is_reversed()], [r for r in split if r.is_reversed()]]
                            for read_group in groups:
                                positions = [increment_by_read(r, L) for r in read_group]
                                if len(positions) > 1: # deal with repeat regions
                                    for locii in zip(*positions): # locii all correlated
                                        locii = set(locii)
                                        locii.discard(-1)
                                        for m, n in combinations(locii, 2):
                                            x = correlations.setdefault(m, {})
                                            x[n] = x.get(n, 0) + 1
                                            x = correlations.setdefault(n, {})
                                            x[m] = x.get(m, 0) + 1
                        read_group = [read]

                    if counter % 10000 == 0:
                        print(counter, "lines of the .sam file are processed!",
                            round(time.clock() - start_time, 2))
                    counter += 1
            # deal with final read group after the loop
            segment_split = [[r for r in read_group if not r.is_first_segment()], [r for r in read_group if r.is_first_segment()]]

            for split in segment_split:
                L = len(split)
                groups = [[r for r in split if not r.is_reversed()], [r for r in split if r.is_reversed()]]
                for read_group in groups:
                    positions = [increment_by_read(r, L) for r in read_group]
                    if len(positions) > 1: # deal with repeat regions
                        for locii in zip(*positions): # locii all correlated
                            locii = set(locii)
                            locii.discard(-1)
                            for m, n in combinations(locii, 2):
                                x = correlations.setdefault(m, {})
                                x[n] = x.get(n, 0) + 1
                                x = correlations.setdefault(n, {})
                                x[m] = x.get(m, 0) + 1
            pickle.dump(correlations, open(working_dir + read_dir + "output/correlations.pickle", "w+"))
        print("Done processing " + sam_address)
        return counts
    
    def dump(counts, file):
        with open(file, 'w+') as outfile:
            def println(i, count, e):
                reference_base = ref[i]
                outfile.write(str(e) + "\t")
                outfile.write(str(i) + "\t") # reference position
                outfile.write(reference_base + "\t")
                outfile.write(str(round(float(D.get(i, "0")), 3))+ "\t") # expression level
                if (count):
                    N_A = sum(1 for c in count if c[0] == "A")
                    N_C = sum(1 for c in count if c[0] == "C")
                    N_G = sum(1 for c in count if c[0] == "G")
                    N_T = sum(1 for c in count if c[0] == "T")
                    outfile.write(str(N_A) + "\t")
                    outfile.write(str(N_C) + "\t")
                    outfile.write(str(N_G) + "\t")
                    outfile.write(str(N_T) + "\t")
                    outfile.write("\t".join([str(i) for i in map(",".join, count)]))
                if j in correlations:
                    correlated = correlations[j]
                    outfile.write("\t|")
                    for q in correlated:
                        outfile.write("\t")
                        k = correlated[q]
                        outfile.write(",".join(map(str, [q, k, round(float(D.get(q, "0")), 3), ref[q]])))
                outfile.write("\n")
            for i, j in enumerate(exon_pos):
                x = counts.get(j, [])
                println(j, x, i)
        print(file + " written")
        # with open('insertions' + file, 'w+') as outfile:
        #     insertions = counts["insertions"]
        #     for pos in insertions.iterkeys():
        #         outfile.write(str(pos) + "\t")
        #         outfile.write(str(round(float(D.get(pos, "0")), 3)))
        #         outfile.write("\t")
        #         outfile.write("\t".join([str(i) for i in map(",".join, insertions[pos])]))
        #         outfile.write("\n")
        # print('insertions' + file + " written")

    counts = count(sorted_sam)
    dump(counts, working_dir + read_dir + "output/debug_count.txt")
    del counts

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
    working_dir = "/Users/Soheil/Dropbox/Transcriptome/SNP-Calling-Summer15/data_2/"
    read_dir = "tophat_out_m/"
    reads = working_dir + read_dir + "accepted_hits.sam"        #"data/G_30000000-1/read_l100_sorted.sam"
    reference = working_dir + "Chr15.fa"
        #"data/G_30000000-1/Chr15.fa"
    coverage = working_dir + "coverage.txt"
        #"data/G_30000000-1/coverage.txt"
    generate_count_file(reads, reference, coverage)
