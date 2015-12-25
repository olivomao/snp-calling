from itertools import izip
from test import *
from read import *
import time
from operator import itemgetter
import sys

def ExpectedCoverage(exp_address, N):
    # Reads Per Kilobases transcript per Millions of reads
    # RPKM(i) = C8(i) * 10^9 / ( C2(i) * sum_j C8(j) )
    # N = total number of reads
    # Expected number of reads from transcript i with length Li = RPKM(i) * (N/10^6) * (Li/10^3)
    #     = N * C8(i) / sum_j C8(j)
    line_counter=0;
    Explevel=[]
    Length = []
    with open(exp_address) as exp_file:
        for line in exp_file:
            if line_counter==0:
                line_counter = 1
            else:
                line_splitted = line.split()
                Explevel.append( float( line_splitted[7] ) ) #7 #10
                Length.append( int(line_splitted[1] ) ) #1 #2
                
    Coverage = [ 0 for i in range(len(Explevel)) ]
    Explevel_sum = sum(Explevel)
    for i in range(0,len(Explevel)):
        Coverage[i] = float(N) * Explevel[i] / Explevel_sum
        
    return Explevel_sum

def ExpressionLevel2Coverage(BED_address, exp_address, N, L, file_type = "rsem"):
    #b=re.search('([\S]+)/([\S]+)\.([\S]+)', BED_address)
    #coverage_address = b.group(1) + '/coverage.txt'
    coverage_address = "cuff_coverage.txt"
    
    Explevel_file = open(exp_address, 'rU')
    Explevel_file.readline() # because the first line is just header
    BED_file = open(BED_address, 'rU')
    
    Total_number_segments=0

    vector = []
    
    if file_type == "rsem":
        Explevel_sum = ExpectedCoverage(exp_address, N)
        for B_line, E_line in zip(BED_file,Explevel_file):
            x = B_line.split()
            tr_start = int( x[1] )
            number_exon = int( x[9] )
            exon_len = x[10].split(',')
            exon_start = x[11].split(',')
            
            
            y = E_line.split()
            explevel = float( y[7] ) #7 #10
            transcript_len = int (y[1]) #1 #2
            
            line_cover = L * float(N) * explevel / (Explevel_sum * float(transcript_len))
            for i in range(number_exon):
                vector.append([tr_start + int(  exon_start[i] ) ,  Total_number_segments + i, line_cover, -1] )
                vector.append([tr_start + int(  exon_start[i] ) + int(exon_len[i]) , Total_number_segments + i, line_cover, 1] )
            Total_number_segments += number_exon
            
    if file_type == "cuff":
        bed_info = {}
        for B_line in BED_file:
            x = B_line.split()
            bed_info[x[3]] = x

        for E_line in Explevel_file:
            y = E_line.split()
            explevel = float( y[9] )
            transcript_len = int (y[7])
            bed_id = y[0]

            x = bed_info[bed_id]
            tr_start = int( x[1] )
            number_exon = int( x[9] )
            exon_len = x[10].split(',')
            exon_start = x[11].split(',')
            
            line_cover = L * float(N) * explevel / (1 * float(transcript_len))
            for i in range(number_exon):
                vector.append([tr_start + int(  exon_start[i] ) ,  Total_number_segments + i, line_cover, -1] )
                vector.append([tr_start + int(  exon_start[i] ) + int(exon_len[i]) , Total_number_segments + i, line_cover, 1] )
            Total_number_segments += number_exon

    vector_sorted = sorted(vector, key=itemgetter(0,3) )
    vector_sorted_updated = [[vector_sorted[0][0], vector_sorted[0][2] ]]
    index = 0
    for i in range(1,len(vector_sorted)):
        if vector_sorted[i][0]> vector_sorted_updated[index][0]:
            vector_sorted_updated.append([ vector_sorted[i][0], vector_sorted_updated[index][1] - vector_sorted[i][3]* vector_sorted[i][2]] )
            index += 1
        else:
            vector_sorted_updated[index][1] -= vector_sorted[i][3]* vector_sorted[i][2]


    counter = 0
    Accumul_length = 0
    Cov_file = open( coverage_address, 'w+')
    for i in range(len(vector_sorted_updated)):
        if abs(vector_sorted_updated[i][1])>0.0000001:
            counter +=1
            Accumul_length +=vector_sorted_updated[i+1][0] - vector_sorted_updated[i][0]
            Cov_file.write( repr(counter) + '\t' + repr(vector_sorted_updated[i][0]) + '\t' + repr(vector_sorted_updated[i+1][0]) + '\t' + repr(Accumul_length) + '\t' + repr(round(vector_sorted_updated[i][1], 3)) +'\n')   
         
    Cov_file.close()
    return coverage_address

bed_file = "13-4_transcripts-chr15.bed"
exp_file = "13-4_isoforms-chr15.fpkm_tracking"
N = sys.argv[3] # 1
L = sys.argv[4] # 1
file_type = sys.argv[5] # "cuff" or "rsem"
print(ExpressionLevel2Coverage(bed_file, exp_file, float(N), float(L), file_type))