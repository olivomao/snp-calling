import csv
import pdb
from progress.bar import Bar
from operator import itemgetter

def process(curr_grp):
    if curr_grp == []:
        #pdb.set_trace()
        return [0, 0]
        
    #pdb.set_trace()
    lid = curr_grp[0][0] #lambda id
    s = 0
    n = 0
    for itm in curr_grp:
        n = n+1
        s = s+itm[1]
    #pdb.set_trace()
    return [lid, float(s)/n]

def test():
    #working_dir = '/home/sreeramkannan/singleCell/SNP_Calling_Summer15/data_0814/'
    working_dir = '/home/sreeramkannan/singleCell/SNP_Calling_Summer15/data_0827_SNP1k_Reads10M/'
    read_dir = "/tophat_out/"
    
    #count_fn = "/count_alt_mapping_debug.txt"
    count_fn = "/count.txt"
    
    #count_fn_stat = '/count_stat_alt_mapping_debug.txt'
    count_fn_stat = '/count_stat.txt'
    
    Count_file = working_dir + read_dir + count_fn
    list_stat = []
    
    num_lines = sum(1 for line in open(Count_file))
    bar = Bar('Processing (1/3) read count...', max=num_lines)
    with open(Count_file) as counts:
            
            reader = csv.reader(counts, delimiter='\t')
            
            for row in reader:
                bar.next()
                if len(row) >= 8:
                    #pdb.set_trace()
                    Lamda_Knot = int(float(row[3])) # expression level
                    [Na,Nc,Ng,Nt] = map(int, row[4:8]) # counts
                    
                    N = Na + Nc + Ng + Nt
                    if N>0:
                        list_stat.append([Lamda_Knot, N])
                        #Count_file_stat.write('%d\t%d\n'%(Lamda_Knot, N))
    bar.finish()
    
    list_stat_sorted = sorted(list_stat, key=itemgetter(0))
    Count_file_stat = open(working_dir + read_dir + count_fn_stat, 'w+')
    
    # merge (l,n) of same l
    #pdb.set_trace()
    curr_id = -1
    curr_grp = []
    list_res = []
    
    num_lines = len(list_stat_sorted)
    bar = Bar('Processing (2/3) merge same lambda val...', max=num_lines)
    for itm in list_stat_sorted:
        bar.next()
        if itm[0] != curr_id:
            list_res.append(process(curr_grp))        
            curr_id = itm[0]
            curr_grp =[itm]
        else:
            curr_grp.append(itm)
    list_res.append(process(curr_grp))     
    bar.finish()
        
    
    num_lines = len(list_res)
    bar = Bar('Processing (3/3) write list_res...', max=num_lines)
    for itm in list_res:
        bar.next()
        Count_file_stat.write('%d\t%d\n'%(itm[0], itm[1]))
    bar.finish()
    Count_file_stat.close()
    
    print('compare_lambda0_and_reads.py: exits')   
    pdb.set_trace()
    return

if __name__ == "__main__":
    
    test()