import sys, pdb
from intervaltree import Interval, IntervalTree 

from tool_address import FA2FQ_address

from old_code.batch_run_case1plus6 import enforce_unique_readID
from old_code.util import run_cmd
from old_code.snp_oper import load_snps, gen_genome_of_snps
from old_code.Synthesis import calc_exp_sum, ExpressionLevel2Coverage
from old_code.ReadProcess import BED2Exon1

from old_code_multiShannon.sim_reads import sim_reads_g_b

def main():

    pdb.set_trace()
    
    return

'''
usage:

python sim_data_generator.py --mark_homozygous -m snp_m file -p snp_p file
                                               [-m2 snp_m2_file] [-p2 snp_p2_file]
                                               [--no_homozygous]
                                               [--outputHomo snp_homo_file]

description:
- if --no_homozygous:
  homozygous snps will not be output
  from snp_m.txt and snp_p.txt, output snp_m_noHomo.txt and snp_p_noHomo.txt
- if no --no_homozygous (default)
  homozygous will be marked as: gPos rB --> tB (*)
  from snp_m.txt and snp_p.txt, output snp_m_markHomo.txt and snp_p_markHomo.txt
- if --outputHomo is specified, we output homozygous snps into specified snp_homo_file

'''
def mark_homozygous(args):

    if '--no_homozygous' in args:
        noHomo = True
    else:
        noHomo = False

    m1 = args[args.index('-m')+1]
    if '-m2' in args:
        m2 = args[args.index('-m2')+1]
    else:
        if noHomo == True:
            m2 = m1[:-4]+'_noHomo.txt'
        else:
            m2 = m1[:-4]+'_markHomo.txt'

    p1 = args[args.index('-p')+1]
    if '-p2' in args:
        p2 = args[args.index('-p2')+1]
    else:
        if noHomo == True:
            p2 = p1[:-4]+'_noHomo.txt'
        else:
            p2 = p1[:-4]+'_markHomo.txt'

    snps_m = load_snps(m1)
    snps_p = load_snps(p1)

    with open(m2, 'w') as of:
        itms = sorted(snps_m.items(), key=lambda x:x[0])
        for k, v in itms:
            if k in snps_p:
                if noHomo == False:
                    st = '%d\t%s\t-->\t%s\t(*)\n'%(k, v[0], v[1])
                    of.write(st)
            else:
                st = '%d\t%s\t-->\t%s\n'%(k, v[0], v[1])
                of.write(st)
    print('%s written'%m2)

    with open(p2, 'w') as of:
        itms = sorted(snps_p.items(), key=lambda x:x[0])
        for k, v in itms:
            if k in snps_m:
                if noHomo == False:
                    st = '%d\t%s\t-->\t%s\t(*)\n'%(k, v[0], v[1])
                    of.write(st)
            else:
                st = '%d\t%s\t-->\t%s\n'%(k, v[0], v[1])                
                of.write(st)
    print('%s written'%p2)

    if '--outputHomo' in args:
        #pdb.set_trace()
        snp_homo_file = args[args.index('--outputHomo')+1]

        with open(snp_homo_file, 'w') as shf:
            itms = sorted(snps_m.items(), key=lambda x:x[0])
            for k, v in itms:
                if k in snps_p:
                    #pdb.set_trace()
                    st = '%d\t%s\t-->\t%s\n'%(k, v[0], v[1])
                    shf.write(st)

        print('%s written'%snp_homo_file)

    return

'''
usage:

python sim_data_generator.py --sel_snps
                             -i input_snp_file
                             -o output_snp_file
                             -e region_file

description:

filter input_snp_file into output_snp_file, such that:
- snps are restricted by region_file

'''
def sel_snps(args):

    input_snp_file = args[args.index('-i')+1]
    output_snp_file = args[args.index('-o')+1]
    region_file = args[args.index('-e')+1]

    tree = IntervalTree()

    with open(region_file, 'r') as rf:

        #pdb.set_trace()

        for line in rf:
            if line[0]=='#': continue
            tokens = line.split()
            if len(tokens)>=2:
                tree.add(Interval(int(tokens[0]), int(tokens[1])))
    #pdb.set_trace()

    with open(input_snp_file, 'r') as f1, open(output_snp_file, 'w') as f2:

        for line in f1:
            if line[0]=='#': continue
            gPos = int(line.split()[0])-1 #1-based into 0-based
            query_res = tree[gPos] #a set of Intervals containing gPos
            if len(query_res)!=0:
                #pdb.set_trace()
                query_res_int = sorted(query_res)[0]
                #st = '%s\t[%d,%d)\n'%(line.strip(), query_res_int.begin, query_res_int.end)
                st = line
                f2.write(st)

    return

'''
generate line coverage (i.e. abundance per exonic location)

usage:
python sim_data_generator.py --gen_coverage 
                             -b sorted_bed -e exp_file -o coverage_file
                             -L readLength -N numReads

note:
we make readLength, numReads and exp_sum more flexible here, instead of fixing them to be 10, 1e6 and 10 as before,
which may be the reason why prev ideal coverage is worse than rsem coverage

'''
def gen_coverage(args):

    #pdb.set_trace()

    sorted_bed = args[args.index('-b')+1]
    exp_file = args[args.index('-e')+1]
    coverage_file = args[args.index('-o')+1]
    cov_fn = '/'+coverage_file.split('/')[-1]

    readLength = int(args[args.index('-L')+1])
    numReads = int(args[args.index('-N')+1])    
    exp_sum = calc_exp_sum(exp_file)

    ExpressionLevel2Coverage(BED_sorted_address=sorted_bed, 
                             exp_address=exp_file,
                             cov_fn=cov_fn, Stat=None,
                             Lr=readLength, tot_N=numReads, exp_sum=exp_sum)
    return

'''
usage:
python sim_data_generator.py --snp_cov -s snp_file -m cov_file_m -p cov_file_p -o snp_cov_file

description:
- find per snp's exon interval and related abundance, helpful to see how snp is covered
- snp_cov_file: gPos(1-based), exon_stt(0-based), exon_stp(excluded), abundance of cov_file_m, abundance of cov_file_p

'''
def snp_cov(args):

    snp_file = args[args.index('-s')+1]
    cov_file_m = args[args.index('-m')+1]
    cov_file_p = args[args.index('-p')+1]
    snp_cov_file = args[args.index('-o')+1]

    tree = IntervalTree()

    #pdb.set_trace()
    with open(cov_file_m, 'r') as cf_m, open(cov_file_p, 'r') as cf_p:

        for line, line2 in zip(cf_m, cf_p):
            if line[0]=='#': continue
            tokens = line.split()
            tokens2 = line2.split()
            if len(tokens)>=5:
                tree.add(Interval(int(tokens[1]), int(tokens[2]), [float(tokens[4]), float(tokens2[4])]))
    
    #pdb.set_trace()
    snps = load_snps(snp_file)
    snp_cov = []
    for k, v in snps.items():
        gPos = k-1 
        query_res = tree[gPos]
        if len(query_res)!=0:
            query_res_int =  sorted(query_res)[0]
            snp_cov.append([k, v[0], v[1], query_res_int.begin, query_res_int.end, query_res_int.data[0], query_res_int.data[1]])
        else:
            snp_cov.append([k, v[0], v[1], -1, -1, 0, 0])

    #pdb.set_trace()
    snp_cov = sorted(snp_cov, key=lambda x:x[0])
    with open(snp_cov_file, 'w') as of:
        st = '#gPos\texonStart\texonStop\tmCoverage\tpCoverage\n'
        of.write(st)
        for k, rB, tB, exon_stt, exon_stp, abundance_m, abundance_p in snp_cov:
            st = '%d\t%s\t-->\t%s\t%d\t%d\t%f\t%f\n'%(k, rB, tB, exon_stt, exon_stp, abundance_m, abundance_p)
            of.write(st)
    #pdb.set_trace()

    print('%s written'%snp_cov_file)
    return

'''
find snp and whether it's covered by reads

usage:
python sim_data_generator.py --snp_read_cov -s snp_file -m m_read_bed -p p_read_bed -o snp_read_cov_file

format of snp_read_cov_file:
col-0: gPos (1-based)
col-1: # of m reads covering the snp
col-2: # of p reads covering the snp
'''
'''
def snp_read_cov(args):

    snp_file = args[args.index('-s')+1]
    m_read_bed = args[args.index('-m')+1]
    p_read_bed = args[args.index('-p')+1]
    read_beds = [m_read_bed, p_read_bed]
    trees = [IntervalTree(), IntervalTree()]

    snp_read_cov_file = args[args.index('-o')+1]

    snps = load_snps(snp_file)

    intermediate_files = []

    #build trees of read bed files for query purpose    
    for i in range(2):
        tmp_exon_file = read_beds[i] + '.tmp_exon_file'; intermediate_files.append(tmp_exon_file)
        BED2Exon1(read_beds[i], tmp_exon_file)
        print('%s written'%tmp_exon_file)
        #pdb.set_trace()

        tree = trees[i]
        with open(tmp_exon_file, 'r') as tef:
            for line in tef:
                if line[0]=='#' or len(line.split())<3: continue
                tokens = line.split()
                e_stt = int(tokens[0]) #0-based
                e_stp = int(tokens[1]) #exclusive
                data = [] #read ids
                for t in xrange(2,len(tokens)):
                    rid = tokens[t].split(',')[0]
                    data.append(rid)
                tree.add(Interval(e_stt, e_stp, data))
        print('%d-th tree built'%i)
        pdb.set_trace()

    #check snp read coverage and write to output
    snp_cnt = {} #key gPos (1-based), val=[# of m reads covering the snp, # of p reads covering the snp]
    for k, v in snps.items():
        gPos = k-1
        mp_sum = [0,0]
        for j in range(2):
            query_res = trees[j][gPos]
            if len(query_res) != 0:
                for interval in sorted(query_res):
                    mp_sum[j] += len(interval.data)
        snp_cnt[k]=mp_sum
    #pdb.set_trace()
    snp_cnt = snp_cnt.items()
    snp_cnt = sorted(snp_cnt, key=lambda x:x[0])
    pdb.set_trace()


    with open(snp_read_cov_file, 'w') as of:

        st = '#gPos\t#_m_reads\t#_t_reads\n'
        of.write(st)

        for gPos, mn_cnts in snp_cnt:
        
            st = '%d\t%d\t%d\n'%(gPos, mn_cnts[0], mn_cnts[1])
            of.write(st)
        #pdb.set_trace()
        print('%s written'%snp_read_cov_file)

    #clear intermediate files
    for i in range(len(intermediate_files)):
        cmd = 'rm %s'%intermediate_files[i]
        run_cmd(cmd)
    #pdb.set_trace()

    return
'''

def bed2cov(bed):
    res = {} #gPos (0-based): coverage

    with open(bed, 'r') as inf:

        nLines=sum([1 for l in open(bed,'r')]); T=nLines/100; p=0; q=0;

        for line in inf:
            p += 1
            if p>=T: p=0; q+=1; sys.stdout.write('\r'); sys.stdout.write('%d %% processed (bed2cov)'%q); sys.stdout.flush()
            if line[0]=='#' or len(line.split())<4: continue

            x = line.split()
            tr_ID = x[3]
            tr_start = int( x[1] )
            number_exon = int( x[9] )
            EXON_len_strng = x[10].split(',') 
            EXON_len = [ int(EXON_len_strng[i]) for i in range(number_exon) ]
            EXON_start_strng = x[11].split(',')
            EXON_start = [ int(EXON_start_strng[i]) for i in range(number_exon) ]
        
            for i in range(number_exon):
                exon_start = tr_start + int(EXON_start[i]) #0-based,
                exon_end = exon_start + int(EXON_len[i]) #exclusive
                for j in xrange(exon_start, exon_end):
                    if j in res:
                        res[j]+=1
                    else:
                        res[j]=1
    return res

'''
find snp and whether it's covered by reads

usage:
python sim_data_generator.py --snp_read_cov1 -s snp_file -m m_read_bed -p p_read_bed -o snp_read_cov_file

format of snp_read_cov_file:
col-0: gPos (1-based)
col-1: # of m reads covering the snp
col-2: # of p reads covering the snp
col-3: (*) if no m/p reads covering the snps
'''
def snp_read_cov1(args): #try to avoid BED2Exon1

    snp_file = args[args.index('-s')+1]
    m_read_bed = args[args.index('-m')+1]
    p_read_bed = args[args.index('-p')+1]
    snp_read_cov_file = args[args.index('-o')+1]

    snps = load_snps(snp_file)

    #bed --> {gPos:cov}
    gPos_cov_m = bed2cov(m_read_bed)
    gPos_cov_p = bed2cov(p_read_bed)

    #pdb.set_trace()
    snp_cnt = {}
    for k, v in snps.items():
        gPos = k-1
        mp_sum = [0,0]
        if gPos in gPos_cov_m:
            mp_sum[0]=gPos_cov_m[gPos]
        if gPos in gPos_cov_p:
            mp_sum[1]=gPos_cov_p[gPos]
        snp_cnt[k]=mp_sum
    snp_cnt = snp_cnt.items()
    snp_cnt = sorted(snp_cnt, key=lambda x:x[0])
    #pdb.set_trace()

    num_no_coverage = 0
    with open(snp_read_cov_file, 'w') as of:

        st = '#gPos\t#_m_reads\t#_t_reads\t(no coverage)\n'
        of.write(st)

        for gPos, mn_cnts in snp_cnt:
        
            if sum(mn_cnts)==0:
                num_no_coverage += 1
                st = '%d\t%d\t%d\t(*)\n'%(gPos, mn_cnts[0], mn_cnts[1])
            else:
                st = '%d\t%d\t%d\n'%(gPos, mn_cnts[0], mn_cnts[1])
            of.write(st)
        #pdb.set_trace()
        print('%d out of %d snps not covered'%(num_no_coverage, len(snp_cnt)))
        print('%s written'%snp_read_cov_file)

    return


'''
usage:
python sim_data_generator.py --reasign_exp 
                             -s snp_file
                             -b sorted_bed_file
                             -e old_exp_file
                             -o new_exp_file

description:

- we collect exp vals from old_exp_file, and reasign them from high to low to the transcripts that contain highest snps to lowest snps
- for example, for snps at data_real/ref_snp/ in chr15 coding region w/o homozygous snps, 
               and based on bed file gencode.v19.annotation.target.chr15.sorted.bed,
               and generated exp files at data_real_abSNP_batch/reads_N100000_L100_Err0.00/intermediate/,
  we find for m allele, there're 1001 heterozygous snps. among 7314 trids, about 1200 (16.4%) trids contain at least 1 snp. we assign top 16.4% highest exp vals to these trids.
          for p allele, there're 899  heterozygous snps, among 7314 trids, about 1100 (15.0%) trids contain at least 1 snp. we assign top 15.0% highest exp vals to these trids.
          for m & p alleles, there're 1462 homozygous snps, among 7314 trids, about 1700 (23%) trids contain at least 1 snp. we assign top 23% highest exp vals to these trids.
'''
def reasign_exp(args): 

    snp_file = args[args.index('-s')+1]
    sorted_bed_file = args[args.index('-b')+1]
    old_exp_file = args[args.index('-e')+1]
    new_exp_file = args[args.index('-o')+1]

    intermediate_files = []
    tmp_exon_file = sorted_bed_file + '.tmp_exon_file'; intermediate_files.append(tmp_exon_file)
    BED2Exon1(sorted_bed_file, tmp_exon_file)
    #pdb.set_trace()

    #interval tree, each interval=exon, corresponding data=list of trs
    tree = IntervalTree()
    with open(tmp_exon_file, 'r') as ef:
        for line in ef:
            if line[0]=='#' or len(line.split())<3: continue
            tokens = line.split()
            e_stt = int(tokens[0]) #0-based
            e_stp = int(tokens[1]) #exclusive
            data = []
            for t in xrange(2,len(tokens)):
                trid = tokens[t].split(',')[0]
                data.append(trid)
            tree.add(Interval(e_stt, e_stp, data))
    #pdb.set_trace()

    #tr_snp_cnt = {} #key - trid; val - # of snps
    tr_snp_cnt = {}
    with open(sorted_bed_file, 'r') as bf:
        for line in bf:
            if line[0] == '#' or len(line.split())<4: continue
            tokens = line.split()
            trid = tokens[3]
            tr_snp_cnt[trid] = 0
    #pdb.set_trace()

    #mark transcripts that contain snps
    snps = load_snps(snp_file)
    for gPos, rB_tB in snps.items():
        gPos = gPos-1 
        query_res = tree[gPos]
        if len(query_res)!=0:
            intervals =  sorted(query_res)
            for interval in intervals:
                interval_trids = interval.data 
                for trid in interval_trids:
                    if trid in tr_snp_cnt:
                        tr_snp_cnt[trid] += 1
    tr_snp_cnt = tr_snp_cnt.items()
    tr_snp_cnt = sorted(tr_snp_cnt, key=lambda x:-x[1])
    #pdb.set_trace()

    #re-asign exp values
    exp_vals = []
    with open(old_exp_file, 'r') as oef:
        for line in oef:
            if line[0]=='#' or len(line)<8: continue
            exp_vals.append(float(line.split()[7]))
    exp_vals = sorted(exp_vals, key=lambda x:-x)
    #pdb.set_trace()

    tr_exp = {} #trid, new exp
    for i in range(len(tr_snp_cnt)):
        tr_exp[tr_snp_cnt[i][0]] = exp_vals[i]
    #pdb.set_trace()

    with open(old_exp_file, 'r') as oef, open(new_exp_file, 'w') as nef:
        for line in oef:
            if line[0]=='#' or len(line)<8:
                nef.write(line)
            else:
                tokens = line.split()
                tokens[7] = str(tr_exp[tokens[0]])
                st = '\t'.join(tokens)+'\n'
                nef.write(st)

    #clear intermediate files
    for i in range(len(intermediate_files)):
        cmd = 'rm %s'%intermediate_files[i]
        run_cmd(cmd)
    #pdb.set_trace()

    return

'''
example:

python sim_data_generator.py --read_generation -g /data1/shunfu1/SNPCalling/snp20_reads100k_10/Tar_m.txt 
                             -b /data1/shunfu1/SNPCalling/snp20_reads100k_10/hg19_chr15-UCSC-sorted.bed
                             -O /data1/shunfu1/SNPCalling/snp20_reads100k_10_abSNPcode/reads_N100k_L100_Err0 
                             -n 100000 -l 100 -r 0 
                             --suffix _m --toFq
'''
def read_generation(args):

    sim_reads_g_b(args)

    #add suffix
    out_dir = args[args.index('-O')+1]

    suf = args[args.index('--suffix')+1]      
    
    src_dst_pairs = []

    reads = '%s/reads.fa'%(out_dir); readsNew = '%s/reads%s.fa'%(out_dir, suf)
    src_dst_pairs.append([reads, readsNew])
    
    if '--exp_path' not in args:
        exp = '%s/intermediate/exp.txt'%(out_dir); expNew = '%s/intermediate/exp%s.txt'%(out_dir, suf)
        src_dst_pairs.append([exp, expNew])

    reads_bed = '%s/intermediate/reads.bed'%(out_dir); reads_bedNew = '%s/intermediate/reads%s.bed'%(out_dir, suf)
    src_dst_pairs.append([reads_bed, reads_bedNew])

    for i in range(len(src_dst_pairs)):
        run_cmd('mv %s %s'%(src_dst_pairs[i][0], src_dst_pairs[i][1]))

    #fa to fq
    readsFq = '%s/reads%s.fq'%(out_dir, suf)
    run_cmd( 'perl ' + FA2FQ_address + ' ' + readsNew + ' > '  + readsFq )

    return

def merge_reads(args):

    reads_m = args[args.index('-m')+1]
    reads_p = args[args.index('-p')+1]
    merged_reads = args[args.index('-o')+1]

    run_cmd('cat '  + reads_m + ' ' + reads_p + ' > ' + merged_reads)

    if '--uniqID' in args:
        enforce_unique_readID(merged_reads, flag=True)
        
    return

'''
usage:

#filter input_snp_file into output_snp_file, such that:
# - snps are restricted by region_file

python sim_data_generator.py --sel_snps
                             -i input_snp_file
                             -o output_snp_file
                             -e region_file

#mark or filter snps
#- if --no_homozygous:
#  homozygous snps will not be output
#  from snp_m.txt and snp_p.txt, output snp_m_noHomo.txt and snp_p_noHomo.txt
#- if no --no_homozygous (default)
#  homozygous will be marked as: gPos rB --> tB (*)
#  from snp_m.txt and snp_p.txt, output snp_m_markHomo.txt and snp_p_markHomo.txt

python sim_data_generator.py --mark_homozygous -m snp_m file -p snp_p file
                                               [-m2 snp_m2_file] [-p2 snp_p2_file]
                                               [--no_homozygous]
                                               [--outputHomo snp_homo_file]


# we collect exp vals from old_exp_file, and reasign them from high to low to the transcripts that contain highest snps to lowest snps

python sim_data_generator.py --reasign_exp 
                             -s snp_file
                             -b sorted_bed_file
                             -e old_exp_file
                             -o new_exp_file

#generate target genome (m or p) w/ extracted snp info

python sim_data_generator.py --gen_genome_of_snps -r ref_genome -c target_chr -s snp_file -t path/to/target_genome


# use RNASeqReadSimulator to generate reads per target allele
#
# we re-use the simulator written in multi shannon
# for snp calling purpose, we need to add suffix (e.g. _m, _p) and convert fa to fq
#
# note:
#
# .fa, \.sorted.bed
#                               --> out_dir/intermediate/exp[_suffix].txt (if --exp_path unspecified), reads[_suffix].bed 
#                               --> out_dir/reads[_suffix].[fa or fq] (SE)
# 
# .sorted.bed needs 1st line to be "# trNum (trNum) avgTrLen (avgTrLen)" for '-c readCoverage'
#

python sim_data_generator.py --read_generation -g genomeFile -b trBedFile -O out_dir
                             (-n numReads or -c readCoverage) -l readLength -r errRate
                             --suffix suf --toFq
                             [--exp_path exp_path]

# pool reads sampled from maternal and paternal alleles together
#

python sim_data_generator.py --merge_reads -m reads_m.fq -p reads_p.fq -o merged_reads.fq [--uniqID]

'''
if __name__ == "__main__":

    args = sys.argv

    # select snps (true snps --> restrict to snps of interest e.g. coding regions)
    if '--sel_snps' in args:
        sel_snps(args)
    #from snp_m.txt and snp_p.txt, output snp_m_markHomo.txt and snp_p_markHomo.txt
    elif '--mark_homozygous' in args: 
        mark_homozygous(args)
    #generate coverage file
    elif '--gen_coverage' in args:
        gen_coverage(args)
    #find snp and its coverage
    elif '--snp_cov' in args:
        snp_cov(args)
    #find snp and whether it's covered by reads -- not efficient, replaced by snp_read_cov1
    #elif '--snp_read_cov' in args:
    #    snp_read_cov(args)
    elif '--snp_read_cov1' in args:
        snp_read_cov1(args)
    #re-asign exp file so that snp contained transcript gets high expression (e.g. check multiple alignments)
    elif '--reasign_exp' in args:
        reasign_exp(args)
    #generate target individual genome
    elif '--gen_genome_of_snps' in args:
        gen_genome_of_snps(args)
    # generate simulated reads
    elif '--read_generation' in args:
        read_generation(args)
    elif '--merge_reads' in args:
        merge_reads(args)
    else:
        main()