from util import *
import sys, pdb, copy

'''
convert filtered vcf (e.g. chr15, snv only) into our target snps (e.g. SNP_m -- 1st allele & SNP_p -- 2nd allele)

usage:
python batch_run_realData.py --vcf2snp -i vcf_file -o path/to/<snp prefix>

note:
- scope:
  currently we consider GT patterns 0/1 0|1 1|0 1|1 etc (no allele 2, 3 etc)

- vcf_file format:
  #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  <sample_name>

'''
def vcf2snp(args):

    vcf_file = args[args.index('-i')+1]
    snp_prefix = args[args.index('-o')+1]
    snp_files = {} # key: 0/1 val: file path
    snp_files[0] = snp_prefix + '_m.txt'
    snp_files[1] = snp_prefix + '_p.txt'

    num_lines = sum([1 for line in open(vcf_file)])
    print('%d lines at %s'%(num_lines, vcf_file))
    i=0; j=0; T=num_lines/100;

    with open(vcf_file, 'r') as f_vcf, open(snp_files[0], 'w') as f_snp_m, open(snp_files[1], 'w') as f_snp_p:

        cnt_m = 0
        cnt_p = 0

        for line in f_vcf:

            i += 1
            if i>=T: i=0; j+=1; sys.stdout.write('\r'); sys.stdout.write('%d %% processed'%j); sys.stdout.flush()

            if line[0]=='#': continue

            tokens = line.split()

            chrom = int(tokens[0])
            pos = int(tokens[1])
            rB = (tokens[3])
            tB = (tokens[4])

            dst_st = '%d\t%s\t-->\t%s\n'%(pos, rB, tB)

            GT_info = tokens[9].split(':')[0] # 0/1, 0|1, 1|0, 1|1

            allele0 = int(GT_info[0])
            allele1 = int(GT_info[2])

            if allele1<=1 and allele0 == 1:
                cnt_m += 1
                f_snp_m.write(dst_st)

            if allele0<=1 and allele1 == 1:
                cnt_p += 1
                f_snp_p.write(dst_st)

    print('\n%s m (tot %d) & p (tot %d) files written'%(snp_prefix, cnt_m, cnt_p))

    return

'''
generate snps from ref genome, target genome, as a cross check for gen_genome_of_snps

usage:
--fa2snps -r ref_fa -c target_ch -m m_fa -p p_fa -s path/to/snp_prefix

note:
- genome m & p generated from snp m & p cross checked w/ ref genome (3/18) by test code:

        ref_fa = '/data1/shunfu1/SNPCalling/data/genome/chr15.fa'
        tgt_ch = 'chr15'
        m_fa = '/data1/shunfu1/SNPCalling/data/ref_snp/chr15_[from_snp_m].fa'
        p_fa = '/data1/shunfu1/SNPCalling/data/ref_snp/chr15_[from_snp_p].fa'
        s_pre = '/data1/shunfu1/SNPCalling/data/ref_snp/chr15_snp2'
        args = '--fa2snps -r %s -c %s -m %s -p %s -s %s'%(ref_fa, tgt_ch, m_fa, p_fa, s_pre)

        fa2snps(args.split())

'''
def fa2snps(args):

    c = args[args.index('-c')+1]

    ref_fa = args[args.index('-r')+1]
    r_genome = from_fasta(ref_fa)[c]

    m_fa = args[args.index('-m')+1]
    m_genome = from_fasta(m_fa)[c]

    p_fa = args[args.index('-p')+1]
    p_genome = from_fasta(p_fa)[c]

    snp_pre = args[args.index('-s')+1]
    snp_m = snp_pre + '_m_[fa2snps].txt'
    snp_p = snp_pre + '_p_[fa2snps].txt'

    if len(r_genome)!=len(m_genome) or len(r_genome)!=len(p_genome):
        print('fa2snps exception: unequal len of genome')
        pdb.set_trace()

    #pdb.set_trace()

    cnt_m = 0; cnt_p = 0

    with open(snp_m, 'w') as mf, open(snp_p, 'w') as pf:

        N = len(r_genome)

        for i in range(N):
            rB = r_genome[i].upper()
            mB = m_genome[i].upper()
            pB = p_genome[i].upper()

            if rB != mB:
                cnt_m += 1
                mf.write('%d\t%s\t-->\t%s\n'%(i+1, rB, mB))

            if rB != pB:
                cnt_p += 1
                pf.write('%d\t%s\t-->\t%s\n'%(i+1, rB, pB))

        print('m: %d snps and p: %d snps'%(cnt_m, cnt_p))

    return

'''
generate target genome (m or p) w/ extracted snp info

usage:
python batch_run_realData.py --gen_genome_of_snps -r ref_genome -c target_chr -s snp_file -t path/to/target_genome

note:
- gPos of snps file is 1-based

- test code:

        ref_fa = '/data1/shunfu1/SNPCalling/data/genome/chr15.fa'
        #ref_fa = '/data1/shunfu1/SNPCalling/data/ref_snp/test_ref_fa.fa' #chr15.fa'
        tgt_ch = 'chr15'
        #snp_f = '/data1/shunfu1/SNPCalling/data/ref_snp/chr15_snp_m.txt'
        snp_f = '/data1/shunfu1/SNPCalling/data/ref_snp/chr15_snp_p.txt'
        #snp_f = '/data1/shunfu1/SNPCalling/data/ref_snp/test_snp.txt' #chr15_snp_m.txt'
        #tgt_fa = '/data1/shunfu1/SNPCalling/data/ref_snp/chr15_[from_snp_m].fa'
        tgt_fa = '/data1/shunfu1/SNPCalling/data/ref_snp/chr15_[from_snp_p].fa'
        #tgt_fa = '/data1/shunfu1/SNPCalling/data/ref_snp/test_tgt_fa.fa' #chr15_[from_snp_m].fa'
        args = '-r %s -c %s -s %s -t %s'%(ref_fa, tgt_ch, snp_f, tgt_fa)
        gen_genome_of_snps(args.split())

'''
def gen_genome_of_snps(args):

    #pdb.set_trace()

    ref_genome = args[args.index('-r')+1]
    target_chr = args[args.index('-c')+1]
    snp_file = args[args.index('-s')+1]
    target_genome = args[args.index('-t')+1]

    #pdb.set_trace()

    #load snps
    snps = load_snps(snp_file) #dic key=gPos val=[rB, tB]

    #load genome
    ref_genome_seq = from_fasta(ref_genome)[target_chr]

    #add snps
    target_genome_seq = add_snps(ref_genome_seq, snps)

    #output target genome
    to_fasta(target_genome, target_genome_seq, target_chr)

    return

#snp_file: e.g. 102399211       A       -->     C
#gpos: 1 based
def load_snps(snp_file):

    snps = {}

    with open(snp_file, 'r') as sf:

        cnt = 0
        for line in sf:
            if line[0]=='#' or len(line.split())<4:
                continue
            tokens = line.split()
            gPos = int(tokens[0])
            rB = tokens[1]
            tB = tokens[3]
            if gPos in snps:
                print('load_snps exception -- existing gPos: %s'%line)
                #pdb.set_trace()
            else:
                cnt += 1
                snps[gPos] = [rB, tB]

        print('%d snps loaded'%cnt)

    #pdb.set_trace()

    return snps

'''
add snps (dic, key=gPos 1-based, val = [rB, tB]) onto ref genome seq

test code:

    rs = 'ABCDEFG'
    snps = {}
    snps[1] = ['A', 'a']
    snps[2] = ['B', 'b']
    snps[4] = ['D', 'd']
    #snps[6] = ['F','f']
    snps[7] = ['G', 'g']
    ts = add_snps(rs, snps)
    print(rs); print(ts)

'''
def add_snps(ref_genome_seq, snps):

    target_genome_seq = list(copy.copy(ref_genome_seq))

    cnt = 0

    for gPos, rB_tB in snps.items():

        cnt += 1
        if cnt >= 5000:
            sys.stdout.write('\r'); sys.stdout.write('%d processed (add snps)'%cnt); sys.stdout.flush()

        rB = rB_tB[0]; tB = rB_tB[1]

        if ref_genome_seq[gPos-1]==rB:
            #target_genome_seq[gPos] = rB_tB[1]
            #target_genome_seq = target_genome_seq[:gPos-1]+tB+target_genome_seq[gPos:]
            target_genome_seq[gPos-1]=tB
        else:
            print('add_snps exception -- rB inconsistency at gPos=%d'%gPos)
            pdb.set_trace()

    #pdb.set_trace()
    target_genome_seq = ''.join(target_genome_seq)
    print('\nadd snps done')

    return target_genome_seq