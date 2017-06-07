import pdb, sys, os
from old_code.util import run_cmd, run_cmds, from_fasta
from intervaltree import Interval, IntervalTree

'''
usage:

python test_varscan.py --run_varscan
                       [--vdir vdir]
                       [--mode mode]
                       [--filter]
                       --refFa refFa
                       --bamFile bamFile
                       --resFile resFile

description:
run varscan and store res in resFile

--mode: 0,1 ;default 0; 1: tune param down (max sens)
--filter: filter false positives by varscan using default parameters
'''
def run_varscan(args):

    '''    
    #ref_file = '/data1/shunfu1/SNPCalling//data_large_0_idealCov/Chr15.fa' 
    ref_file = '/data1/shunfu1/SNPCalling/SimSNPs_MultiRun_5/reads_N100K_L100_Err0.00/rsem/Chr15.transcripts.fa'
    #bam_file = '/data1/shunfu1/SNPCalling/SimSNPs_MultiRun_5/reads_N100K_L100_Err0.00/gatk/2pass/dedupped.bam'
    #bam_file = '/data1/shunfu1/SNPCalling/SimSNPs_MultiRun_5/reads_N100K_L100_Err0.00/rsem/Chr15.genome.sorted.bam'
    bam_file = '/data1/shunfu1/SNPCalling/SimSNPs_MultiRun_5/reads_N100K_L100_Err0.00/rsem/Chr15.transcript.sorted.bam'

    res_file_name = 'run_5_n100k_err000_rsem_transcriptome_res.txt'
    '''

    if '--vdir' in args:
        varscan_dir = args[args.index('--vdir')+1]
    else:
        varscan_dir = '/home/shunfu1/software/VarScan/'

    ref_file = args[args.index('--refFa')+1]
    bam_file = args[args.index('--bamFile')+1]
    res_file = args[args.index('--resFile')+1]

    if '--mode' in args:
        mode = int(args[args.index('--mode')+1])
    else:
        mode = 0

    if mode == 0:
        mode_str = ''
    elif mode == 1:
        mode_str = ' --min-coverage 1 --min-reads2 1 --min-avg-qual 1 '
    else:
        print('unknown mode for varscan'); pdb.set_trace()

    if '--filter' in args:
        cmd = 'samtools mpileup -f %s %s | java -jar %s/VarScan.v2.3.9.jar pileup2snp %s | java -jar %s/VarScan.v2.3.9.jar filter --output-file %s'%\
          (ref_file, bam_file, varscan_dir, mode_str, varscan_dir, res_file)
    else:
        cmd = 'samtools mpileup -f %s %s | java -jar %s/VarScan.v2.3.9.jar pileup2snp %s > %s'%\
          (ref_file, bam_file, varscan_dir, mode_str, res_file)
    #pdb.set_trace()
    run_cmd(cmd)

    return

'''
usage:
python test_varscan.py --run_varscan_batch

description:
run varscan for a batch of data/modes
'''
def run_varscan_batch(args):

    #configuration parameters
    reads = [] #list of [readLen, read err, read number, read number description]
    #reads.append([100, 0.00, 100000, '100K'])
    #reads.append([100, 0.01, 100000, '100K'])
    #reads.append([100, 0.00, 1000000, '1m'])
    reads.append([100, 0.01, 1000000, '1m'])
    #reads.append([100, 0.00, 10000000, '10m'])
    #reads.append([100, 0.01, 10000000, '10m'])
    
    N_run_stt = 5
    N_run_stp = 9 #inclusive

    modes = [0, 1] #[0,1]

    #filter_str = ''; filter_nm = ''
    filter_str = '--filter'; filter_nm = '_filt'

    sam2bam = 0
    run_varscan = 1
    convert_varscan = 1

    num_p = 10
    if num_p > 1:
        sam2bam_cmds = []
        run_varscan_cmds = []
        convert_varscan_cmds = []
    
    #default parameters
    SrcDir = '/data1/shunfu1/SNPCalling/data_large_0_idealCov/'
    RootFolder = '/data1/shunfu1/SNPCalling/'

    refGenome = '%s/Chr15.fa'%SrcDir
    bedSorted = '%s/hg19_chr15-UCSC-sorted.bed'%SrcDir

    #auto configuration
    #RootFolder/SimSNPs_MultiRun_<i>/expression & snps & tar genome files & cov files
    #                               /<readsLabel>/
    #                                            /intermediate/ read bed files
    #                                            /read files
    #                                            /rsem/ quantification files
    #                                            /gatk/ gatk files
    #                                                 /2pass/Aligned.out.sorted.bam (*)
    #                                            /snp_res/ snp res files
    #                                            /snp_res_gatk/ snp res (vs gatk) files
    #                                            /snp_res_varscan/ snp res (vs varscan) files: snp_varscan_mode<0 or 1>_varscanFormat.txt snp_varscan_mode<0 or 1>.txt

    for ith_run in xrange(N_run_stt, N_run_stp+1):

        #pdb.set_trace()

        RunFolder = '%s/SimSNPs_MultiRun_%d/'%(RootFolder, ith_run)        
        
        for read_config in reads:
            readLen = read_config[0]
            readErr = read_config[1]
            readNum = read_config[2]
            readNumDescription = read_config[3]

            readsLabel = 'reads_N%s_L%s_Err%.2f'%(readNumDescription, readLen, readErr) 

            bamFile = '%s/%s/gatk/2pass/Aligned.out.sorted.bam'%(RunFolder, readsLabel)
            if os.path.exists(bamFile)==False:
                samFile = '%s/%s/gatk/2pass/Aligned.out.sam'%(RunFolder, readsLabel)
                cmd = 'samtools sort -@ 20 -o %s %s'%(bamFile, samFile)
                if sam2bam==1:
                    if num_p==1:
                        run_cmd(cmd)
                    else:
                        sam2bam_cmds.append(cmd)

            resFolder = '%s/%s/snp_res_varscan/'%(RunFolder, readsLabel)
            cmd = 'mkdir -p %s'%resFolder
            run_cmd(cmd)

            #run varscan
            for mode in modes:
                #pdb.set_trace()
                resFile_a = '%s/snp_varscan_mode%d_varscanFormat%s.txt'%(resFolder, mode, filter_nm)
                resFile_b = '%s/snp_varscan_mode%d%s.txt'%(resFolder, mode, filter_nm) 

                cmd = 'python test_varscan.py --run_varscan ' +\
                                              '--mode %d  %s '%(mode, filter_str) +\
                                              '--refFa %s '%refGenome +\
                                              '--bamFile %s '%bamFile +\
                                              '--resFile %s '%resFile_a
                if run_varscan==1:
                    if num_p==1:
                        run_cmd(cmd)
                    else:
                        run_varscan_cmds.append(cmd)

                cmd = 'python test_varscan.py --convertVarScanRes2 '+\
                                              '--snpIn %s '%resFile_a+\
                                              '--snpOut %s '%resFile_b
                if convert_varscan==1:
                    if num_p==1:
                        run_cmd(cmd)
                    else:
                        convert_varscan_cmds.append(cmd)

    pdb.set_trace()
    if sam2bam==1 and num_p>1:
        run_cmds(sam2bam_cmds, num_p)
    if run_varscan==1 and num_p>1:
        run_cmds(run_varscan_cmds, num_p)
    if convert_varscan==1 and num_p>1:
        run_cmds(convert_varscan_cmds, num_p)

    return


'''
description:
build a dic, with key as tr id, val as [stt0, interval tree of tr, tr_seq, dir]
interval tree represents the transcript exon structure

refSeq is the seq of rel chrom
'''
def build_trTree(trBed, refSeq):

    trTree = {}

    with open(trBed, 'r') as f:

        for line in f:
            if line[0]=='#': continue
            tokens = line.split()

            trid = tokens[3]
            stt0 = int(tokens[1])
            direction = tokens[5]
            n_block = int(tokens[9])
            block_lens = [int(i) for i in tokens[10].split(',') if i != '']
            block_stts = [int(i) for i in tokens[11].split(',') if i != '']

            tree = IntervalTree()
            trTree[trid] = [stt0, tree, '', direction]

            for i in range(len(block_stts)):
                t_stt = block_stts[i]
                t_stp = block_stts[i]+block_lens[i]
                tree.add(Interval(t_stt, t_stp))

                trTree[trid][2] += refSeq[stt0+t_stt:stt0+t_stp]

    return trTree


'''
usage:
python test_varscan.py --convertVarScanRes
                       --snpIn snpInput
                       --snpOut snpOut

description:
extract snp calls into our snp res format for evaluation 
'''
def convertVarScanRes(args):

    snpIn = args[args.index('--snpIn')+1]
    snpOut = args[args.index('--snpOut')+1]

    with open(snpIn, 'r') as fin, open(snpOut, 'w') as fout:

        line = fin.readline()

        for line in fin:
            tokens = line.split()
            gPos1 = int(tokens[1])
            rB = tokens[2]
            tB = tokens[3]
            fout.write('%d\t%s\t-->\t%s\n'%(gPos1-1, rB, tB))

    print('%s written'%snpOut)

    return

'''
usage:
python test_varscan.py --convertVarScanRes2
                       --snpIn snpInput
                       --snpOut snpOut

description:
extract snp calls into our snp res format for evaluation 
- check tB not in {A,C,T,G} ==> IUPAC codes, need to further unpack (e.g. tB is K ==> G, T, if rB is G then tB is T)
- merge same pos (choose the one with highest var freq)
'''
def convertVarScanRes2(args):

    iupac = {} #key: iupac code val: possible a/c/g/t
    iupac['R']=['A','G']
    iupac['Y']=['C','T']
    iupac['S']=['G','C']
    iupac['W']=['A','T']
    iupac['K']=['G','T']
    iupac['M']=['A','C']

    snpIn = args[args.index('--snpIn')+1]
    snpOut = args[args.index('--snpOut')+1]

    with open(snpIn, 'r') as fin, open(snpOut, 'w') as fout:

        snps = {} #key - pos val - list of (rB, tB, varFreq)

        line = fin.readline()

        for line in fin:
            tokens = line.split()
            gPos1 = int(tokens[1])
            rB = tokens[2]
            tB = tokens[3]
            varFreq = float(tokens[6][:-1])/100 #skip %

            if tB not in ['A', 'C', 'T', 'G']:
                #pdb.set_trace()
                if tB in iupac:
                    if rB == iupac[tB][0]:
                        tB = iupac[tB][1]
                    elif rB == iupac[tB][1]:
                        tB = iupac[tB][0]
                    else:
                        print('tB=%s, rB=%s not in iupac[tB], skipped'%(tB, rB)); pdb.set_trace()
                        continue
                else:
                    print('tB=%s not in iupac, skipped'%tB); pdb.set_trace()
                    continue
            
            #if gPos1==25342544:
            #    pdb.set_trace()
            if gPos1 in snps:
                snps[gPos1].append((rB, tB, varFreq))#; pdb.set_trace()
            else:
                snps[gPos1]=[(rB, tB, varFreq)]

        itms = sorted(snps.items(), key=lambda x:x[0])
        #pdb.set_trace()
        for gPos1, rB_tB_varFreq_list in itms:
            if len(rB_tB_varFreq_list)>1:
                rtv_list = sorted(rB_tB_varFreq_list, key=lambda x:-x[2])# sort by VarFreq, decreasing
                print('gPos1=%d, rB_tB_varFreq_list=%s; pick one w/ highest VarFreq'%(gPos1, str(rB_tB_varFreq_list)))#; pdb.set_trace()
                rB, tB, varFreq = rtv_list[0]
                #pdb.set_trace()
            else:
                rB, tB, varFreq = rB_tB_varFreq_list[0]
            fout.write('%d\t%s\t-->\t%s\n'%(gPos1-1, rB, tB))

        #pdb.set_trace()

    print('%s written'%snpOut)

    return


'''
usage:
python test_varscan.py  --tloc2gloc
                        --refFa refFa           #reference genome in fasta file
                        --chrom chrom
                        --trBed trBed
                        --snpIn snpInputFile    #(varscan output)
                        [--log]                 #if enabled, debug purpose
                        --snpOut snpOutputFile  #(our format, in genome coordinates)

description:
convert snps called by aligning reads onto transcriptome into genome based coordinates

input:
--refFa refFa           #reference genome in fasta file
--chrom chrom
--trBed trBed
--snpIn snpInputFile    #(varscan output)
[--log]                 #if enabled, debug purpose, path = snpOutputFile + '.log'

output:
--snpOut snpOutputFile  #(our format, in genome coordinates)
'''

reverse_complement = lambda x: ''.join([{'A':'T','C':'G','G':'C','T':'A','N':'N'}[B] for B in x][::-1])

def tloc2gloc(args):

    refFa = args[args.index('--refFa')+1]
    chrom = args[args.index('--chrom')+1] #e.g. chr15
    trBed = args[args.index('--trBed')+1]
    snpIn = args[args.index('--snpIn')+1]
    snpOut = args[args.index('--snpOut')+1]
    
    if '--log' in args:
        log = snpOut+'.log'
    else:
        log = ''

    refGenome = from_fasta(refFa) #dic {chrom:seq}
    #pdb.set_trace()

    trTree = build_trTree(trBed, refGenome[chrom]) #dic {tr_id:[stt0, interval tree of tr, tr_seq, dir]}
    #pdb.set_trace()

    if log=='':
        logF = None
    else:
        logF = open(log, 'w')


    snps = {} #key - pos0 : [rB, tB] to avoid duplicate same snps from different transcripts sharing same exon

    with open(snpIn, 'r') as fin, open(snpOut, 'w') as fout:
        line = fin.readline() #skip the header
        for line in fin:
            tokens = line.split()
            trid = tokens[0]
            if trid not in trTree:
                print('trid not in trTree')
                pdb.set_trace()
            else:
                t_snpPos1 = int(tokens[1])
                rB = tokens[2]
                tB = tokens[3]
                if trTree[trid][3]=='+':
                    #pdb.set_trace()
                    g_snpPos0 = trTree[trid][0]+t_snpPos1-1                    
                    if g_snpPos0 in snps:
                        if rB != snps[g_snpPos0][0] or tB != snps[g_snpPos0][1]: pdb.set_trace()
                    else:
                        snps[g_snpPos0]=[rB, tB]
                        fout.write('%d\t%s\t-->\t%s\n'%(g_snpPos0, rB, tB))

                    '''print('tokens: %s'%str(tokens))
                    print('g_stt: %d'%trTree[trid][0])
                    print('intervals: %s'%str(trTree[trid][1]))
                    ss = trTree[trid][2]
                    print('seq: %s...%s...%s'%(ss[0:4], ss[t_snpPos1-2:t_snpPos1+1],ss[-4:]))
                    print('dir: %s'%trTree[trid][3])

                    pdb.set_trace()'''
                else:

                    #print('tokens: %s'%str(tokens))
                    #print('g_stt0: %d'%trTree[trid][0])
                    #print('intervals: %s'%str(trTree[trid][1]))
                    #ss = trTree[trid][2]
                    #ss = reverse_complement(ss)
                    #print('seq (rev comp): %s...%s...%s'%(ss[0:4], ss[t_snpPos1-2:t_snpPos1+1],ss[-4:]))
                    #print('dir: %s'%trTree[trid][3])

                    #print('modi to 1st strand')
                    ss = trTree[trid][2]
                    x0 = len(ss)-t_snpPos1
                    try:
                        rB_modi = reverse_complement(tokens[2])
                        tB_modi = reverse_complement(tokens[3])
                        g_snpPos0 = trTree[trid][0]+x0

                        if g_snpPos0 in snps:
                            if rB_modi != snps[g_snpPos0][0] or tB_modi != snps[g_snpPos0][1]: pdb.set_trace()
                        else:
                            snps[g_snpPos0]=[rB_modi, tB_modi]
                            fout.write('%d\t%s\t-->\t%s\n'%(g_snpPos0, rB_modi, tB_modi))
                    except:
                        pass; #pdb.set_trace()
                    #print('x0=%d, rB(1st strand)=%s, tB(1st strand)=%s'%(x0, rB_modi, tB_modi))
                    #print('seq: %s...%s...%s'%(ss[0:4], ss[x0-1:x0+2],ss[-4:]))
                    

                    #pdb.set_trace()

    print('%s written'%snpOut)
    #pdb.set_trace()

    if logF is not None:
        logF.close()

    return

'''
def analyze():

    #cmp_res_file = '/home/shunfu1/software/VarScan/run_5_n100k_err000_dedupped_res.txt'
    #cmp_res_file = '/home/shunfu1/software/VarScan/run_5_n100k_err000_rsem_res.txt'
    cmp_res_file = '/home/shunfu1/software/VarScan/run_5_n100k_err000_rsem_transcriptome_res.txt'

    ref = '/data1/shunfu1/SNPCalling/SimSNPs_MultiRun_5/SNP_p.txt'
    pos_idx = 0

    cmp_res = set()
    with open(cmp_res_file, 'r') as f:
        f.readline()
        for line in f:
            if line[0]=='#': continue
            tokens = line.split()
            cmp_res.add(int(tokens[1])-1)

    pos_ref = set()
    with open(ref, 'r') as f:
        for line in f:
            if line[0]=='#': continue
            tokens = line.split()
            pos_ref.add(int(tokens[pos_idx]))

    #pdb.set_trace()

    res = cmp_res.intersection(pos_ref)
    print(len(res))

    pdb.set_trace()

    return
'''


if __name__ == "__main__":

    #run_varscan()

    #analyze()

    args = sys.argv

    if '--tloc2gloc' in args:
        tloc2gloc(args)
    elif '--convertVarScanRes' in args:
        convertVarScanRes(args)
    elif '--convertVarScanRes2' in args:
        convertVarScanRes2(args)
    elif '--run_varscan' in args:
        run_varscan(args)
    elif '--run_varscan_batch' in args:
        run_varscan_batch(args)
    else:
        pdb.set_trace()



