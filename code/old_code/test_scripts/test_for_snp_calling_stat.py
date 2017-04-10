'''
call case6 (ours) and case2 (GATK) and snp analysis multiple times in order to get snp calling statistics in terms of mis-detection and false positive
'''

from util import *
from Address import *
import pdb

def main():

    #pdb.set_trace()

    condition = 'snp20_reads100k_moreResIncluded'

    a = 1 #inclusive
    b = 30 #exclusive
    cnt = 0

    for n in xrange(a,b):

        case = '%s_%d'%(condition, n)

        #prepare data
        cmd = 'mkdir -p %s'%Default_Ref_Path
        run_cmd(cmd)

        cmd = 'cp %s/* %s'%(Default_Genome_Sources, Default_Ref_Path)
        run_cmd(cmd)

        #run case 6
        cmd = 'python batch_run_case1plus6.py'
        run_cmd(cmd)

        #run case 2
        cmd = 'python batch_run_case2.py'
        run_cmd(cmd)

        #test snp analysis, res kept in code/tmp/
        cmd = 'python test_snp_analysis.py --condition %s '+\
              ' --round %d '%(condition, n)+\
              ' --copyRes '+\
              ' --thre 1 --compare 1'

        run_cmd(cmd)

        cnt += 1

        if False: #cnt < 10:
            #organize data
            cmd = 'mv %s %s/%s/'%(Default_Ref_Path, Default_Ref_Path_Root, case)
            run_cmd(cmd)
        else:
            cmd = 'rm -r %s'%Default_Ref_Path
            run_cmd(cmd)

    return

if __name__ == "__main__":
    main()