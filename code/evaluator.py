import sys
from old_code.snp_analysis import SnpNodeDic

'''
Description:
- use snp_analysis.py
- compare snp res among different callers and true snps
- loadSnpInfo: read snp info from true snps or called snps, add count/countAlt info if related count/countAlt files supplied (to check uniqueness)

Usage:

python evaluator.py    --loadSnpInfo 
                       -L1 label1 -F1 file1 [-L2 label2 -F2 file2 ...]
                       [-C1 countFile] [-C2 countAltFile]
                       [--snpLog outFile]
                       [--snpSum sumFile]
                       [--snpSum2 sumFile2]

'''

if __name__ == "__main__":

    args = sys.argv

    if '--loadSnpInfo' in args:

        snDic = SnpNodeDic()

        snDic.loadSnpInfo(args)

        if '--snpLog' in args:
            snDic.writeSnpLog(args)

        if '--snpSum' in args and '--snpSum2' in args:
            snDic.writeSnpSummary2(args) #snDic.writeSnpSummary(args)
