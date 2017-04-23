#!/usr/bin/python
#coding:utf-8

'''
Description:
- compare snp res among different callers and true snps

- loadSnpInfo: read snp info from true snps or called snps, add count/countAlt info if related count/countAlt files supplied (to check uniqueness)

Usage:

python snp_analysis.py --loadSnpInfo -L1 label1 -F1 file1 [-L2 label2 -F2 file2 ...] [-C1 countFile] [-C2 countAltFile] [--snpLog outFile] [--snpSum sumFile] [--snpSum2 sumFile2]

python snp_analysis.py --loadSnpLog --InSnpLog inputSnpLog [--snpLog outFile] [--snpSum sumFile] [--snpSum2 sumFile2]

'''

import sys, pdb, subprocess, re
try:
    from itertools import izip
except ImportError:  #python3.x
    izip = zip

'''
example: 
gPos = 28950303
rB = 'A'
tBs = {'T':'C', 'O':'', 'G':'C'} here 'T' means true snp is 'C'; 'O' means our caller miss it; 'G' means 'GATK' calls it as 'C'
'''
class SnpNode:

    def __init__(self, gPos, rB, labs):

        self.gPos = gPos #genome loc
        self.rB = rB #ref base
        self.tBs = {} #key: lab (e.g. true snp 'T'), val: target base
        for lab in labs:
            self.tBs[lab]=''
        self.countStr = ''
        self.countAltStr = ''
        self.ReadsUniqMapped = -1 #-1: undefined; 0: reads that aligned over this snp is non-uniq mapped; 1: uniq mapped
        self.notCapturedByCount = 0 #1: this snp info is not captured by count/countAlt

    def add(self, currLab, tB):

        if currLab in self.tBs:
            self.tBs[currLab] = tB
        else:
            print('SnpNode add: unknown currLab %s'%currLab)
            pdb.set_trace()

    '''
    count example:
    2440687 82798325    G   1.874   5   0   49  0   G,I,55.05,56.93 G,I,55.05,56.93 G,I,45.94,47.82 ...
    '''
    def addCountStr(self, countFile):
        if countFile=='':
            return
        try:
            res = subprocess.check_output('grep \'%d\' %s'%(self.gPos, countFile), shell=True)
            res = res.split('\n')
            for r in res:
                if r!='' and self.gPos==int(r.split()[1]):
                    self.countStr = r
        except:
            pass

        return

    def addCountStr2(self, countStr):
        self.countStr = countStr.strip()
        return


    '''
    count example:
    2440687 82798325,1.874  [-A:85777969,0.0,83023818,0.0,83182903,45.9440162226,84898863,9.10984643735,]...
    '''
    def addCountAltStr(self, countAltFile):
        if countAltFile=='':
            return
        try:
            res = subprocess.check_output('grep \'%d\' %s'%(self.gPos, countAltFile), shell=True)
            res = res.split('\n')
            for r in res:
                if r!= '' and self.gPos==int((r.split()[1]).split(',')[0]):
                    self.countAltStr = r
        except:
            pass

        return

    def addCountAltStr2(self, countAltStr):
        self.countAltStr = countAltStr.strip()
        return


    def updateUniqMapping(self):

        if self.countAltStr == '':
            return

        #pdb.set_trace()

        if len(self.countAltStr.split())>=3:
            self.ReadsUniqMapped = 0
        else:
            self.ReadsUniqMapped = 1

        return



    #check (grep) gPos from list of files
    #
    #e.g. check at gPos, if SNP is called in different files
    #i.e. print: grep 'gPos' each_file
    #
    #return a dic of grep res (key: file, val: grep str)
    '''
    def check(self, files):

        resDic = {}

        for file in files:

            try:
                res = subprocess.check_output('grep \'%d\' %s'%(self.gPos, file), shell=True)
                resDic[file] = res
            except:
                resDic[file] = ''
                continue

        return resDic
    '''

    def __str__(self):

        st = ''
        st += 'gPos=%d '%self.gPos
        st += 'rB=%s '%self.rB
        st += 'tBs=['
        itms = list(self.tBs.items())
        itms.sort(key=lambda x: x[0])
        for lab, tB in itms:
            st += '%s:%s, '%(lab, tB)
        st += ']'        
        return st

    def fulStr(self): #return a full description of SnpNode

        st = self.__str__()
        st += ' ReadsUniqMapped=%d'%self.ReadsUniqMapped
        if self.notCapturedByCount==1:
            st += ' [not captured by count]'
        st += '\n--------------------\ncount:\n%s\n--------------------\n'%(self.countStr)
        st += 'countAlt:\n%s\n--------------------\n'%(self.countAltStr)

        return st

class SnpNodeDic():

    def __init__(self):
        self.snDic = {} #key - genome location val - SnpNode

        self.labs = []
        self.files = []

        self.countFile = ''# \n deleted from count/ countAlt str
        self.countAltFile = ''# \n deleted from count/ countAlt str

        self.snpLog = ''
        self.snpSum = ''
        self.snpSum2 = ''

        return

    #
    # load snps from different files (e.g. true snps, called snps) into a SnpNodeDic object
    # add count/countAlt info if related count/countAlt files supplied (to check uniqueness)
    # output file onto outFile
    #
    def loadSnpInfo(self, args):

        #load snp info
        i = 1
        while (i>0):
            lab = '-L%d'%i
            if lab in args:
                lab = args[args.index(lab)+1]
                self.labs.append(lab)
                i = i+1
            else:
                break

        for i in range(len(self.labs)):
            lab = self.labs[i]
            file = args[args.index('-F%d'%(i+1))+1]
            self.files.append(file)
            self.loadOneSnpFile(lab, file)

        #load count/ countAltInfo

        if '-C1' in args:
            self.countFile = args[args.index('-C1')+1]

        if '-C2' in args:
            self.countAltFile = args[args.index('-C2')+1]

        #pdb.set_trace()
        if '-C1' in args and '-C2' in args:
            nLines=sum([1 for l in open(self.countFile,'r')]); T=nLines/100; p=0; q=0;
            with open(self.countFile) as c_F, open(self.countAltFile) as cA_F:
                
                for xline, yline in izip(c_F, cA_F):
                    p += 1
                    if p>=T: p=0; q+=1; sys.stdout.write('\r'); sys.stdout.write('%d %% processed (loadSnpInfo countFile)'%q); sys.stdout.flush()

                    x = xline.split()
                    y = yline.split()
                    gPos = int(x[1])
                    gPos_cA = int((y[1]).split(',')[0])
                    if gPos != gPos_cA:
                        print('unexpected count/count_altInfo (unequal gPos)')
                        pdb.set_trace()
                    elif gPos in self.snDic:
                        #pdb.set_trace()
                        sn = self.snDic[gPos]
                        sn.addCountStr2(xline.strip())# \n deleted from count/ countAlt str
                        sn.addCountAltStr2(yline.strip())
                        sn.updateUniqMapping()
                    else:
                        continue

        for gPos, sn in self.snDic.items(): #it's possible that sn is not captured by count/countAlt
            if sn.ReadsUniqMapped==-1:
                #print('gPos not captured by count: %s (ReadsUniqMapped set as 1)'%str(sn))
                sn.ReadsUniqMapped = 1
                sn.notCapturedByCount = 1


        # old method
        '''
        for gPos, sn in self.snDic.items():
            print('add count/countAlt for %s'%str(sn))
            #if gPos == 82798325: pdb.set_trace()
            sn.addCountStr(self.countFile)
            sn.addCountAltStr(self.countAltFile)
            sn.updateUniqMapping()
        '''

        return [self.snDic, self.labs, self.files]

    #
    # load snps from a snpLog file (which contains snps, count/countAlt info) into a SnpNodeDic object
    # 
    # mainly used to improve snp_snalysis.py (this file), on local pc instead of on server
    #
    def loadSnpLog(self, args):

        snpLog = args[args.index('--InSnpLog')+1]

        #load snp info
        with open(snpLog, 'r') as sL:
        #with open(snpLog, 'r') as sL, open('tmp/loadSnpLog/test_snpLog_dmp.txt', 'w') as sLdmp:
            for line in sL:
                if len(line)>=4 and line[0:4]=='gPos':
                    line_snp = line #e.g. line_snp is: gPos=42103371 rB=G tBs=[p:, m:, G:C, O:, ] ReadsUniqMapped=1
                    sL.readline(); sL.readline(); line_count = sL.readline().strip() #no \n
                    sL.readline(); sL.readline(); line_countAlt = sL.readline().strip()
                    '''sLdmp.write(line)
                    sLdmp.write(line_count)
                    sLdmp.write(line_countAlt)
                    sLdmp.write('\n')'''

                    b=re.search('gPos=([\S]+)', line_snp)
                    if b is not None:
                        gPos = int(b.group(1))
                    else:
                        print('search for gPos fail in line_snp: %s'%line_snp)
                        pdb.set_trace()

                    b=re.search('rB=([\S]+)', line_snp)
                    if b is not None:
                        rB = (b.group(1))
                    else:
                        print('search for rB fail in line_snp: %s'%line_snp)
                        pdb.set_trace()

                    b=re.search('tBs=\[(.+)\]', line_snp)
                    if b is not None:
                        tBs = (b.group(1))
                        tBs = [[itm.split(':')[0].strip(), itm.split(':')[1].strip()] for itm in tBs.split(',') if itm.strip()!='']
                    else:
                        print('search for tBs fail in line_snp: %s'%line_snp)
                        pdb.set_trace()

                    b=re.search('ReadsUniqMapped=([\S]+)', line_snp)
                    if b is not None:
                        uniq = int(b.group(1))
                    else:
                        print('search for ReadsUniqMapped fail in line_snp: %s'%line_snp)
                        pdb.set_trace()

                    #sLdmp.write(line.strip()); sLdmp.write('\n'+str(gPos)); sLdmp.write('\n'+str(rB)); sLdmp.write('\n'+str(tBs)); sLdmp.write('\n'+str(uniq)+'\n\n')

                    labs = [tB[0] for tB in tBs]
                    self.snDic[gPos] = SnpNode(gPos, rB, labs)
                    for currLab, tB in tBs:
                        self.snDic[gPos].add(currLab, tB)
                    self.snDic[gPos].addCountStr2(line_count)# \n deleted from count/ countAlt str
                    self.snDic[gPos].addCountAltStr2(line_countAlt)
                    self.snDic[gPos].updateUniqMapping()

                else:
                    continue

                #pdb.set_trace()

        for gPos, sn in self.snDic.items(): #it's possible that sn is not captured by count/countAlt
            if sn.ReadsUniqMapped==-1:
                print('gPos not captured by count: %s (ReadsUniqMapped set as 1)'%str(sn))
                sn.ReadsUniqMapped = 1
                sn.notCapturedByCount = 1

        return [self.snDic, self.labs, self.files]


    # file format: e.g. "22925895        C       -->     A"
    def loadOneSnpFile(self, currLab, file):

        with open(file, 'r') as f:

            for line in f:

                tokens = line.split()

                gPos = int(tokens[0])
                rB = tokens[1].upper() # ref Base
                tB = tokens[3].upper() # target SNP Base

                if gPos in self.snDic:

                    self.snDic[gPos].add(currLab, tB)

                else:

                    self.snDic[gPos] = SnpNode(gPos, rB, self.labs)
                    self.snDic[gPos].add(currLab, tB)
        return

    # check different types of snps (e.g. md by caller1 but not caller 2; fp by caller1 and caller2 etc)
    '''
    def checkSnpNodeDic(self):

        print('\ncheckSnpNodeDic:')

        itms = self.snDic.items()
        itms.sort(key=lambda x:x)

        i=0
        for gPos, sn in itms:
            #if sn.tBs['p']!='' and sn.tBs['m']=='' and sn.tBs['O']=='' and sn.tBs['G']!='':# p md by O but not G
            #if sn.tBs['p']!='' and sn.tBs['m']=='' and sn.tBs['O']!='' and sn.tBs['G']=='':# p md by G but not O
            #if sn.tBs['p']!='' and sn.tBs['m']=='' and sn.tBs['O']=='' and sn.tBs['G']=='':# p md by G and O
            #if sn.tBs['p']!='' and sn.tBs['m']=='' and sn.tBs['O']!='' and sn.tBs['G']!='':# p detected by G and O
            #if sn.tBs['p']=='' and sn.tBs['m']=='' and sn.tBs['O']!='' and sn.tBs['G']=='':# p fp by O but not G
            if sn.tBs['p']=='' and sn.tBs['m']=='' and sn.tBs['O']=='' and sn.tBs['G']!='':# p fp by G but not O
            #if sn.tBs['p']=='' and sn.tBs['m']=='' and sn.tBs['O']!='' and sn.tBs['G']!='':# p fp by G and O
                i+=1
                print('%d\t%s'%(i,sn))
                print('\n--------------------')
                check_res = sn.check(self.files)
                print('\n'.join([s[0]+'\n'+s[1] for s in check_res.items()]))
                print('--------------------\n')

        return
    '''


    def writeSnpLog(self, args):

        if '--snpLog' in args:
            outFile = args[args.index('--snpLog')+1]
            with open(outFile, 'w') as oF:
                itms = list(self.snDic.items())
                #pdb.set_trace()
                itms.sort(key=lambda x : x[0])
                for gPos, sn in itms:
                    oF.write(sn.fulStr()+'\n')
        return

    #get mis-detection/false-positive summary
    #
    # O - our snp caller
    # G - GATK or other caller
    #
    #esp, 
    #
    #in snpSum file:
    #snp w/ only uniq mapped reads: O-md O-fp G-md G-fp
    #snp w/ non-uniq mapped reads:  O-md O-fp G-md G-fp
    #overall:                       O-md O-fp G-md G-fp
    #
    #in snpSum2 file:
    #snp w/ only uniq mapped reads: O-md\G-md: list of gPos
    #                               O-md&G-md: list of gPos
    #                               G-md\O-md: list of gPos
    #                               O-fp\G-fp: list of gPos
    #                               O-fp&G-fp: list of gPos
    #                               G-fp\O-fp: list of gPos
    #
    #snp w/ non-uniq mapped reads:  O-md\G-md: list of gPos
    #                               O-md&G-md: list of gPos
    #                               G-md\O-md: list of gPos
    #                               O-fp\G-fp: list of gPos
    #                               O-fp&G-fp: list of gPos
    #                               G-fp\O-fp: list of gPos
    #
    # snpSum and snpSum2 files need both to be specified by [--snpSum sumFile] [--snpSum2 sumFile2]
    '''
    def writeSnpSummary(self, args):

        snpSum = ''
        if '--snpSum' in args:
            snpSum = args[args.index('--snpSum')+1]

        snpSum2 = ''
        if '--snpSum2' in args:
            snpSum2 = args[args.index('--snpSum2')+1]

        if snpSum=='' or snpSum2=='':
            print('need to specify paths for snpSum and snpSum2')
            pdb.set_trace()
            return

        #pdb.set_trace()

        k_dscrpt_pairs = []

        k_dscrpt_pairs.append([(1, 1, 0, 1), 'U[O-md\G-md]']) #'U[O-md\G-md]'： snp has only uniq mapped reads, snp is true snp, mis-detected by our caller, called by GATK
        k_dscrpt_pairs.append([(1, 1, 0, 0), 'U[O-md&G-md]'])
        k_dscrpt_pairs.append([(1, 1, 1, 0), 'U[G-md\O-md]'])
        k_dscrpt_pairs.append([(1, 0, 1, 0), 'U[O-fp\G-fp]'])
        k_dscrpt_pairs.append([(1, 0, 1, 1), 'U[O-fp&G-fp]'])
        k_dscrpt_pairs.append([(1, 0, 0, 1), 'U[G-fp\O-fp]'])

        k_dscrpt_pairs.append([(0, 1, 0, 1), 'N[O-md\G-md]']) #'N[O-md\G-md]'： snp has non uniq mapped reads, snp is true snp, mis-detected by our caller, called by GATK
        k_dscrpt_pairs.append([(0, 1, 0, 0), 'N[O-md&G-md]'])
        k_dscrpt_pairs.append([(0, 1, 1, 0), 'N[G-md\O-md]'])
        k_dscrpt_pairs.append([(0, 0, 1, 0), 'N[O-fp\G-fp]'])
        k_dscrpt_pairs.append([(0, 0, 1, 1), 'N[O-fp&G-fp]'])
        k_dscrpt_pairs.append([(0, 0, 0, 1), 'N[G-fp\O-fp]'])

        k2dscrpt = {}
        dscrpt2k = {}

        for k, dscrpt in k_dscrpt_pairs:
            k2dscrpt[k] = dscrpt
            dscrpt2k[dscrpt] = k

        #pdb.set_trace()  

        res = {} #key: (uniq_reads_mapping, true_snp, caller by our caller, called by GATK) val=list of gPos of related snps
        for k, dscrpt in k2dscrpt.items():
            res[k] = [] # list of str pos (if pos is not cap by count, denoted as pos*)

        #pdb.set_trace()  

        for gPos, sn in self.snDic.items():

            key = [0,0,0,0]

            if sn.ReadsUniqMapped==-1:
                print('%s undefined read uniq mapping or not'%str(sn))
                continue #pdb.set_trace()

            elif sn.ReadsUniqMapped==1:
                key[0]=1

            else:
                pass

            if sn.tBs['p']!='' or sn.tBs['m']!='':
                key[1]=1

            if sn.tBs['O']!='':
                key[2]=1

            if sn.tBs['G']!='':
                key[3]=1

            key = tuple(key)
            if key[1]==1 and key[2]==1 and key[3]==1: #correct detected
                pass
            elif key not in res:
                print('undefined key=%s, sn=%s'%(str(key), str(sn)))
            else:
                if sn.notCapturedByCount==0:
                    res[key].append(str(gPos))
                else:
                    res[key].append(str(gPos)+'*')

        #pdb.set_trace()  

        with open(snpSum, 'w') as sF:

            st = 'O-md\tO-fp\tG-md\tG-fp\t(snps of uniq mapped reads)\n'
            omd = len(res[dscrpt2k['U[O-md\G-md]']]+res[dscrpt2k['U[O-md&G-md]']])
            ofp = len(res[dscrpt2k['U[O-fp\G-fp]']]+res[dscrpt2k['U[O-fp&G-fp]']])
            gmd = len(res[dscrpt2k['U[G-md\O-md]']]+res[dscrpt2k['U[O-md&G-md]']])
            gfp = len(res[dscrpt2k['U[G-fp\O-fp]']]+res[dscrpt2k['U[O-fp&G-fp]']])
            st += '%d\t%d\t%d\t%d\n'%(omd, ofp, gmd, gfp)
            sF.write(st)

            st = 'O-md\tO-fp\tG-md\tG-fp\t(snps of NON uniq mapped reads)\n'
            omd2 = len(res[dscrpt2k['N[O-md\G-md]']]+res[dscrpt2k['N[O-md&G-md]']])
            ofp2 = len(res[dscrpt2k['N[O-fp\G-fp]']]+res[dscrpt2k['N[O-fp&G-fp]']])
            gmd2 = len(res[dscrpt2k['N[G-md\O-md]']]+res[dscrpt2k['N[O-md&G-md]']])
            gfp2 = len(res[dscrpt2k['N[G-fp\O-fp]']]+res[dscrpt2k['N[O-fp&G-fp]']])
            st += '%d\t%d\t%d\t%d\n'%(omd2, ofp2, gmd2, gfp2)
            sF.write(st)

            st = 'O-md\tO-fp\tG-md\tG-fp\t(tot)\n'
            st += '%d\t%d\t%d\t%d\n'%(omd+omd2, ofp+ofp2, gmd+gmd2, gfp+gfp2)
            sF.write(st)

        #pdb.set_trace()

        with open(snpSum2, 'w') as sF:

            st = ''
            for k, dscrpt in k_dscrpt_pairs:
                st += dscrpt + ':' + ','.join([p for p in res[k]])
                st += '\n'
            st += '\n'

            sF.write(st)

        #pdb.set_trace()  

        return
    '''
    
    #
    #content_idx: 0,1,2 corresponding to snp w/ only uniq mapped reads, w/o uniq mapped reads and overall
    #caller_idx: 0,1 corresponding to O, G
    #returns related [cd, md, fp, toterr]
    #
    #in case of err (e.g. file not exist), return []
    def extract_snpSum(self, file, content_idx, caller_idx):
        res = []
        try:
            with open(file, 'r') as f:
                line = ''
                for i in range(content_idx):
                    line = f.readline()
                    line = f.readline()
                line = f.readline()
                line = f.readline()
                if len(line.strip().split())>=8:
                    vals = [int(v) for v in line.strip().split()]
                    res = vals[0+4*caller_idx:4+4*caller_idx]
                    return res
                else:
                    print('unexpected line: %s in extract_snpSum'%line)
                    pdb.set_trace()
                    return []
        except:
            return []


    #get mis-detection/false-positive summary
    #modified based on writeSnpSummary(self, args) to add stat for correct detection, tot error
    #
    # O - our snp caller
    # G - GATK or other caller
    #
    #esp,
    #
    #in snpSum file:
    #snp w/ only uniq mapped reads: O-cd O-md O-fp O-toterr G-cd G-md G-fp G-toterr
    #snp w/ non-uniq mapped reads:  O-cd O-md O-fp O-toterr G-cd G-md G-fp G-toterr
    #overall:                       O-cd O-md O-fp O-toterr G-cd G-md G-fp G-toterr
    #
    #in snpSum2 file:
    #snp w/ only uniq mapped reads: O-cd\G-cd: list of gPos ==> discarded, same as G-md\O-md
    #                               O-cd&G-cd: list of gPos
    #                               G-cd\O-cd: list of gPos ==> discarded
    #                               O-md\G-md: list of gPos
    #                               O-md&G-md: list of gPos
    #                               G-md\O-md: list of gPos
    #                               O-fp\G-fp: list of gPos
    #                               O-fp&G-fp: list of gPos
    #                               G-fp\O-fp: list of gPos
    #
    #snp w/ non-uniq mapped reads:  O-cd\G-cd: list of gPos ==> discarded
    #                               O-cd&G-cd: list of gPos
    #                               G-cd\O-cd: list of gPos ==> discarded
    #                               O-md\G-md: list of gPos
    #                               O-md&G-md: list of gPos
    #                               G-md\O-md: list of gPos
    #                               O-fp\G-fp: list of gPos
    #                               O-fp&G-fp: list of gPos
    #                               G-fp\O-fp: list of gPos
    #
    # snpSum and snpSum2 files need both to be specified by [--snpSum sumFile] [--snpSum2 sumFile2]

    def writeSnpSummary2(self, args):

        snpSum = ''
        if '--snpSum' in args:
            snpSum = args[args.index('--snpSum')+1]

        snpSum2 = ''
        if '--snpSum2' in args:
            snpSum2 = args[args.index('--snpSum2')+1]

        if snpSum=='' or snpSum2=='':
            print('need to specify paths for snpSum and snpSum2')
            pdb.set_trace()
            return

        #pdb.set_trace()

        k_dscrpt_pairs = []

        #k_dscrpt_pairs.append([(1, 1, 1, 0), 'U[O-cd\G-cd]'])
        k_dscrpt_pairs.append([(1, 1, 1, 1), 'U[O-cd&G-cd]'])
        #k_dscrpt_pairs.append([(1, 1, 0, 1), 'U[G-cd\O-cd]'])
        k_dscrpt_pairs.append([(1, 1, 0, 1), 'U[O-md\G-md]']) #'U[O-md\G-md]'： snp has only uniq mapped reads, snp is true snp, mis-detected by our caller, called by GATK
        k_dscrpt_pairs.append([(1, 1, 0, 0), 'U[O-md&G-md]'])
        k_dscrpt_pairs.append([(1, 1, 1, 0), 'U[G-md\O-md]'])
        k_dscrpt_pairs.append([(1, 0, 1, 0), 'U[O-fp\G-fp]'])
        k_dscrpt_pairs.append([(1, 0, 1, 1), 'U[O-fp&G-fp]'])
        k_dscrpt_pairs.append([(1, 0, 0, 1), 'U[G-fp\O-fp]'])

        #k_dscrpt_pairs.append([(0, 1, 1, 0), 'N[O-cd\G-cd]'])
        k_dscrpt_pairs.append([(0, 1, 1, 1), 'N[O-cd&G-cd]'])
        #k_dscrpt_pairs.append([(0, 1, 0, 1), 'N[G-cd\O-cd]'])
        k_dscrpt_pairs.append([(0, 1, 0, 1), 'N[O-md\G-md]']) #'N[O-md\G-md]'： snp has non uniq mapped reads, snp is true snp, mis-detected by our caller, called by GATK
        k_dscrpt_pairs.append([(0, 1, 0, 0), 'N[O-md&G-md]'])
        k_dscrpt_pairs.append([(0, 1, 1, 0), 'N[G-md\O-md]'])
        k_dscrpt_pairs.append([(0, 0, 1, 0), 'N[O-fp\G-fp]'])
        k_dscrpt_pairs.append([(0, 0, 1, 1), 'N[O-fp&G-fp]'])
        k_dscrpt_pairs.append([(0, 0, 0, 1), 'N[G-fp\O-fp]'])

        k2dscrpt = {}
        dscrpt2k = {}

        for k, dscrpt in k_dscrpt_pairs:
            k2dscrpt[k] = dscrpt
            dscrpt2k[dscrpt] = k

        #pdb.set_trace()

        res = {} #key: (uniq_reads_mapping, true_snp, caller by our caller, called by GATK) val=list of gPos of related snps
        for k, dscrpt in k2dscrpt.items():
            res[k] = []

        #pdb.set_trace()

        for gPos, sn in self.snDic.items():

            key = [0,0,0,0]

            if sn.ReadsUniqMapped==-1:
                print('%s undefined read uniq mapping or not'%str(sn))
                continue #pdb.set_trace()

            elif sn.ReadsUniqMapped==1:
                key[0]=1

            else:
                pass

            if sn.tBs['p']!='' or sn.tBs['m']!='':
                key[1]=1

            if sn.tBs['O']!='':
                key[2]=1

            if sn.tBs['G']!='':
                key[3]=1

            key = tuple(key)
            #if key[1]==1 and key[2]==1 and key[3]==1: #correct detected
            #    pass
            #elif key not in res:
            if key not in res:
                print('undefined key=%s, sn=%s'%(str(key), str(sn)))
            else:
                if sn.notCapturedByCount==0:
                    res[key].append(str(gPos))
                else:
                    res[key].append(str(gPos)+'*')

        #pdb.set_trace()

        with open(snpSum, 'w') as sF:

            st = 'O-cd\tO-md\tO-fp\tO-toterr\tG-cd\tG-md\tG-fp\tG-toterr\t(snps of uniq mapped reads)\n'
            ocd = len(res[dscrpt2k['U[O-cd&G-cd]']]+res[dscrpt2k['U[G-md\O-md]']])
            omd = len(res[dscrpt2k['U[O-md\G-md]']]+res[dscrpt2k['U[O-md&G-md]']])
            ofp = len(res[dscrpt2k['U[O-fp\G-fp]']]+res[dscrpt2k['U[O-fp&G-fp]']])
            gcd = len(res[dscrpt2k['U[O-cd&G-cd]']]+res[dscrpt2k['U[O-md\G-md]']])
            gmd = len(res[dscrpt2k['U[G-md\O-md]']]+res[dscrpt2k['U[O-md&G-md]']])
            gfp = len(res[dscrpt2k['U[G-fp\O-fp]']]+res[dscrpt2k['U[O-fp&G-fp]']])
            st += '%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n'%(ocd, omd, ofp, omd+ofp, gcd, gmd, gfp, gmd+gfp)
            sF.write(st)

            st = 'O-cd\tO-md\tO-fp\tO-toterr\tG-cd\tG-md\tG-fp\tG-toterr\t(snps of NON uniq mapped reads)\n'
            ocd2 = len(res[dscrpt2k['N[O-cd&G-cd]']]+res[dscrpt2k['N[G-md\O-md]']])
            omd2 = len(res[dscrpt2k['N[O-md\G-md]']]+res[dscrpt2k['N[O-md&G-md]']])
            ofp2 = len(res[dscrpt2k['N[O-fp\G-fp]']]+res[dscrpt2k['N[O-fp&G-fp]']])
            gcd2 = len(res[dscrpt2k['N[O-cd&G-cd]']]+res[dscrpt2k['N[O-md\G-md]']])
            gmd2 = len(res[dscrpt2k['N[G-md\O-md]']]+res[dscrpt2k['N[O-md&G-md]']])
            gfp2 = len(res[dscrpt2k['N[G-fp\O-fp]']]+res[dscrpt2k['N[O-fp&G-fp]']])
            st += '%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n'%(ocd2, omd2, ofp2, omd2+ofp2, gcd2, gmd2, gfp2, gmd2+gfp2)
            sF.write(st)

            st = 'O-cd\tO-md\tO-fp\tO-toterr\tG-cd\tG-md\tG-fp\tG-toterr\t(tot)\n'
            st += '%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n'%(ocd+ocd2, omd+omd2, ofp+ofp2, omd+omd2+ofp+ofp2,
                                                      gcd+gcd2, gmd+gmd2, gfp+gfp2, gmd+gmd2+gfp+gfp2)
            sF.write(st)

        #pdb.set_trace()

        with open(snpSum2, 'w') as sF:

            st = ''
            for k, dscrpt in k_dscrpt_pairs:
                st += dscrpt + ':' + ','.join([p for p in res[k]])
                st += '\n'
            st += '\n'

            sF.write(st)

        #pdb.set_trace()

        return

if __name__ == "__main__":

    args = sys.argv

    #test loadSnpLog
    #python snp_analysis.py --loadSnpLog --InSnpLog inputSnpLog
    ''' test_data1
    args = ['--loadSnpLog',
            '--InSnpLog', 'tmp/loadSnpLog/test_snpLog.txt',
            '--snpLog', 'tmp/loadSnpLog/test_snpLog_dumpAfterLoad.txt',
            '--snpSum', 'tmp/loadSnpLog/test_snpLog_sum.txt',
            '--snpSum2', 'tmp/loadSnpLog/test_snpLog_sum2.txt']
    '''
    ''' test_data2
    args = ['--loadSnpLog',
            '--InSnpLog', 'tmp/loadSnpLog/snpLog_snp20_reads100k_4.txt',
            '--snpLog',   'tmp/loadSnpLog/snpLog_snp20_reads100k_4_readAndDump.txt',
            '--snpSum', 'tmp/loadSnpLog/snpSum_snp20_reads100k_4.txt',
            '--snpSum2', 'tmp/loadSnpLog/snpSum2_snp20_reads100k_4.txt']
    '''
    ''''
    args = ['--loadSnpLog',
            '--InSnpLog', 'tmp/loadSnpLog/test_data3_snp1k_reads10m_varyT_rsemCount2/snpLog_T1_c1.txt',
            '--snpLog',   'tmp/loadSnpLog/snpLog_T1_c1_readAndDump.txt',
            '--snpSum', 'tmp/loadSnpLog/snpLog_T1_c1_snpSum.txt',
            '--snpSum2', 'tmp/loadSnpLog/snpLog_T1_c1_snpSum2.txt']
    '''

    '''args = '--loadSnpInfo '+\
           '-L1 p -F1 /data1/shunfu1/SNPCalling/data/ref_snp/chr15_snp_p.txt '+\
           '-L2 m -F2 /data1/shunfu1/SNPCalling/data/ref_snp/chr15_snp_m.txt '+\
           '-L3 O -F3 /data1/shunfu1/SNPCalling/snp20_reads100k_10/data_GATK_flexiblePipeline2_RsemSimReadsEncMod/GATK_out/raw_variants.vcf.txt '+\
           '-L4 G -F4 /data1/shunfu1/SNPCalling/snp20_reads100k_10/data_GATK_flexiblePipeline2_RsemSimReadsEncMod/GATK_out/raw_variants.vcf.txt '+\
           '--snpLog /data1/shunfu1/SNPCalling/snp20_reads100k_10/data_GATK_flexiblePipeline2_RsemSimReadsEncMod/GATK_out/snp_res.log '+\
           '--snpSum /data1/shunfu1/SNPCalling/snp20_reads100k_10/data_GATK_flexiblePipeline2_RsemSimReadsEncMod/GATK_out/snp_res.sum '+\
           '--snpSum2 /data1/shunfu1/SNPCalling/snp20_reads100k_10/data_GATK_flexiblePipeline2_RsemSimReadsEncMod/GATK_out/snp_res.sum2'
    args = args.split()
    print('modified args: %s'%str(args)); pdb.set_trace()
    '''

    if '--loadSnpInfo' in args:

        snDic = SnpNodeDic()

        snDic.loadSnpInfo(args)

        if '--snpLog' in args:
            snDic.writeSnpLog(args)

        if '--snpSum' in args and '--snpSum2' in args:
            snDic.writeSnpSummary2(args) #snDic.writeSnpSummary(args)

    elif '--loadSnpLog' in args: #test for writeSnpSummary2

        #pdb.set_trace()

        snDic = SnpNodeDic()

        snDic.loadSnpLog(args)

        if '--snpLog' in args:
            snDic.writeSnpLog(args)

        if '--snpSum' in args and '--snpSum2' in args:
            snDic.writeSnpSummary2(args)

        #pdb.set_trace()



