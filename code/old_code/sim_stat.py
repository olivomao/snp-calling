# -*- coding: utf-8 -*-
"""
Created on Thu Aug 13 18:50:05 2015

@author: olivo
"""
import pdb
import numpy as np
#from scipy import stats
from Address import *
import subprocess

if isServerSide==False:
    import pylab as pl

class sim_stat:
    """
    object to
    * store the data of interest during simulation
    in order to do statistical analysis
    
    * be careful about memory
    * FILE I/O    
    """

    def __init__(self, dmp_fn='/sim_stat_dmp.txt'):
        print('sim_stat object created')
        # file i/o
        self.dmp_address = Default_Ref_Path + dmp_fn
        
        print('File %s'%self.dmp_address)
        c = 0 #int(input('Remove it? 1-yes 0-no:'))
        if c==1:
            subprocess.call( 'rm ' +  self.dmp_address,  shell=True ) #clear old file
        
        # track statistics during Synthesis->ExpressionLevel2Coverage
        self.acc_sign=[]        
        self.acc_cover=[]
        
        self.qt=[] # quantile
        self.acc_cover_qt=[] # acc_cover wrt quantile
        
        # track if SNP is covered by reads
        self.num_snp_m_covered_by_reads=0
        self.num_snp_p_covered_by_reads=0
        
    def set_dmp_addr(self, addr):
        self.dmp_address = addr
        
    
    """
    append info in obj, based on description
    """
    def dmp(self, obj, description): #obj: list of list
        dmp_file = open(self.dmp_address, 'a+') #append
        
        dmp_file.write(description+'\n')
        
        for i in range(len(obj)):
            item = obj[i]
            
            for j in range(len(item)):
                dmp_file.write(str(item[j])+'\t')
            
            dmp_file.write('\n')
        
        dmp_file.write('\n')
        
        dmp_file.close()
        
        
    """
    append info in obj, based on description
    """
    def dmp_vec(self, obj, description): #obj: list/array of str/int/float
        dmp_file = open(self.dmp_address, 'a+') #append
        
        dmp_file.write(description+'\n')
        
        for i in range(len(obj)):
            item = obj[i]
            
            dmp_file.write(str(item))
            
            dmp_file.write('\n')
        
        dmp_file.write('\n')
        
        dmp_file.close()
        
    """
    append info in obj, based on description
    """
    def dmp_scaler(self, item, description): #obj: list/array of str/int/float
        dmp_file = open(self.dmp_address, 'a+') #append
        
        dmp_file.write(description+'\n')
        
        dmp_file.write(str(item))
        
        dmp_file.write('\n')
        
        dmp_file.close()
        
    def set_qt(self, qt):
        self.qt = qt
        
    def set_acc_cov_qt(self):        
        self.acc_cover_qt = np.percentile(self.acc_cover, self.qt)
        
    def get_acc_cov_hd(self, COV_fn='/coverage.txt'):
        if len(self.acc_cover)>0:
            print('Stat:get_acc_cov_hd:acc_cover already init')
            return
            
        cov_address = Default_Ref_Path + COV_fn
        with open(cov_address) as cov_file:
            for line in cov_file:
                x = line.split()
                self.acc_cover.append(float(x[4]))
    
    """
    def dump_stats_Synthesis(self):
        
        #pdb.set_trace()
        # check acc_sign
        if len(self.acc_sign) > 0:
            acc_sign_arr = stats.itemfreq(np.asarray(self.acc_sign))
            L = len(acc_sign_arr[:,1])
            S = sum(acc_sign_arr[:,1])
            acc_sign_list = [[int(acc_sign_arr[i][0]), int(acc_sign_arr[i][1]), float(acc_sign_arr[i][1])/S] for i in range(L)]
            self.dmp(acc_sign_list, 'acc_sign vals:freq:percentage')
        else:
            self.dmp([], 'acc_sign vals:freq:percentage')
        
        #pdb.set_trace()
        # check acc_cover
        #if len(self.acc_cover) == 0:
        #    self.acc_cover = self.acc_cover_hd
        
        #pdb.set_trace()
        # check acc_cover
        if len(self.acc_cover) > 0:
            
            if isServerSide==False:
                #import pylab as pl
                acc_cover_log10 = np.log10(self.acc_cover)
                pl.hist(acc_cover_log10, bins=200, normed=True, cumulative=True)
                pl.gca().set_xlabel('log10(acc_cover)')
                pl.gca().set_ylabel('cdf')
                #pl.show()
                
                pl.gcf().savefig(Default_Ref_Path+'/sim_stat_dmp_fig_cdf_of_acc_cover_log10.png')
                #acc_cover_qt_log10 = np.percentile(acc_cover_log10, q=[50, 60, 70, 80, 90, 100])
                #pdb.set_trace()
        
        # check qt and acc_cover_qt
        if len(self.qt)>0 and len(self.acc_cover_qt)>0:
            self.dmp_vec(self.qt, 'quantile (%)')
            self.dmp_vec(self.acc_cover_qt, 'acc cover wrt quantile(%)')
        else:
            self.dmp_vec([], 'quantile (%)')
            self.dmp_vec([], 'acc cover wrt quantile(%)')
        
        #pdb.set_trace()
        
    """    

        
        
        
