# -*- coding: utf-8 -*-
"""
Created on Sun Oct  4 16:27:23 2015

@author: olivo
"""

import subprocess

cross_check=[1]
Threshold_num_reads=[1, 100, 200, 300]

for j in cross_check:
    for i in Threshold_num_reads:
        cmd = 'python batch_run_parallel.py '+\
              '%d '%j +\
              '%d '%i
              
        subprocess.call(cmd, shell=True)
