#!/usr/bin/python -tt


import subprocess
import re
import numpy as np
from numpy import *
import pdb


def index(BowtiePath, Ref_address, aligner):
 
    #pdb.set_trace()

    if (aligner == 'bowtie2') or (aligner =='tophat'):
        #if ServerNum==2:
        #indexCommand = 'bowtie-build' # shunfu '/bowtie2-build'
        #else:
        indexCommand = 'bowtie2-build' # shunfu '/bowtie2-build'
    elif aligner == 'bowtie':
        indexCommand = 'bowtie-build' #'/bowtie-build'
    
    b=re.search('([-/\w\s]+)/([-\w.]+)\.([-\w.]+)', Ref_address)
    Path = b.group(1)
    Index_address = Path + '/index/' + b.group(2)

    IndexProgram = BowtiePath + indexCommand
    #subprocess.call([IndexProgram , Ref_address , Index_address] )
    subprocess.call(IndexProgram + ' ' + Ref_address + ' ' + Index_address, shell = True)
    return Index_address
    
def Bowtie(BowtiePath, Index_address, Read_address, Map_address, args, aligner):
    if aligner == 'bowtie':
        alignCommand = '/bowtie'
    elif aligner == 'bowtie2':
        alignCommand = '/bowtie2'
    
    AlignProgram = BowtiePath + alignCommand
    subprocess.call(AlignProgram + ' ' + args + ' ' + Index_address + ' ' +  Read_address + ' ' + Map_address, shell=True )
    
def SamRefIndex(SamPath, Ref_address):
    IndexProgram = SamPath + '/samtools'
    subprocess.call([IndexProgram , 'faidx', Ref_address] )
    
def SamAlignIndex(SamPath, sam_address):
    SamProgram =  SamPath + '/samtools'
    b = re.search('([/\s\w.-]+).sam', sam_address)
    bam_address = b.group(1) + '.bam' 
    
    subprocess.call( SamProgram + ' view -bS ' + sam_address + ' > ' + bam_address , shell=True )
    

    
    sorted_address = b.group(1) + '.sorted' 
    subprocess.call(SamProgram + ' sort ' + bam_address + ' ' + sorted_address, shell=True )

    
    sorted_address = bam_address = b.group(1) + '.sorted.bam'
    subprocess.call([SamProgram , 'index' , sorted_address] )
    return sorted_address
    
    
def TopHat(TopHatPath, Index_address, Read_address, out_address, args):
    #pdb.set_trace() #debug
    TopHatProgram = TopHatPath + '/tophat' 
    subprocess.call(TopHatProgram +  ' -o ' + out_address + args + Index_address + ' ' + Read_address , shell=True )
