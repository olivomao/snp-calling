import subprocess, re, time, sys, os, shutil, pdb
import multiprocessing

def run_cmd(cmd, shell=True):
    #print(cmd)
    subprocess.call(cmd, shell=shell)

def run_cmds(cmds,noJobs=20):

    if noJobs==1:
        #pdb.set_trace()
        for cmd in cmds:
            run_cmd(cmd)
    else:
        #cmds is a tuple of strings
        #noJobs is a an integer with the no of Jobs to run in parallel
        p = multiprocessing.Pool(noJobs)
        p.map(run_cmd,cmds)

'''
read fasta file and output sequences (seq_name & seq)

copied from refS util
'''
def from_fasta(filename):
    sequences = {} #[]
    seq = []
    next_name = '_'
    with open(filename) as f:
        for line in f:
            if line[0] == '>':
                #sequences.append((next_name, ''.join(seq)))
                sequences[next_name]=''.join(seq)
                seq = []
                next_name = line.split()[0][1:]
            else:
                seq.append(line.strip().upper())
    #sequences.append((next_name, ''.join(seq)))
    sequences[next_name]=''.join(seq)
    del sequences['_']
    #pdb.set_trace()
    return sequences #[1:]

def to_fasta(filename, seq, seq_name):
    seqs_to_fasta(filename, [(seq_name, seq)])
    print('%s written'%filename)

def seqs_to_fasta(filename, seqs):
    LINE_LENGTH = 50
    with open(filename, 'w') as f:
        for seq_name, seq in seqs:
            f.write(">{}\n".format(seq_name))
            for i in xrange(0, len(seq), LINE_LENGTH):
                f.write("{}\n".format(seq[i:i+LINE_LENGTH]))

def parent_dir(dir_path):
    dir_path = dir_path.split('/')
    dir_path = [itm for itm in dir_path if itm != '']
    return '/'+'/'.join(dir_path[:-1])+'/'