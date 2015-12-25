from itertools import izip
import re

# precompile regular expression patterns
num_pattern = re.compile('[\d]+')
chr_pattern = re.compile('[^\d]+')
dna_pattern = re.compile(r'([acgtnACGTN]+)')
BASES = {"A": 0, "C": 1, "G": 2, "T": 3}
NUM_TO_BASE = {0: "A", 1: "C", 2: "G", 3: "T"}

def load_exon(exon_address):
    with open(exon_address) as exon_file:
        return exon_file.readline().strip()

def dump_exon(exon_address, exon):
    with open(exon_address, 'w+') as exon_file:
        exon_file.write(exon)