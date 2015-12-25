from count import *
from write_correlated_regions import *

def test01():
    reads = "data/test_01/reads.sam"
    reference = "data/G_30000000-1/Chr15.fa"
    coverage = "data/G_30000000-1/coverage.txt"
    generate_count_file(reads, reference, coverage)
    correlations = "output/correlations.pickle"
    countfile = "output/count.txt"
    outfile = "output/correlated_regions.txt"
    write_correlations(countfile, correlations, outfile)

def test02():
    reads = "data/test_02/mini_sorted.sam"
    reference = "data/test_02/hg19_chr15.fa"
    coverage = "data/test_02/rsem_coverage.txt"
    generate_count_file(reads, reference, coverage)
    correlations = "output/correlations.pickle"
    countfile = "output/count.txt"
    outfile = "output/correlated_regions.txt"
    write_correlations(countfile, correlations, outfile)

if __name__ == "__main__":
    test02()