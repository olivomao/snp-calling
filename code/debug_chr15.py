import count
from write_correlated_regions import write_correlations

prefix = "/data/soheil/Chr15/"

# reads = prefix + "Bowtie/tophat-GTF/with-T-paired/accepted_hits.sam"
reads = "data/accepted_hits_sorted.sam"
reference = prefix + "hg19_chr15.fa"
coverage = prefix + "rsem_coverage.txt"
count.generate_count_file(reads, reference, coverage)

correlations = "output/correlations.pickle"
countfile = "output/count.txt"
outfile = "output/correlated_regions.txt"
write_correlations(countfile, correlations, outfile)