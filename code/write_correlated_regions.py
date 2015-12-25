import csv
import operator
import pickle

def write_correlations(count_file, correlations_file, out_file):
    correlations = pickle.load(open(correlations_file))
    print("correlations loaded")

    out = open(out_file, "w+")

    with open(count_file) as counts:
        _counter = 0
        reader = csv.reader(counts, delimiter='\t')
        for row in reader:
            if len(row) >= 8:
                _counter += 1
                if _counter % 100000 == 0:
                    print(_counter, "lines processed!")

                j = int(row[1]) # ref position

                if j in correlations: # repeat regions with enough info
                    # i = int(row[0]) # exon position
                    # r = row[2].upper() # ref base
                    # x = float(row[3]) # expression level
                    
                    N = map(int, row[4:8]) # counts
                    y = sorted(N)[::-1]
                    if sum(N) > 10 and y[0] < 5 * y[1]:
                        # data = (i, j, r, x, N, correlations[j])
                        row.append(str(correlations[j]))
                        out.write("\t".join(row))
                        out.write("\n")

if __name__ == "__main__":
    working_dir = "/Users/Soheil/Dropbox/Transcriptome/SNP-Calling-Summer15/data_2/"
    correlations = working_dir + "output/correlations.pickle"
    countfile = working_dir + "output/debug_count.txt"
    outfile = working_dir + "output/correlated_regions.txt"
    write_correlations(countfile, correlations, outfile)