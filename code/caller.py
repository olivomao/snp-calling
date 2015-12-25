import csv
import operator
import pickle

acgt = ["A", "C", "G", "T"]
correlations = pickle.load(open("output/correlations.txt"))

def max_caller(N, x):
    max_index, max_value = max(enumerate(N), key=operator.itemgetter(1))
    return acgt[max_index]

def exp_caller(N, x):
    pass

def run_caller(input, call):
    with open(input) as counts:
        reader = csv.reader(counts, delimiter='\t')

        data = {}
        for row in reader:
            if len(row) >= 8:
                i = int(row[0]) # exon position
                j = int(row[1]) # ref position
                r = row[2].upper() # ref base
                x = float(row[3]) # expression level
                N = map(int, row[4:8]) # counts
                c = call(N, x)
                data[j] = (i, j, r, x, N, c, correlations[j] if j in correlations else None)

    with open(input) as counts:
        reader = csv.reader(counts, delimiter='\t')
        s = 0
        for row in reader:
            if len(row) >= 8:
                i = int(row[0]) # exon position
                j = int(row[1]) # ref position
                r = row[2].upper() # ref base
                x = float(row[3]) # expression level
                N = map(int, row[4:8]) # counts
                c = call(N, x)
                if sum(N) > 10 and j in correlations: # repeat regions with enough info
                    c = call(N, x)
                    y = sorted(N)[::-1]
                    if y[0] < 5 * y[1]:
                        s = s + 1
                        print "-----------------------"
                        print i, j, r, x, N, c
                        for d in correlations[j]:
                            if d in data:
                                print data[d]
                            else:
                                print "missing!"
        print s

if __name__ == "__main__":
    run_caller("output/count.txt", max_caller)
