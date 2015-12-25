import csv
import operator
import time
import math

ACGT = ["A", "C", "G", "T"]
BASES = {"A": 0, "C": 1, "G": 2, "T": 3}
NUM_TO_BASE = {0: "A", 1: "C", 2: "G", 3: "T"}

ERROR = 0.1

def quality(character):
    """Using SAM specs: https://samtools.github.io/hts-specs/SAMv1.pdf
    Maybe this should be the Pr{base is wrong}? Some exp function."""
    Q = ord(character) - 33
    return 1.0 - (10.0 ** (-Q / 10.0))

def caller(i, ref, l_0, counts, reads):
    """Follows the paper (note-4.pdf) closely"""
    # counts_sorted, bases = zip(*sorted(zip(counts, ACGT))[::-1])
    # majority_base = bases[0]

    # # of_interest: boolean, is this an interesting call to consider?
    # of_interest = counts_sorted[1] > ERROR * counts_sorted[0]

    # if not of_interest:
    #     return majority_base, False, False

    # # this call is of_interest
    # if counts_sorted[0] < exp * (1 - ERROR):
    #     return bases[1], True, True
    
    # return majority_base, False, False



    # all variables from paper
    b, q, ly, lsum = map(list, zip(*reads))
    # todo: r is missing in description
    q = map(quality, q)
    ly = map(float, ly)
    lsum = map(float, lsum)

    # filter reads for deletions
    j = 0
    while j < len(b):
        if b[j] == 'D':
            del b[j]
            del q[j]
            del ly[j]
            del lsum[j]
            j -= 1
        j += 1

    n = len(b) # == sum(counts) # should we call it N?

    # b_j is the same as y?
    def P(x, j):
        """Probability {read j = y | X = x}"""
        y = b[j]
        if y == x and x == ref:
            z = l_0 * q[j] + ly[j]
        elif y == x and x != ref:
            z = l_0 * q[j] / 2 + ly[j]
        elif y == ref and ref != x:
            z = l_0 * q[j] / 2 + ly[j]
        else:
            z = l_0 * (1 - q[j]) / 3 + ly[j]
        s = lsum[j] + 10**(-10)
        return z / s if s else 0

    probabilities = []
    for base in BASES:
        p = 0
        for j in range(n):
            p += math.log(P(base, j))
        probabilities.append((p, base))

    return max(probabilities)[1] + ref

def split(read_base):
    """Split up a string like T,I,1,1"""
    x = read_base.split(",")

    # special case when the quality score is a comma
    if len(x) == 5:
        del x[1]
        x[1] = ","

    return x

def decode(row):
    global_pos = int(row[1])
    ref = row[2].upper()
    exp = float(row[3])

    counts = map(int, row[4:8])
    reads = map(split, row[8:])

    return global_pos, ref, exp, counts, reads

def main():
    start_time = time.clock()

    calls = open("output/calls.txt", "w+")
    with open("output/count.txt") as count_file:
        reader = csv.reader(count_file, delimiter='\t')

        progress = 0
        for row in reader:
            progress += 1
            if progress % 100000 == 0:
                print("positions called:", str(progress),
                            round(time.clock() - start_time, 2))

            if len(row) <= 5:
                continue
            
            call = caller(*decode(row))
            info = [row[1]] + row[4:8] + [call]
            outrow = "\t".join(info) + "\n"
            calls.write(outrow)

if __name__ == "__main__":
    main()