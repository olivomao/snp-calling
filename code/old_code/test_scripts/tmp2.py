try:
    from itertools import izip as zip
except ImportError: # will be 3.x series
    pass

def merge_coverage(coverage_address_m, coverage_address_p):

    coverage_address = coverage_address_m[:-6]+'.txt' #e.g. coverage_m.txt --> coverage.txt

    with open(coverage_address_m, 'r') as cm, open(coverage_address_p, 'r') as cp, open(coverage_address, 'w') as c:

        for mline, pline in zip(cm, cp):

            if mline=='' or pline=='':
                continue

            mtoks = mline.split()
            ptoks = pline.split()
            cov_sum = float(mtoks[4])+float(ptoks[4])

            newline = '\t'.join(mtoks[0:4])+'\t%f\n'%cov_sum

            c.write(newline)    

    return coverage_address

merge_coverage('tmp/merge_coverage/cov_m.txt', 'tmp/merge_coverage/cov_p.txt')
