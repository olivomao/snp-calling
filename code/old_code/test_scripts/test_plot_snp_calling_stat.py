
'''
snpSum example:

O-md    O-fp    G-md    G-fp    (snps of uniq mapped reads)
11  0   12  0
O-md    O-fp    G-md    G-fp    (snps of NON uniq mapped reads)
0   0   4   0
O-md    O-fp    G-md    G-fp    (tot)
11  0   16  0

'''

import pdb
import matplotlib.pyplot as plt
import numpy as np

def extract_snpSum(file):

    with open(file, 'r') as f:

        f.readline()
        U = f.readline() #snps with reads uniq mapped statistics (e.g. O-md, O-fp, G-md, G-fp)
        tokens = [int(i) for i in U.split() if i != '']
        U_O = tokens[0:2] #our caller [md, fp]
        U_G = tokens[2:4] #GATK etc caller [md, fp]

        f.readline()
        N = f.readline() #snps with reads non uniq mapped statistics (e.g. O-md, O-fp, G-md, G-fp)
        tokens = [int(i) for i in N.split() if i != '']
        N_O = tokens[0:2] #our caller [md, fp]
        N_G = tokens[2:4] #GATK etc caller [md, fp]

    return [U_O, U_G, N_O, N_G]

if __name__ == "__main__":

    prefix = 'snp20_reads100k_'

    cases = list(range(12))+[i+28 for i in list(range(5))]

    #pdb.set_trace()

    U_O_mul = [] #list of md & fp
    U_G_mul = []
    N_O_mul = []
    N_G_mul = []

    for i in cases:

        f = 'tmp/%s%d/snpSum.txt'%(prefix, i)

        [U_O, U_G, N_O, N_G] = extract_snpSum(f)

        U_O_mul.append(U_O)
        U_G_mul.append(U_G)
        N_O_mul.append(N_O)
        N_G_mul.append(N_G)

    #pdb.set_trace()

    fig, axarr = plt.subplots(3,2)

    #md-fp curve

    ax = axarr[0,0]
    omd = [i[0] for i in U_O_mul]
    ofp = [i[1] for i in U_O_mul]
    gmd = [i[0] for i in U_G_mul]
    gfp = [i[1] for i in U_G_mul]

    ax.scatter(omd, ofp, marker='+', label='Our caller')
    ax.scatter(gmd, gfp, marker='o',label='GATK caller')
    ax.set_title('snps w/ only uniquely mapped reads')
    ax.set_xlim([-1,50])
    ax.set_ylim([-1,50])
    ax.set_xlabel('mis detection')
    ax.set_ylabel('false positive')
    ax.legend()
    ax.grid(True)

    #md histogram

    ax = axarr[1,0]
    ax.set_xlabel('mis detection')
    ax.set_ylabel('histogram')
    #[hist, bins]=np.histogram(omd, bins=10)
    ax.hist([omd, gmd])

    #fp histogram

    ax = axarr[2,0]
    ax.set_xlabel('false positive')
    ax.set_ylabel('histogram')
    #[hist, bins]=np.histogram(omd, bins=10)
    ax.hist([ofp, gfp])


    omd = [i[0] for i in N_O_mul]
    ofp = [i[1] for i in N_O_mul]
    gmd = [i[0] for i in N_G_mul]
    gfp = [i[1] for i in N_G_mul]

    ax = axarr[0,1]
    ax.scatter([i[0] for i in N_O_mul], [i[1] for i in N_O_mul], marker='+', label='Our caller')
    ax.scatter([i[0] for i in N_G_mul], [i[1] for i in N_G_mul], marker='o',label='GATK caller')
    ax.set_title('snps w/ only multiply mapped reads')
    ax.set_xlim([-1,50])
    ax.set_ylim([-1,50])
    ax.set_xlabel('mis detection')
    ax.set_ylabel('false positive')
    ax.legend()
    ax.grid(True)

    #md histogram

    ax = axarr[1,1]
    ax.set_xlabel('mis detection')
    ax.set_ylabel('histogram')
    #[hist, bins]=np.histogram(omd, bins=10)
    ax.hist([omd, gmd])

    #fp histogram

    ax = axarr[2,1]
    ax.set_xlabel('false positive')
    ax.set_ylabel('histogram')
    #[hist, bins]=np.histogram(omd, bins=10)
    ax.hist([ofp, gfp])



    plt.show()

