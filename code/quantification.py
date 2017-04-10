
import sys, pdb

from tool_address import rsemPath, STAR_path
from old_code.util import run_cmd
from old_code.ReadProcess import BED2Exon1, RSEM2Coverage1

def main():

    pdb.set_trace()

    return

#modified from ReadProcess.py
#SE reads
#external
#estimate abundance level using RSEM (http://deweylab.biostat.wisc.edu/rsem/README.html)
#python quantification.py --RSEM1 --ref fa -r reads.fq -g gtf -O rsemOutDir -p rsemPrefixName

def RSEM1(args): #(ref_address, readFQ_address, BED_address, gtf_address):

    #pdb.set_trace()

    ref_address = args[args.index('--ref')+1]
    readFQ_address = args[args.index('-r')+1]
    gtf_address =args[args.index('-g')+1]
    rsemOutDir = args[args.index('-O')+1]
    rsemPrefixName = args[args.index('-p')+1]
    name = '%s/%s'%(rsemOutDir, rsemPrefixName)
    
    cmd = 'mkdir -p %s'%rsemOutDir
    run_cmd(cmd)

    rsem_index_command = rsemPath + '/rsem-prepare-reference'
    rsem_estim_command = rsemPath + '/rsem-calculate-expression'
   
    #cmd = rsem_index_command + ' --star -p 20 --star-path %s --gtf '%STAR_path + gtf_address + ' ' + ref_address + ' '  + name
    cmd = rsem_index_command + ' --star -p 20 --star-path %s --gtf '%STAR_path + gtf_address  + ' ' + ref_address + ' '  + name
    run_cmd(cmd)

    cmd = rsem_estim_command + ' --star -p 20 --star-path %s '%STAR_path + readFQ_address+ ' ' + name + ' ' + name
    run_cmd(cmd)

    RSEM_result_address = name + '.isoforms.results'

    return RSEM_result_address 

#generated rsem coverage based on rsem quantification results, output: exon_address (by product), and estimated rsem Coverage
#
#python quantification.py --rsemCoverage -i rsemRes (.isoforms.results) -b bedFile -e exon_address -c rsem coverage -L readLength
def rsemCoverage(args):

    #pdb.set_trace()

    rsemRes = args[args.index('-i')+1]
    bedFile = args[args.index('-b')+1]
    exonFile = args[args.index('-e')+1]
    rsemCov = args[args.index('-c')+1]
    readLength = int(args[args.index('-L')+1])

    #bed to exonFile
    BED2Exon1(bedFile, exonFile)
    print('%s written'%exonFile)

    #estimated coverage by rsem
    RSEM2Coverage1(rsemRes, exonFile,  rsemCov, readLength)
    print('%s written'%rsemCov)

    return

'''
python quantification.py --tbam2gbam -p rsemGenomePrefix (e.g. path/to/Chr15) -t tbam -g outputGbam (and Gsam, sorted by n)
example output:
Chr15.genome.bam, Chr15.genome.sorted_n.sam
'''
def tbam2gbam(args):

    rsemGenomePrefix = args[args.index('-p')+1]
    tbam = args[args.index('-t')+1]
    gbam = args[args.index('-g')+1]
    gsam = gbam[:-3]+'sorted_n.sam' #replace bam with sorted_n.sam

    #pdb.set_trace()

    cmd = '%s/rsem-tbam2gbam %s %s %s'%(rsemPath, rsemGenomePrefix, tbam , gbam)
    run_cmd(cmd)

    cmd = 'samtools sort -n -o %s -@ 20 %s'%(gsam, gbam)
    run_cmd(cmd)

    return

'''
usage: 

#quantify rna-seq onto transcriptome, output: rsemOutDir/rsemPrefixName.isoforms.results
#
python quantification.py --RSEM1 --ref fa -r reads.fq -g gtf -O rsemOutDir -p rsemPrefixName (e.g. chr15)

#generate rsem coverage based on rsem quantification results, output: exon_address (by product), and estimated rsem Coverage
#
python quantification.py --rsemCoverage -i rsemRes (.isoforms.results) -b bedFile -e exon_address -c rsem coverage -L readLength

#generate genome based read alignment based on rsem quantification results
#example output: Chr15.genome.bam, Chr15.genome.sorted_n.sam
#
python quantification.py --tbam2gbam -p rsemGenomePrefix (e.g. path/to/Chr15) -t tbam -g outputGbam (and Gsam, sorted by n)

'''
if __name__ == "__main__":

    args = sys.argv

    if '--RSEM1' in args:# reads --> transcriptome quantification
        RSEM1(args)
    elif '--rsemCoverage' in args:
        rsemCoverage(args)
    elif '--tbam2gbam' in args:
        tbam2gbam(args)
    else:
        main()