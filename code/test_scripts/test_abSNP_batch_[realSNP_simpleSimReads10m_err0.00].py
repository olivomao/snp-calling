import pdb
from old_code.util import run_cmd

#use prev data (e.g. fa, snps, exp, target alleles), to generate new reads and do snp calls by abSNP code

def main():

    #snp1k_reads10m_abSNP_batch
    #'''
    srcDir = '/data1/shunfu1/SNPCalling/data_real/'
    dstDir = '/data1/shunfu1/SNPCalling/data_real_abSNP_batch/'
    numReads = 10000000
    L = 100
    errRate = 0.00
    readsLabel = 'reads_N%s_L%s_Err%.2f'%(numReads, L, errRate)
    #'''
    
    chrom = 'chr15'

    refGenome = '%s/genome/%s.fa'%(srcDir, chrom)

    genomeFile = ['%s/ref_snp/chr15_[from_snp_m].fa'%srcDir,
                  '%s/ref_snp/chr15_[from_snp_p].fa'%srcDir]

    trBedFile = '%s/gene_annotation/gencode.v19.annotation.target.chr15.sorted.bed'%srcDir
    gtfFile = '%s/gene_annotation/gencode.v19.annotation.gtf'%srcDir     

    read_out_dir = '%s/%s/'%(dstDir, readsLabel)
    rsem_out_dir = '%s/rsem_%s/'%(dstDir, readsLabel)

    exon_address = '%s/exon.txt'%rsem_out_dir #intermediate output
    rsemCov = '%s/rsemCoverage.txt'%rsem_out_dir #intermediate output   

    rsemPrefixName = chrom
    rsemGenomePrefix = '%s/%s'%(rsem_out_dir, rsemPrefixName)
    tbam = '%s.transcript.bam'%rsemGenomePrefix
    gbam = '%s.genome.bam'%rsemGenomePrefix
    gsam = '%s.genome.sorted_n.sam'%rsemGenomePrefix    

    suf = ['_m', '_p']
    exp_path = ['%s/intermediate/exp_m.txt'%read_out_dir,  #intermediate output
                '%s/intermediate/exp_p.txt'%read_out_dir]  #intermediate output

    pdb.set_trace()    

    #read generation
    for i in range(2):

        cmd = 'python sim_data_generator.py --read_generation '+\
              '-g %s -b %s -O %s '%(genomeFile[i], trBedFile, read_out_dir)+\
              '-n %d -l %d -r %f '%(numReads, L, errRate)+\
              '--suffix %s --toFq '%suf[i] #+'--exp_path %s'%(exp_path[i])
        #pdb.set_trace()
        run_cmd(cmd)

    #merge reads
    cmd = 'python sim_data_generator.py --merge_reads -m %s/reads_m.fq -p %s/reads_p.fq -o %s/merged_reads.fq --uniqID'% \
          (read_out_dir, read_out_dir, read_out_dir)
    #pdb.set_trace()
    run_cmd(cmd)

    #quantification

    cmd = 'python quantification.py --RSEM1 --ref %s -r %s/merged_reads.fq -g %s -O %s -p %s'% \
          (refGenome, read_out_dir, gtfFile, rsem_out_dir, rsemPrefixName)
    #pdb.set_trace()
    run_cmd(cmd)

    rsemIsoformsResults = '%s/%s.isoforms.results'%(rsem_out_dir, rsemPrefixName)

    cmd = 'python quantification.py --rsemCoverage -i %s -b %s -e %s -c %s -L %d'% \
          (rsemIsoformsResults, trBedFile, exon_address, rsemCov, L)
    #pdb.set_trace()
    run_cmd(cmd)

    cmd = 'python quantification.py --tbam2gbam -p %s -t %s -g %s'% \
          (rsemGenomePrefix, tbam, gbam)
    #pdb.set_trace()
    run_cmd(cmd)

    #snp call - tbd

    return

if __name__ == "__main__":

    main()