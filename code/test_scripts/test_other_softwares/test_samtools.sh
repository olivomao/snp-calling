ref_file=/data1/shunfu1/SNPCalling/data_large_0_idealCov/Chr15.fa
bam_file=/data1/shunfu1/SNPCalling/SimSNPs_MultiRun_5/reads_N100K_L100_Err0.00/rsem/Chr15.genome.sorted.bam
out_file=/data1/shunfu1/SNPCalling/SimSNPs_MultiRun_5/reads_N100K_L100_Err0.00/samtools/res.vcf.gz

samtools mpileup -f $ref_file $bam_file | bcftools call -o $out_file