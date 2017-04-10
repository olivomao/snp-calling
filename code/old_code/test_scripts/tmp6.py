#utilize rsem bam
#
#rsem-tbam2gbam Chr15 Chr15.transcript.bam Chr15.genome.bam -p 20
#samtools sort -n -o Chr15.genome.sorted_n.sam -@ 20 Chr15.genome.bam