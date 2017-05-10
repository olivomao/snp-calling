demo_dir=demo

snp_dir=$demo_dir/snp_dir/
read_dir=$demo_dir/reads/
res_dir=$demo_dir/res/
num_cpu=1

python run_abSNP.py --snp_simulator \
                    --bedSorted $demo_dir/hg19_chr15-UCSC-sorted.bed \
                    --outDir $snp_dir \
                    --qt 99 \
                    --refGenome $demo_dir/Chr15.fa \
                    --NumSNP 100

python run_abSNP.py --read_simulator \
                    --tarGenome_m $snp_dir/Tar_m.txt \
                    --tarGenome_p $snp_dir/Tar_p.txt \
                    --exp_m $snp_dir/exp_m.txt \
                    --exp_p $snp_dir/exp_p.txt \
                    --bedSorted $demo_dir/hg19_chr15-UCSC-sorted.bed \
                    --outDir $read_dir/ \
                    --numReads 10000 \
                    --readLength 100 \
                    --errRate 0.00

python run_abSNP.py --call_snps \
                    --refGenome $demo_dir/Chr15.fa \
                    --readFile $read_dir/merged_reads.fq \
                    --gtfFile $demo_dir/hg19_chr15-UCSC.gtf \
                    --bedSorted $demo_dir/hg19_chr15-UCSC-sorted.bed \
                    --readLength 100 \
                    --alpha 0.0 \
                    --outDir $res_dir/ \
                    --p $num_cpu

python run_abSNP.py --check \
                    --snp_res $res_dir/snp_candidates.txt \
                    --snp_m $snp_dir/SNP_m.txt \
                    --snp_p $snp_dir/SNP_p.txt