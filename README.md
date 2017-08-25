# Topics

<a href='#intro'> Introduction </a>

<a href='#cite'> Citation </a>

<a href='#install'> Installation </a>

<a href='#demo'> Demo </a>

<a> Usages </a>

- <a href='#snp_sim'> Usage 1: SNP Simulator </a>

- <a href='#read_sim'> Usage 2: Read Simulator </a>

- <a href='#snp_call'> Usage 3: SNP Calling </a>

- <a href='#eval'> Usage 4: Evaluation </a>

<a href='#history'> History </a>

---

# Introduction <a id='intro'></a>

abSNP is an RNA-Seq SNP calling software. It takes raw reads or read alignment as input, and output SNP candidates. Unlike existing SNP callers, it explores abundance estimation and is able to detect SNPs at repetitive genomic regions.

This software is developed by [Shunfu Mao](shunfu@uw.edu) and [Sreeram Kannan](ksreeram@uw.edu) at University of Washington (UW), and [Soheil Mohajer](soheil@umn.edu) at University of Minnesota, and is currently maintained at UW. 

---

# Citation <a id='cite'></a>

If you would like to cite abSNP, please find more information [here](http://dx.doi.org/10.4230/LIPIcs.WABI.2017.15).

---

# Work Flow

The software of abSNP mainly contains three parts:

- A SNP simulator which generates random SNPs and target maternal and paternal alleles from a reference genome.

- A read simulator which samples reads from target alleles containing SNP information.

- A SNP Caller which calls SNPs based on raw reads and output SNP candidates.

We also provide a simple and quick check of performance in terms of mis-detection and false positive.

---

# Installation <a id='install'></a>

## System Requirement

Unix system (e.g. Ubuntu).

We have not tested on Mac system but it should work if you change STAR_path in tool_address.py for Mac option.

## Steps

1) At terminal, create a folder to download codes:

~~~
mkdir -p "<path/to/abSNP>"
~~~

2)  Go to this folder as your current dir and download code. You'll find a folder named 'abSNP' at your current dir.

~~~
cd <path/to/abSNP>

git clone "https://github.com/shunfumao/abSNP"
~~~

3) Unzip necessary tools

~~~
cd "<path/to/abSNP>/abSNP/code/tools/"

unzip tools.zip
~~~

Tools include:

tools | usage
---- | ----
FA2FQ | used in read simulator
RNASeqReadSimulator | used in snp simulator and read simulator
rsem | used in SNP calling
STAR | used in SNP calling
GATK, picard | not required here. They're only used to run GATK best practice for performance comparison purpose. 


4) Dependent Python packages.
In addition to the tools installed at 3), we also use some Python packages.

Packages include:

packages | how to obtain
---- | ----
intervaltree | pip install intervaltree

---

# Demo <a id='demo'></a>

To check the installation or to see how the abSNP is used, you can go to the main code folder as your current working directory, and run dem. For example:

~~~
cd "</path/to/abSNP>/abSNP/code/"

sh demo.sh
~~~

The whole process takes about 8 ~ 10 minutes to run on a lab server using 1 CPU core.

Initially we have a reference genome (e.g. Chr15.fa) and transcriptome annotations (e.g. hg19_chr15-UCSC.gtf, and hg19_chr15-UCSC-sorted.bed). Based on these, 

- simulated snps as well as target genomes (both maternal and paternal) will be generated and stored at demo/snp_dir/. 

- Then reads are sampled from the alleles and stored at demo/reads/.

- Finally, based on reads and reference genome and annotated transcriptome, SNPs are called and stored at demo/res/.

Simple statistics will be printed on screen by comparing the ground truth SNPs in demo/snp_dir/ with our called SNPs in demo/res/.

---

# Usage 1 : SNP Simulator <a id='snp_sim'></a>

Based on a reference transcriptome (a set of RNA transcripts), the simulator assigns rand expression levels for every RNA transcript independently for each maternal and paternal allele. Then random snps restricted to these RNA transcripts are generated for each allele. The genomes of target alleles that contain SNPs are also generated.

## Command

~~~
python run_abSNP.py --snp_simulator
                    --bedSorted <bedSorted>
                    --outDir <outDir>
                    --qt <quantile>
                    --refGenome <refGenomeFile>
                    --NumSNP <numSNP>
~~~

## Input

argument | description
---- | ----
\--bedSorted <bedSorted\> | path of a sorted bed file (e.g. a reference transcriptome) to specify where SNPs shall be generated. To sort a bed file (e.g by its chrom & genome start pos), you can use e.g. "sort -k 1,1 -k 2,2n in.BED > out.BED".
\--outDir <outDir\> |  path of the directory to store result files, including:  exp_m.txt, exp_p.txt, SNP_m.txt, SNP_p.txt, Tar_m.txt, Tar_p.txt. m, p stand for maternal & paternal alleles
\--qt <quantile\> | an integer value in [0,99]; e.g. \--qt 90 means SNPs are generated in top (100-90)% highly expressed genes
\--refGenome <refGenomeFile\> | path of the reference genome containing one chromosome
\--NumSNP <numSNP\> | an integer value in [1, inf) to specify the number of SNPs generated per allele

## Output

files | description
---- | ----
exp_m.txt, exp_p.txt | expression levels of RNA transcript for maternal  & paternal alleles 
SNP_m.txt, SNP_p.txt | simulated snps, in format of col-0 as 0-based locus, col-1 as reference base, col-2 as '\-->', and col-3 as target/mutated base.
Tar_m.txt, Tar_p.txt | genomes (one chromosome) of target alleles that contain the simulated SNPs. in fasta format.

---

# Usage 2 : Read Simulator <a id='read_sim'></a>

Sample reads from the transcriptome of target maternal and paternal alleles. The amount of reads sampled from per RNA transcript is specified by expression files. The number of reads per allele, the read length and error rate are also required.

## Command
~~~
python run_abSNP.py --read_simulator
                    --tarGenome_m <tarGenome_m>
                    --tarGenome_p <tarGenome_p>
                    --exp_m <exp_m>
                    --exp_p <exp_p>
                    --bedSorted <bedSorted>
                    --outDir <outDir>
                    --numReads <numReads>
                    --readLength <readLength>
                    --errRate <errRate>

~~~

## Input
argument | description
---- | ----
\--tarGenome_m <tarGenome_m\> | path of the target maternal genome, containing one chromosome, in fasta format. Generated from SNP simulator.
\--tarGenome_p <tarGenome_p\> | path of the target paternal genome, containing one chromosome, in fasta format. Generated from SNP simulator.
\--exp_m <exp_m\> | path of the expression level file of RNA transcript for maternal allele. Generated from SNP simulator.
\--exp_p <exp_p\> | path of the expression level file of RNA transcript for paternal allele. Generated from SNP simulator.
\--bedSorted <bedSorted\> | path of the reference transcriptome in sorted bed format. To be consistent with the one used in SNP simulator.
\--outDir <outDir\> | path of the directory to store the generated simulated reads in fastq format (e.g. merged_reads.fq).
\--numReads <numReads\> | integer value. Specify the number of reads to be generated per target allele.
\--readLength <readLength\> | integer value. Specify the read length (e.g. 100).
\--errRate <errRate\> | float value. Specify the read error rate (e.g. 0.01)

## Output

files | description
---- | ----
merged_reads.fq | the generated simulated reads in fastq format, including reads sampled from both alleles.

---

# Usage 3 : SNP Calling <a id='snp_call'></a>

## Command

~~~
python run_abSNP.py --call_snps
                    --refGenome <refGenome>
                    --readFile <readFile>
                    --gtfFile <gtfFile>
                    --bedSorted <bedSorted>
                    --p <numProcess>
                    --readLength <readLength>
                    --alpha <alpha>
                    --outDir <outDir>
                                                
~~~

## Input

argument | description
---- | ----
\--refGenome <refGenome\> | path of the reference genome containing one chromosome
\--readFile <readFile\>  | path of the sampled read file in fastq format
\--gtfFile <gtfFile\> | path of the reference transcriptome (set of RNA transcripts, which are considered as target regions to detect SNPs) in gtf format
\--bedSorted <bedSorted\> | path of the reference transcriptome in sorted bed format. To be consistent with the gtfFile.
\--p <numProcess\> | an integer value specifyng how many CPU processes to use. Default 1. Recommend to use 10 to 20.
\--readLength <readLength\> | integer value. Specify the read length (e.g. 100).
\--alpha <alpha\> | a float value between 0 and 1. 0 corresponds to minimum false positive and 1 corresponds to max sensitivity.
\--outDir <outDir\> | ath of the directory to store the file of called SNP candidates.

## Output

files | description
---- | ----
snp_candidates.txt | SNPs called by abSNP. format, with tab separated columns. col-0: 0-based genome locus; col-1: reference base; col-2: '\-->' (a dummy symbol); col-3: target base

---

# Usage 4 : Evaluation <a id='eval'></a>

We provde a simple and quick check of SNP calling performance in terms of mis-detection and false positive.

## Command

~~~
python run_abSNP.py --check
                    --snp_res <snpRes>
		    --snp_m <snpM>
		    --snp_p <snpP>                     
                                                
~~~

## Input

argument | description
---- | ----
\--snp_res <snpRes\> |  path of SNP candidates called by abSNP
\--snp_m <snpM\> | path of ground truth SNPs from maternal allele of the target individual
\--snp_p <snpP\> | path of ground truth SNPs from the paternal allele of the target individual

## Output

Statistics of mis-detection and false positives will be printed on screen.

---


# History<a id='history'></a>

2017.5.9 initial manual created
