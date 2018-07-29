#PBS -S /bin/bash
#PBS -q batch
#PBS -o /lustre1/jl03308/20xIllumina-lambda/analysis/jobs/
#PBS -e /lustre1/jl03308/20xIllumina-lambda/analysis/jobs/
#PBS -l nodes=1:ppn=12:HIGHMEM
#PBS -l mem=200gb
#PBS -l walltime=200:00:00
#PBS -M jl03308@uga.edu
#PBS -m ae
#PBS -N bwa_mem_1211-1_S2_L001
module load bwa/0.7.15
module load samtools/1.3.1
module load bedtools/2.24.0
cd /lustre1/jl03308/20xIllumina-lambda/analysis/lambda_junction-maize
bwa mem -t 12 /lustre1/jl03308/reference/source/maize_genome_v4/Zea_mays.AGPv4.dna.toplevel.fa /lustre1/jl03308/20xIllumina-lambda/data/1211-trimmed/1211-1_S2_L001_R1_001_val_1.fq.gz /lustre1/jl03308/20xIllumina-lambda/data/1211-trimmed/1211-1_S2_L001_R2_001_val_2.fq.gz > 1211-1_S2_L001_maize.sam
samtools view -b -o 1211-1_S2_L001_maize.bam 1211-1_S2_L001_maize.sam

samtools view -b -f4 1211-1_S2_L001_maize.bam > 1211-1_S2_L001_unmapped.bam
samtools sort 1211-1_S2_L001_unmapped.bam > 1211-1_S2_L001_unmapped.sorted.bam
bedtools bamtofastq -i 1211-1_S2_L001_unmapped.sorted.bam -fq 1211-1_S2_L001_unmapped.fq

bwa mem -t 12 /lustre1/jl03308/reference/source/lambda_genome/lambdax2.fa 1211-1_S2_L001_unmapped.fq > 1211-1_S2_L001_lambdaX2.sam
samtools view -bhq 20 -o 1211-1_S2_L001_lambdaX2_q20.bam 1211-1_S2_L001_lambdaX2.sam

