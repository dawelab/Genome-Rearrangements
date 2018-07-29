#PBS -S /bin/bash
#PBS -q batch
#PBS -o /lustre1/jl03308/20xIllumina-control/analysis/jobs/
#PBS -e /lustre1/jl03308/20xIllumina-control/analysis/jobs/
#PBS -l nodes=1:ppn=12:HIGHMEM
#PBS -l mem=20gb
#PBS -l walltime=10:00:00
#PBS -M jl03308@uga.edu
#PBS -m ae
#PBS -N bwa_mem_1419-10_1
module load bwa/0.7.15
module load samtools/1.3.1
cd /lustre1/jl03308/20xIllumina-control/analysis/10A
mkdir 1419-10
cd 1419-10
bwa mem -t 12 /lustre1/jl03308/reference/concatenated_genome/rice_10A.fa /lustre1/jl03308/20xIllumina-control/data/trimmed/10A/1419-10_S5_L001_R1_001_val_1.fq.gz /lustre1/jl03308/20xIllumina-control/data/trimmed/10A/1419-10_S5_L001_R2_001_val_2.fq.gz > 1419-10_1.sam
samtools view -@ 12 -b -o 1419-10_1.bam 1419-10_1.sam

