#PBS -S /bin/bash
#PBS -q batch
#PBS -o /lustre1/jl03308/20xIllumina-control/analysis/jobs/
#PBS -e /lustre1/jl03308/20xIllumina-control/analysis/jobs/
#PBS -l nodes=1:ppn=4:HIGHMEM
#PBS -l mem=20gb
#PBS -l walltime=10:00:00
#PBS -M jl03308@uga.edu
#PBS -m ae
#PBS -N trimgalore_12A_1
module load python/2.7.8
module load cutadapt/1.9.dev1 fastqc/0.11.3 trimgalore/4.0
module load java/jdk1.8.0_20 fastqc
/usr/local/apps/trimgalore/4.0/trim_galore --fastqc --gzip --paired /lustre1/jl03308/20xIllumina-control/data/1419_12A/1419-8_S5_L001_R1_001.fastq.gz /lustre1/jl03308/20xIllumina-control/data/1419_12A/1419-8_S5_L001_R2_001.fastq.gz -o /lustre1/jl03308/20xIllumina-control/data/trimmed/12A

