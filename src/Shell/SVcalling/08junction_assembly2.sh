#PBS -S /bin/bash
#PBS -q batch
#PBS -N cluster-E24
#PBS -o /lustre1/jl03308/20xIllumina-lambda/analysis/jobs
#PBS -e /lustre1/jl03308/20xIllumina-lambda/analysis/jobs
#PBS -l nodes=1:ppn=12:HIGHMEM
#PBS -l mem=100gb
#PBS -l walltime=480:00:00
#PBS -M jl03308@uga.edu
#PBS -m ae

module load samtools/1.3.1
module load bedtools/2.24.0
module load spades/3.10.0
module load anaconda
cd /lustre1/jl03308/20xIllumina-lambda/analysis/RiceE24/supplement/lumpy
bedtools bamtofastq -i RiceE24.discordants.bam -fq RiceE24.discordants.fq
bedtools bamtofastq -i RiceE24.splitters.bam -fq RiceE24.splitters.fq
cat RiceE24.discordants.fq RiceE24.splitters.fq > RiceE24-sum.fq 
python3 /usr/local/apps/spades/3.10.0/bin/spades.py -s RiceE24-sum.fq  -o sum-junction
