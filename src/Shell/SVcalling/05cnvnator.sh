#PBS -S /bin/bash
#PBS -q batch
#PBS -o /lustre1/jl03308/20xIllumina-control/analysis/10A/jobs/
#PBS -e /lustre1/jl03308/20xIllumina-control/analysis/10A/jobs/
#PBS -l nodes=1:ppn=2:HIGHMEM
#PBS -l mem=20gb
#PBS -l walltime=1:00:00
#PBS -M jl03308@uga.edu
#PBS -m ae
#PBS -N cnvnator_1419-10
module load cnvnator/0.3.3-root5.34.36
module load anaconda3/5.0.0
cd /lustre1/jl03308/20xIllumina-control/analysis/10A/1419-10
#mkdir cnvnator
cd cnvnator
cnvnator -root 500bp.root -unique -chrom 1 2 3 4 5 6 7 8 9 10 11 12 -genome /lustre1/jl03308/reference/concatenated_genome/rice_10A.fa -tree ../1419-10.rmdup.q20.sorted.bam
cnvnator -root 500bp.root -chrom 1 2 3 4 5 6 7 8 9 10 11 12 -genome /lustre1/jl03308/reference/concatenated_genome/rice_10A.fa -his 500 -d /lustre1/jl03308/20xIllumina-control/analysis/cfg/chr_dna_10A
cnvnator -root 500bp.root -chrom 1 2 3 4 5 6 7 8 9 10 11 12 -genome /lustre1/jl03308/reference/concatenated_genome/rice_10A.fa -stat 500
cnvnator -root 500bp.root -chrom 1 2 3 4 5 6 7 8 9 10 11 12 -genome /lustre1/jl03308/reference/concatenated_genome/rice_10A.fa -partition 500
cnvnator -root 500bp.root -chrom 1 2 3 4 5 6 7 8 9 10 11 12 -genome /lustre1/jl03308/reference/concatenated_genome/rice_10A.fa -call 500 > 1419-10.500bp.cnvnator
awk '{ print $2 } END { print "exit" }' 1419-10.500bp.cnvnator | cnvnator -root 500bp.root -genotype 500 > 1419-10.500bp.genotype
