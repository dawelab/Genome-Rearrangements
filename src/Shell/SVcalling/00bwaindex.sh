#PBS -S /bin/bash
#PBS -q batch
#PBS -N BWAindex
#PBS -o /lustre1/jl03308/lambdaX20/analysis/jobs
#PBS -e /lustre1/jl03308/lambdaX20/analysis/jobs
#PBS -l nodes=1:ppn=4:AMD
#PBS -l mem=30gb
#PBS -l walltime=480:00:00
#PBS -M jl03308@uga.edu
#PBS -m ae

# bwa index concatenated genome 
cd /lustre1/jl03308/lambdaX20/reference/source/
module load bwa/0.7.15
bwa index ./maize_genome_v4/Zea_mays.AGPv4.dna.toplevel.fa
bwa index ./rice_genome_v1/Oryza_sativa.IRGSP-1.0.dna.toplevel.fa