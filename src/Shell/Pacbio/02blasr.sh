#PBS -S /bin/bash
#PBS -N blasr-1
#PBS -q batch
#PBS -l nodes=1:ppn=24:AMD
#PBS -l walltime=50:00:00
#PBS -l mem=100gb
#PBS -M jl03308@uga.edu
#PBS -m ae


module load smrtlink/5.1.0.26412

cd /lustre1/jl03308/pacbio-RiceE29/data/1341-r1-1310-r9/r54193_20180224_141936/1_A01
blasr m54193_180224_142901-2.ccs.bam /lustre1/jl03308/reference/source/lambda_genome/NCBI_lambda_genome.fa --bam --clipping soft --out m54193_180224_142901-2_ccs_lambda.aligned.bam --unaligned m54193_180224_142901-2_ccs_lambda.unaligned.bam --nproc 22 
blasr m54193_180224_142901-2.ccs.bam /lustre1/jl03308/reference/source/rice_genome_v1/Oryza_sativa.IRGSP-1.0.dna.toplevel.fa --bam --clipping soft --out m54193_180224_142901-2_ccs_rice.aligned.bam --unaligned m54193_180224_142901-2_ccs_rice.unaligned.bam --nproc 22 
