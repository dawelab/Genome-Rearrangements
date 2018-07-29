
#PBS -S /bin/bash
#PBS -N smrtlink-ccs-1
#PBS -q batch
#PBS -l nodes=1:ppn=12:AMD
#PBS -l walltime=500:00:00
#PBS -l mem=100gb
#PBS -M jl03308@uga.edu
#PBS -m ae


module load smrtlink/5.1.0.26412

cd /lustre1/jl03308/pacbio-RiceE29/data/1341-r1-1310-r9/r54193_20180224_141936/1_A01

ccs --numThreads=12 --polish --minLength=50 --maxLength=30000 \
--minPasses=0 --minPredictedAccuracy=0.8 --minZScore=-3.4 \
--maxDropFraction=0.34 --minPredictedAccuracy=0.8 --minSnr=3.75 \
--reportFile=m54193_180224_142901-2.ccs.report \
m54193_180224_142901.subreads.bam m54193_180224_142901-2.ccs.bam
