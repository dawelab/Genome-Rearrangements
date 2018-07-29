#PBS -S /bin/bash
#PBS -q batch
#PBS -N mummer
#PBS -o /lustre1/jl03308/20xIllumina-control/analysis/jobs/
#PBS -e /lustre1/jl03308/20xIllumina-control/analysis/jobs/
#PBS -l nodes=1:ppn=4:AMD
#PBS -l mem=30gb
#PBS -l walltime=20:00:00
#PBS -M jl03308@uga.edu
#PBS -m ae
module load mummer/3.23

cd /lustre1/jl03308/pacbio-RiceE29/data 
nucmer --prefix=pacbio-sum /lustre1/jl03308/reference/source/lambda_genome/NCBI_lambda_genome.fa RiceE29_prinseq_good_7Dmd.fastq
show-coords -rcl pacbio-sum.delta > pacbio-sum.coords
delta-filter -1 pacbio-sum.delta > pacbio-sum-filtered.delta
#show-aligns pacbio-sum.delta refname qryname > ref_qry.aligns
#delta-filter -1 pacbio-sum.delta > pacbio-filtered-filtered.delta
mummerplot --size large -fat --color -f --png pacbio-sum-filtered.delta -p pacbio-sum-filtered

#show-coords -rcl chr9-filtered.delta > chr9-filtered.coords
#show-tiling chr9-filtered.delta > chr9-filtered.tiling
#mapview -n 1 -f pdf -p mapview chr9-filtered.coords
# nucmer --prefix=pacbio /lustre1/jl03308/reference/source/lambda_genome/NCBI_lambda_genome.fa lambda-pacbio.fq
# delta-filter -1 pacbio.delta > pacbio-filtered.delta
# mummerplot --size large -fat --color -f --png pacbio-filtered.delta -p pacbio-filtered

cd /lustre1/jl03308/pacbio-RiceE29/data 
nucmer --prefix=chr1 Modified_chr1.fa Oryza_sativa.IRGSP-1.0.dna.chromosome.1.fa
show-coords -rcl chr1.delta > chr1.coords
delta-filter -1 chr1.delta > chr1-filtered.delta
#show-aligns chr1.delta chr1modified qryname > ref_qry.aligns
#delta-filter -1 pacbio-sum.delta > pacbio-filtered-filtered.delta
mummerplot --size large -fat --color -f --png chr1-filtered.delta -p chr1-filtered

nucmer --prefix=lambda-1 Modified_chr1-1.fa NCBI_lambda_genome.fa
show-coords -rcl lambda-1.delta > lambda-1.coords
delta-filter -1 lambda-1.delta > lambda-filtered-1.delta
#show-aligns chr1.delta chr1modified qryname > ref_qry.aligns
#delta-filter -1 pacbio-sum.delta > pacbio-filtered-filtered.delta
mummerplot --size large -fat --color -f --png lambda-filtered-1.delta -p lambda-filtered-1

