#PBS -S /bin/bash
#PBS -q batch
#PBS -o /lustre1/jl03308/20xIllumina-control/analysis/jobs/
#PBS -e /lustre1/jl03308/20xIllumina-control/analysis/jobs/
#PBS -l nodes=1:ppn=1:HIGHMEM
#PBS -l mem=30gb
#PBS -l walltime=30:00:00
#PBS -M jl03308@uga.edu
#PBS -m ae
#PBS -N covcal_1419-10
module load bedtools/2.26.0
module load samtools/1.3.1
cd /lustre1/jl03308/20xIllumina-control/analysis/10A/1419-10/
mkdir cov_cal
cd cov_cal
# 10kb segment coverage calculation MAPQ20 (copycat)
/lustre1/jl03308/20xIllumina-control/scr/copycat-master/copycat ../1419-10.rmdup.q20.sorted.bam /lustre1/jl03308/20xIllumina-control/analysis/cfg/10A.len 1419-10_q20

# 10kb segment coverage calculation  (copycat)
/lustre1/jl03308/20xIllumina-control/scr/copycat-master/copycat ../1419-10.rmdup.sorted.bam /lustre1/jl03308/20xIllumina-control/analysis/cfg/10A.len 1419-10

# 1bp position coverage calculation MAPQ20 (bedtools)
bedtools genomecov -ibam  ../1419-10.rmdup.q20.sorted.bam -g /lustre1/jl03308/20xIllumina-control/analysis/cfg/10A.len -bga > 1419-10_q20-coverage-bedtools.txt

# 1bp position coverage calculation (bedtools)
bedtools genomecov -ibam  ../1419-10.rmdup.sorted.bam -g /lustre1/jl03308/20xIllumina-control/analysis/cfg/10A.len -bga > 1419-10-coverage-bedtools.txt

cat 1419-10_q20-coverage-bedtools.txt | awk '{if($1~"10A"){print$0}}' > 1419-10_q20-10A.coverage.1bp
cat 1419-10_q20-10A.coverage.1bp | sed 's/^/os/g' > 1419-10_q20-10A-circos.coverage.1bp
cat 1419-10_q20.coverage.10kb.for_IGV.seg | cut -f2,3,4,5 |sed 's/^/os/g' > 1419-10_q20-10kb-circos.txt
cat 1419-10_q20-10A-circos.coverage.1bp 1419-10_q20-10kb-circos.txt > 1419-10_q20_combined.txt
