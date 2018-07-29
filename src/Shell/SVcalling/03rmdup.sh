#PBS -S /bin/bash
#PBS -q batch
#PBS -o /lustre1/jl03308/20xIllumina-control/analysis/jobs/
#PBS -e /lustre1/jl03308/20xIllumina-control/analysis/jobs/
#PBS -l nodes=1:ppn=12:HIGHMEM
#PBS -l mem=60gb
#PBS -l walltime=80:00:00
#PBS -M jl03308@uga.edu
#PBS -m ae
#PBS -N rmdup_1419-10
module load samtools/1.3.1
module load picard/2.4.1
module load bedtools/2.26.0
cd /lustre1/jl03308/20xIllumina-control/analysis/10A/1419-10/
samtools merge -@ 12 -O BAM 1419-10.bam 1419-10*.bam
samtools sort -o 1419-10.sorted.bam -T 1419-10 -@ 12 1419-10.bam
java -jar /usr/local/apps/picard/2.4.1/picard.jar MarkDuplicates I=/lustre1/jl03308/20xIllumina-control/analysis/10A/1419-10/1419-10.sorted.bam O=/lustre1/jl03308/20xIllumina-control/analysis/10A/1419-10/1419-10.rmdup.bam M=/lustre1/jl03308/20xIllumina-control/analysis/10A/1419-10/1419-10.marked_dup_metrics.txt REMOVE_DUPLICATES= true
samtools sort -o /lustre1/jl03308/20xIllumina-control/analysis/10A/1419-10/1419-10.rmdup.sorted.bam -T 1419-10 -@ 12 /lustre1/jl03308/20xIllumina-control/analysis/10A/1419-10/1419-10.rmdup.bam
samtools index /lustre1/jl03308/20xIllumina-control/analysis/10A/1419-10/1419-10.rmdup.sorted.bam
samtools view -@ 12 -bhq 20 /lustre1/jl03308/20xIllumina-control/analysis/10A/1419-10/1419-10.rmdup.sorted.bam -o /lustre1/jl03308/20xIllumina-control/analysis/10A/1419-10/1419-10.rmdup.q20.bam
samtools sort -o /lustre1/jl03308/20xIllumina-control/analysis/10A/1419-10/1419-10.rmdup.q20.sorted.bam -T 1419-10_q20 -@ 12 /lustre1/jl03308/20xIllumina-control/analysis/10A/1419-10/1419-10.rmdup.q20.bam
samtools index /lustre1/jl03308/20xIllumina-control/analysis/10A/1419-10/1419-10.rmdup.q20.sorted.bam
mkdir /lustre1/jl03308/20xIllumina-control/analysis/10A/1419-10/svdetect
/usr/local/apps/svdetect/latest/scripts/BAM_preprocessingPairs.pl -t 1 -p 1 /lustre1/jl03308/20xIllumina-control/analysis/10A/1419-10/1419-10.rmdup.q20.sorted.bam -o /lustre1/jl03308/20xIllumina-control/analysis/10A/1419-10/svdetect

mkdir /lustre1/jl03308/20xIllumina-control/analysis/10A/1419-10/lumpy
samtools view -b -F 1294 /lustre1/jl03308/20xIllumina-control/analysis/10A/1419-10/1419-10.rmdup.sorted.bam > /lustre1/jl03308/20xIllumina-control/analysis/10A/1419-10/lumpy/1419-10.discordants.unsorted.bam
samtools view -h /lustre1/jl03308/20xIllumina-control/analysis/10A/1419-10/1419-10.rmdup.sorted.bam | /usr/local/apps/lumpy-sv/0.2.13/scripts/extractSplitReads_BwaMem -i stdin > /lustre1/jl03308/20xIllumina-control/analysis/10A/1419-10/lumpy/1419-10.splitters.unsorted.bam
samtools sort -o /lustre1/jl03308/20xIllumina-control/analysis/10A/1419-10/lumpy/1419-10.discordants.bam -T 1419-10  -@ 12 /lustre1/jl03308/20xIllumina-control/analysis/10A/1419-10/lumpy/1419-10.discordants.unsorted.bam
samtools sort -o /lustre1/jl03308/20xIllumina-control/analysis/10A/1419-10/lumpy/1419-10.splitters.bam -T 1419-10  -@ 12 /lustre1/jl03308/20xIllumina-control/analysis/10A/1419-10/lumpy/1419-10.splitters.unsorted.bam
samtools view /lustre1/jl03308/20xIllumina-control/analysis/10A/1419-10/1419-10.rmdup.sorted.bam | tail -n +10000000 |/usr/local/apps/lumpy-sv/0.2.13/scripts/pairend_distro.py -r 75 -X 4 -N 10000000 -o /lustre1/jl03308/20xIllumina-control/analysis/10A/1419-10/lumpy/1419-10.histo

