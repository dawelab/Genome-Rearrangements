#PBS -S /bin/bash
#PBS -q batch
#PBS -o /lustre1/jl03308/20xIllumina-control/analysis/jobs/
#PBS -e /lustre1/jl03308/20xIllumina-control/analysis/jobs/
#PBS -l nodes=1:ppn=12:HIGHMEM
#PBS -l mem=60gb
#PBS -l walltime=25:00:00
#PBS -M jl03308@uga.edu
#PBS -m ae
#PBS -N chr3
module load art/20160605
module load bwa/0.7.15
module load samtools/1.3.1
module load picard/2.4.1
module load bedtools/2.26.0
cd /lustre1/jl03308/20xIllumina-lambda/analysis/simfinal/Chr3
art_illumina -ss NS50 -i Modified_chr1-1copy-2500_diploid.fa -p -l 75 -f 10 -m 300 -s 80 -na -o PE75_art
module load python/2.7.8
module load cutadapt/1.9.dev1 fastqc/0.11.3 trimgalore/0.4.4
module load java/jdk1.8.0_20 fastqc
/usr/local/apps/trimgalore/4.0/trim_galore --fastqc --gzip --paired PE75_art1.fq PE75_art2.fq
bwa mem -t 12 /lustre1/jl03308/20xIllumina-lambda/analysis/simfinal/unmodified/Chr3_chr1.fa PE75_art1_val_1.fq.gz PE75_art2_val_2.fq.gz  > 1copy_20cov.sam

samtools view -b -o 1copy_20cov.bam 1copy_20cov.sam
samtools sort -o 1copy_20cov.sorted.bam -T 1copy_20cov -@ 12 1copy_20cov.bam
samtools index 1copy_20cov.sorted.bam
java -jar /usr/local/apps/picard/2.4.1/picard.jar MarkDuplicates I=1copy_20cov.sorted.bam O=1copy_20cov.rmdup.bam M=1copy_20cov.marked_dup_metrics.txt REMOVE_DUPLICATES= true
samtools sort -o 1copy_20cov.rmdup.sorted.bam -T 1copy_20cov -@ 12 1copy_20cov.rmdup.bam
samtools index 1copy_20cov.rmdup.sorted.bam
samtools view -@ 12 -bhq 20 1copy_20cov.rmdup.sorted.bam -o 1copy_20cov.rmdup.q20.bam
samtools sort -o 1copy_20cov.rmdup.q20.sorted.bam -T 1copy_q20 -@ 12 1copy_20cov.rmdup.q20.bam
samtools index 1copy_20cov.rmdup.q20.sorted.bam
cd svdetect
/usr/local/apps/svdetect/latest/scripts/BAM_preprocessingPairs.pl -t 1 -p 1 ../1copy_20cov.rmdup.q20.sorted.bam
cd ../
cd lumpy
samtools view -b -F 1294 ../1copy_20cov.rmdup.bam > 1copy_20cov.discordants.unsorted.bam
samtools view -h ../1copy_20cov.rmdup.bam | /usr/local/apps/lumpy-sv/0.2.13/scripts/extractSplitReads_BwaMem -i stdin > 1copy_20cov.splitters.unsorted.bam
samtools sort -o 1copy_20cov.discordants.bam -T 1copy_20cov -@ 12 1copy_20cov.discordants.unsorted.bam
samtools sort -o 1copy_20cov.splitters.bam -T 1copy_20cov -@ 12 1copy_20cov.splitters.unsorted.bam
samtools view ../1copy_20cov.rmdup.sorted.bam | tail -n +500000 |/usr/local/apps/lumpy-sv/0.2.13/scripts/pairend_distro.py -r 75 -X 4 -N 500000 -o 1copy_20cov.histo
# module load samtools/1.3.1
# module load bwa/0.7.15
module load svdetect/0.7
module load lumpy-sv/0.2.13
lumpy -mw 2 -tt 0 -pe id:1copy_20cov,bam_file:1copy_20cov.discordants.bam,histo_file:1copy_20cov.histo,mean:300,stdev:80,read_length:75,min_non_overlap:75,discordant_z:5,back_distance:10,weight:1,min_mapping_threshold:20 -sr id:1copy_20cov,bam_file:1copy_20cov.splitters.bam,back_distance:10,weight:1,min_mapping_threshold:20 > 1copy_20cov-2reads.vcf
cd ../svdetect
SVDetect linking -conf 1copy-intra.cfg
SVDetect filtering -conf 1copy-intra.cfg
SVDetect links2SV -conf 1copy-intra.cfg
SVDetect links2circos -conf 1copy-intra.cfg
SVDetect links2bed -conf 1copy-intra.cfg
SVDetect linking -conf 1copy-inter.cfg
SVDetect filtering -conf 1copy-inter.cfg
SVDetect links2SV -conf 1copy-inter.cfg
SVDetect links2circos -conf 1copy-inter.cfg
SVDetect links2bed -conf 1copy-inter.cfg

