#PBS -S /bin/bash
#PBS -q batch
#PBS -N ds20-riceE10
#PBS -o /lustre1/jl03308/20xIllumina-lambda/analysis/jobs
#PBS -e /lustre1/jl03308/20xIllumina-lambda/analysis/jobs
#PBS -l nodes=1:ppn=12:HIGHMEM
#PBS -l mem=50gb
#PBS -l walltime=300:00:00
#PBS -M jl03308@uga.edu
#PBS -m ae
module load gatk/3.6
module load picard/2.4.1
module load samtools/1.3.1 

cd /lustre1/jl03308/20xIllumina-lambda/analysis/RiceE10/supplement/
# create sequence dictionary for reference 
# java -jar /usr/local/apps/picard/2.4.1/picard.jar CreateSequenceDictionary R=/lustre1/jl03308/reference/concatenated_genome/rice_pPvUbi2H_lambdax1.fa O=/lustre1/jl03308/reference/concatenated_genome/rice_pPvUbi2H_lambdax1.dict
# add read groups to input bam file before running downsampling 
## as SAM file doesn't have any read groups defined in the header.  The GATK no longer supports SAM files without read groups
# java -jar /usr/local/apps/picard/2.4.1/picard.jar AddOrReplaceReadGroups I=RiceE10.rmdup.bam O=RiceE10.rmdup.arg.bam RGID=1 RGLB=lib1 RGPL=illumina RGPU=unit2 RGSM=20
# index the bam file with read groups added 
##Cannot process the provided BAM/CRAM file(s) because they were not indexed.  The GATK does offer limited processing of unindexed BAM/CRAMs in --unsafe mode, but this feature is unsupported -- use it at your own risk!
# samtools index RiceE10.rmdup.arg.bam
# downsample with gatk 
java -jar /usr/local/apps/gatk/latest/GenomeAnalysisTK.jar -T PrintReads -R /lustre1/jl03308/reference/concatenated_genome/rice_pPvUbi2H_lambdax1.fa -I RiceE10.rmdup.arg.bam -o RiceE10.rmdup.ds20.bam -dt ALL_READS -dfrac 0.8  

