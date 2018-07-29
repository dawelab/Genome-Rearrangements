#PBS -S /bin/bash
#PBS -q batch
#PBS -N RiceE24-sapdes-blast
#PBS -o /lustre1/jl03308/20xIllumina-lambda/analysis/jobs
#PBS -e /lustre1/jl03308/20xIllumina-lambda/analysis/jobs
#PBS -l nodes=1:ppn=12:HIGHMEM
#PBS -l mem=60gb
#PBS -l walltime=100:00:00
#PBS -M jl03308@uga.edu
#PBS -m ae

module load ncbiblast/2.2.26
module load samtools/1.3.1
module load bedtools/2.24.0
module load spades/3.10.0
module load anaconda

cd /lustre1/jl03308/20xIllumina-lambda/analysis/RiceE24/supplement/spades
#Extract the IDs of reads that are mapped to lambda and plasmid 
samtools view -H /lustre1/jl03308/20xIllumina-lambda/analysis/RiceE24/supplement/RiceE24.rmdup.bam > discordants-header.sam
samtools view /lustre1/jl03308/20xIllumina-lambda/analysis/RiceE24/supplement/RiceE24.rmdup.bam |awk '{if($3~"NC"||$7~"NC"){print$0}}'| cat discordants-header.sam - | samtools view -Sb - > NC_001416.bam
samtools view /lustre1/jl03308/20xIllumina-lambda/analysis/RiceE24/supplement/RiceE24.rmdup.bam |awk '{if($3~"pPvUbi2H"||$7~"pPvUbi2H"){print$0}}'| cat discordants-header.sam - | samtools view -Sb - > pPvUbi2H.bam
samtools merge NC_001416-pPvUbi2H.bam NC_001416.bam pPvUbi2H.bam
samtools view NC_001416-pPvUbi2H.bam | cut -f1 > IDs.txt
#Extract the reads in the original rmdup bam file and output them in the sam file 
samtools view /lustre1/jl03308/20xIllumina-lambda/analysis/RiceE24/supplement/RiceE24.rmdup.sorted.bam | LC_ALL=C grep -w -F -f IDs.txt > NC_001416-pPvUbi2H.sam
#Convert sam file format to bam file format 
cat discordants-header.sam NC_001416-pPvUbi2H.sam | samtools view -Sb - > NC_001416-pPvUbi2H.pe.bam
#Sort bam file 
samtools sort -n -o NC_001416-pPvUbi2H.pe.sortedbyname.bam -T NC_001416-pPvUbi2H NC_001416-pPvUbi2H.pe.bam
#Convert bam format to fastq format 
bedtools bamtofastq -i NC_001416-pPvUbi2H.pe.sortedbyname.bam -fq NC_001416-pPvUbi2H.end1.fq -fq2 NC_001416-pPvUbi2H.end2.fq
#Assemble extracted paired-end reads 
/usr/local/apps/spades/3.10.0/bin/spades.py --pe1-1 NC_001416-pPvUbi2H.end1.fq --pe1-2 NC_001416-pPvUbi2H.end2.fq -o NC_001416-pPvUbi2H-pe


cd /lustre1/jl03308/20xIllumina-lambda/analysis/RiceE24/supplement/spades/blast
# /usr/local/apps/ncbiblast/latest/bin/formatdb -p F -i /lustre1/jl03308/reference/source/lambda_genome/NCBI_lambda_genome.fa 
# /usr/local/apps/ncbiblast/latest/bin/formatdb -p F -i /lustre1/jl03308/reference/source/rice_genome_v1/Oryza_sativa.IRGSP-1.0.dna.toplevel.fa
# /usr/local/apps/ncbiblast/latest/bin/formatdb -p F -i /lustre1/jl03308/reference/source/pPvUbi2H_seq/pPvUbi2H.fa
# Against lambda sequence 
/usr/local/apps/ncbiblast/latest/bin/blastall -p blastn -d /lustre1/jl03308/reference/source/lambda_genome/NCBI_lambda_genome.fa -i /lustre1/jl03308/20xIllumina-lambda/analysis/RiceE24/supplement/spades/NC_001416-pPvUbi2H-pe/scaffolds.fasta -m 9 -b 48000 > RiceE24pe-lambda-junction.txt
# Against rice genome sequence 
/usr/local/apps/ncbiblast/latest/bin/blastall -p blastn -d /lustre1/jl03308/reference/source/rice_genome_v1/Oryza_sativa.IRGSP-1.0.dna.toplevel.fa -i /lustre1/jl03308/20xIllumina-lambda/analysis/RiceE24/supplement/spades/NC_001416-pPvUbi2H-pe/scaffolds.fasta -m 9 -b 5000 > RiceE24pe-rice-junction.txt
# Against plasmid 
/usr/local/apps/ncbiblast/latest/bin/blastall -p blastn -d /lustre1/jl03308/reference/source/pPvUbi2H_seq/pPvUbi2H.fa -i /lustre1/jl03308/20xIllumina-lambda/analysis/RiceE24/supplement/spades/NC_001416-pPvUbi2H-pe/scaffolds.fasta -m 9 -b 5000 > RiceE24pe-plasmid-junction.txt
