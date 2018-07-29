#PBS -S /bin/bash
#PBS -q batch
#PBS -o /lustre1/jl03308/20xIllumina-control/analysis/jobs/
#PBS -e /lustre1/jl03308/20xIllumina-control/analysis/jobs/
#PBS -l nodes=1:ppn=2:HIGHMEM
#PBS -l mem=20gb
#PBS -l walltime=5:00:00
#PBS -M jl03308@uga.edu
#PBS -m ae
#PBS -N sv_1419-10
module load samtools/1.3.1
module load bwa/0.7.15
module load svdetect/0.7
module load lumpy-sv/0.2.13
cd /lustre1/jl03308/20xIllumina-control/analysis/10A/1419-10/lumpy
lumpy -mw 2 -tt 0 -pe id:1419-10,bam_file:1419-10.discordants.bam,histo_file:1419-10.histo,mean:262.939194532,stdev:91.7518785109,read_length:75,min_non_overlap:75,discordant_z:5,back_distance:10,weight:1,min_mapping_threshold:20 -sr id:1419-10,bam_file:1419-10.splitters.bam,back_distance:10,weight:1,min_mapping_threshold:20 > 1419-10-2reads.vcf
cd /lustre1/jl03308/20xIllumina-control/analysis/10A/1419-10/svdetect
SVDetect linking -conf /lustre1/jl03308/20xIllumina-control/analysis/cfg/1419-10-intra.cfg
SVDetect filtering -conf /lustre1/jl03308/20xIllumina-control/analysis/cfg/1419-10-intra.cfg
SVDetect links2SV -conf /lustre1/jl03308/20xIllumina-control/analysis/cfg/1419-10-intra.cfg
SVDetect links2circos -conf /lustre1/jl03308/20xIllumina-control/analysis/cfg/1419-10-intra.cfg
SVDetect links2bed -conf /lustre1/jl03308/20xIllumina-control/analysis/cfg/1419-10-intra.cfg
SVDetect linking -conf /lustre1/jl03308/20xIllumina-control/analysis/cfg/1419-10-inter.cfg
SVDetect filtering -conf /lustre1/jl03308/20xIllumina-control/analysis/cfg/1419-10-inter.cfg
SVDetect links2SV -conf /lustre1/jl03308/20xIllumina-control/analysis/cfg/1419-10-inter.cfg
SVDetect links2circos -conf /lustre1/jl03308/20xIllumina-control/analysis/cfg/1419-10-inter.cfg
SVDetect links2bed -conf /lustre1/jl03308/20xIllumina-control/analysis/cfg/1419-10-inter.cfg
