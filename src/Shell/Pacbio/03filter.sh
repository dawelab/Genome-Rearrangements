#PBS -S /bin/bash
#PBS -q batch
#PBS -o /lustre1/jl03308/20xIllumina-lambda/analysis/jobs/
#PBS -e /lustre1/jl03308/20xIllumina-lambda/analysis/jobs/
#PBS -l nodes=1:ppn=1:HIGHMEM
#PBS -l mem=20gb
#PBS -l walltime=20:00:00
#PBS -M jl03308@uga.edu
#PBS -m ae
#PBS -N filter
module load anaconda3/5.0.0
cd /lustre1/jl03308/pacbio-RiceE29/data/1341-r1-1310-r9/r54193_20180224_141936/1_A01
python /lustre1/jl03308/20xIllumina-control/scr/E29-filter2.py lambda_aligned.txt lambda_aligned_filtered.txt
cat /lustre1/jl03308/pacbio-RiceE29/data/1310-1341/r54193_20180305_215220/2_B01/lambda_aligned_filtered.txt /lustre1/jl03308/pacbio-RiceE29/data/1310-r11-1341-r3/r54193_20180307_143208/2_B01/lambda_aligned_filtered.txt /lustre1/jl03308/pacbio-RiceE29/data/1341-r1-1310-r9/r54193_20180224_141936/1_A01/lambda_aligned_filtered.txt > /lustre1/jl03308/pacbio-RiceE29/data/lambda_aligned_filtered_sum.txt