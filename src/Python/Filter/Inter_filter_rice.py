# Goal of this script: Extract the reads that support inter-chromosomal translocations between lambda and chromosomes as well as between chromosomes around lambda insertion regions

# In lumpy file, extract the inter-chromosomal calls between lambda/plasmid and chormosomes

## load lumpy file, extract inter-chromosomal translcoations between lambda/plasmid and chormosomes
import csv
import sys

csv.field_size_limit(1024 * 1024 * 1024)


def lumpy_input(lumpy_file):
    lumpy_inter_lambda = []
    lumpy_inter_chromosomes = []
    lumpy_sum = []
    with open(lumpy_file) as lumpy:
        lum = csv.reader(lumpy, delimiter="\t")
        for lines in lum:
            if lines[0][0] != "#":
                lumpy_sum.append(lines)
    for i in range(len(lumpy_sum) - 1):
        if lumpy_sum[i][7].split(";")[0] == "SVTYPE=BND" and lumpy_sum[i][0] != lumpy_sum[i + 1][0] and \
                        lumpy_sum[i][2].split("_")[0] == lumpy_sum[i+1][2].split("_")[0]:
            if lumpy_sum[i][0] != "NC_001416.1" and lumpy_sum[i][0] != "pPvUbi2H" and lumpy_sum[i + 1][
                0] != "NC_001416.1" and lumpy_sum[i + 1][0] != "pPvUbi2H":
                lumpy_inter_chromosomes.append(lumpy_sum[i] + lumpy_sum[i + 1])
            else:
                lumpy_inter_lambda.append(lumpy_sum[i] + lumpy_sum[i + 1])
    return lumpy_inter_lambda, lumpy_inter_chromosomes


def unique(sample_inter, wt_inter):
    uniq = []
    common = []
    for i in range(len(sample_inter)):
        for j in range(len(wt_inter)):
            if sample_inter[i][0] == wt_inter[j][0] and sample_inter[i][10] == wt_inter[j][10] and (
                        int(sample_inter[i][1]) - 10000) < int(wt_inter[j][1]) < (int(sample_inter[i][1]) + 10000) and (
                        int(sample_inter[i][11]) - 10000) < int(wt_inter[j][11]) < (int(sample_inter[i][11]) + 10000):
                common.append(sample_inter[i])
            elif sample_inter[i][0] == wt_inter[j][10] and sample_inter[i][10] == wt_inter[j][0] and (
                        int(sample_inter[i][11]) - 10000) < int(wt_inter[j][1]) < (int(sample_inter[i][11]) + 10000) and (
                        int(sample_inter[i][1]) - 10000) < int(wt_inter[j][11]) < (int(sample_inter[i][1]) + 10000):
                common.append(sample_inter[i])
    for x in sample_inter:
        if x not in common:
            uniq.append(x)
    return uniq

def filter_input(filter_file):  # filter file is the svdetect output file
    filter_inter = []
    with open(filter_file) as raw:
        filter_sum = csv.reader(raw, delimiter="\t")
        for lines in filter_sum:
            filter_inter.append(lines)
    return filter_inter

## Compare lumpy file inter-chromosomal output with filtered file from svdetect, extract the reads id from svdetect file
def compare_filter(lumpy_trans, filter_inter):
    common = []
    common_filter = []
    for i in range(len(lumpy_trans)):
        for j in range(len(filter_inter)):
            if filter_inter[j][0] == lumpy_trans[i][0] and filter_inter[j][3] == lumpy_trans[i][10] and (
                        int(filter_inter[j][1]) - 500) < int(lumpy_trans[i][1]) < (int(filter_inter[j][2]) + 500) and (
                        int(filter_inter[j][4]) - 500) < int(lumpy_trans[i][11]) < (int(filter_inter[j][5]) + 500) and \
                            lumpy_trans[i] not in common and filter_inter[j] not in common_filter:
                common.append(lumpy_trans[i])  # common is the output for lumpy format
                common_filter.append(filter_inter[j])  # common_filter is the output for filtered format
            elif filter_inter[j][0] == lumpy_trans[i][10] and filter_inter[j][3] == lumpy_trans[i][0] and (
                        int(filter_inter[j][1]) - 500) < int(lumpy_trans[i][11]) < (int(filter_inter[j][2]) + 500) and (
                        int(filter_inter[j][4]) - 500) < int(lumpy_trans[i][1]) < (int(filter_inter[j][5]) + 500) and \
                            lumpy_trans[i] not in common and filter_inter[j] not in common_filter:
                common.append(lumpy_trans[i])
                common_filter.append(filter_inter[j])
    return common, common_filter

def main():
    sample_lumpy_file = sys.argv[1]
    wt_lumpy_file = sys.argv[2]
    filter_file = sys.argv[3]
    circos_file = sys.argv[4]
    filter_inter = filter_input(filter_file)  # Load in the filtered output of svdetect, which contain inter-chromosomal translocations and read IDs that support specific translocations
    sample_inter_chromosomes = lumpy_input(sample_lumpy_file)[
        1]  # lumpy translocations between chromosomes and chromosomes for sample
    # Read lumpy file for wt
    wt_inter_chromosomes = lumpy_input(wt_lumpy_file)[
        1]  # lumpy translocations between chromosomes and chromosomes for wt
    common = compare_filter(sample_inter_chromosomes, filter_inter)[
        0]  # common links between lumpy and svdetect in total (including lambda-chromosomal and chromosomal-chromosomal), in svdetect filtered file format
    common_sum = unique(common, wt_inter_chromosomes)
    with open(circos_file, "w") as circos:
        for i in range(len(common_sum)):
            if 3 < abs(int(str(common_sum[i][7].split(";")[1].split(":")[1].split(",")[0]))) < 20: 
                circos.write(str(int(i + 1)) + "\t" + "os" + common_sum[i][0] + "\t" + common_sum[i][1] + "\t" + str(int(common_sum[i][1]) + 1) + "\n" + str(int(i + 1)) + "\t" + "os" + common_sum[i][10] + "\t" +common_sum[i][11] + "\t" + str(int(common_sum[i][11]) + 1) + "\n")

    print(len(sample_inter_chromosomes))
    print(len(common))
    print(len(common_sum))
if __name__ == "__main__":
    main()
