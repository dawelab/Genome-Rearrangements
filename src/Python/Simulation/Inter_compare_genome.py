
import csv

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

def filter_input(filter_file):  # filter file is the svdetect output file
    filter_inter = []
    with open(filter_file) as raw:
        filter_sum = csv.reader(raw, delimiter="\t")
        for lines in filter_sum:
            filter_inter.append(lines)
    return filter_inter

def rtrans_input(filter_file): 
    filter_inter = []
    with open(filter_file) as raw:
        filter_sum = csv.reader(raw)
        for lines in filter_sum:
            filter_inter.append(lines)
    return filter_inter

def compare_rtrans(lumpy_trans, rtrans):
    common = []
    common_filter = []
    for i in range(len(lumpy_trans)):
        for j in range(len(rtrans)):
            if rtrans[j][0] == lumpy_trans[i][0] and rtrans[j][2] == lumpy_trans[i][10] and (
                        int(rtrans[j][1]) - 20000) < int(lumpy_trans[i][1]) < (int(rtrans[j][1]) + 20000) and (
                        int(rtrans[j][3]) - 500) < int(lumpy_trans[i][11]) < (int(rtrans[j][3]) + 500) and \
                            lumpy_trans[i] not in common and rtrans[j] not in common_filter:
                common.append(lumpy_trans[i])  # common is the output for lumpy format
                common_filter.append(rtrans[j])  # common_filter is the output for rtrans format
            elif rtrans[j][0] == lumpy_trans[i][10] and rtrans[j][2] == lumpy_trans[i][0] and (
                        int(rtrans[j][1]) - 20000) < int(lumpy_trans[i][11]) < (int(rtrans[j][1]) + 20000) and (
                        int(rtrans[j][3]) - 500) < int(lumpy_trans[i][1]) < (int(rtrans[j][3]) + 500) and \
                            lumpy_trans[i] not in common and rtrans[j] not in common_filter:
                common.append(lumpy_trans[i])
                common_filter.append(rtrans[j])
    return common, common_filter

def compare_rtrans2(rtrans, filter_inter):
    common = []
    common_filter = []
    for i in range(len(rtrans)):
        for j in range(len(filter_inter)):
            if filter_inter[j][0] == rtrans[i][0] and filter_inter[j][3] == rtrans[i][2] and (
                        int(filter_inter[j][1]) - 20000) < int(rtrans[i][1]) < (int(filter_inter[j][2]) + 20000) and (
                        int(filter_inter[j][4]) - 500) < int(rtrans[i][3]) < (int(filter_inter[j][5]) + 500) and \
                            rtrans[i] not in common and filter_inter[j] not in common_filter:
                common.append(rtrans[i])  # common is the output for rtrans format
                common_filter.append(filter_inter[j])  # common_filter is the output for filtered format
            elif filter_inter[j][0] == rtrans[i][2] and filter_inter[j][3] == rtrans[i][0] and (
                        int(filter_inter[j][1]) - 500) < int(rtrans[i][3]) < (int(filter_inter[j][2]) + 500) and (
                        int(filter_inter[j][4]) - 20000) < int(rtrans[i][1]) < (int(filter_inter[j][5]) + 20000) and \
                            rtrans[i] not in common and filter_inter[j] not in common_filter:
                common.append(rtrans[i])
                common_filter.append(filter_inter[j])
    return common, common_filter

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
	rtrans = rtrans_input(sim_trans)
	lumpy = lumpy_input(lumpyvcf_file)
	svdetect = filter_input(svdetect_filter_file)
	lumpy_TP = compare_rtrans(lumpy, rtrans)[1]
	lumpy_sensitivity = len(compare_rtrans(lumpy, rtrans)[1])/len(rtrans)
	lumpy_precision = len(compare_rtrans(lumpy, rtrans)[0])/len(lumpy)
	svdetect_TP = compare_rtrans2(rtrans, svdetect)[1]
	svdetect_sensitivity = len(compare_rtrans2(rtrans, svdetect)[0])/len(rtrans)
	svdetect_precision = len(compare_rtrans2(rtrans, svdetect)[1])/len(svdetect)
	lumpy_svdetect_common = compare_filter(lumpy, svdetect)[0]
	len(lumpy_svdetect_common)
	lumpy_svdetect_TP = compare_rtrans(lumpy_svdetect_common, rtrans)[1]
	lumpy_svdetect_sensitivity = len(compare_rtrans(lumpy_svdetect_common, rtrans)[1])/len(rtrans)
	lumpy_svdetect_precision = len(compare_rtrans(lumpy_svdetect_common, rtrans)[0])/len(lumpy_svdetect_common)
