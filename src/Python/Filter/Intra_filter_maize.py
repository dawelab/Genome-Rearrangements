# SVdetect link filter
import csv
import numpy as np
import sys

csv.field_size_limit(1024 * 1024 * 1024)

def filter_input(filter_file):
    filter_sum = []
    with open(filter_file) as raw:
        filter = csv.reader(raw, delimiter="\t")
        for lines in filter:
            filter_sum.append(lines)
    return filter_sum


def lumpy_input(lumpy_file):
    lumpy_sum = []
    lumpy_intra = []
    with open(lumpy_file) as lumpy:
        lum = csv.reader(lumpy, delimiter="\t")
        for lines in lum:
            if lines[0][0] != "#":
                lumpy_sum.append(lines)
    for i in range(len(lumpy_sum) - 1):
        #print(i)
        if lumpy_sum[i][7].split(";")[0] == "SVTYPE=BND" and lumpy_sum[i][2].split("_")[1] == "1" and lumpy_sum[i + 1][2].split("_")[1] == "2" and \
                        lumpy_sum[i][2].split("_")[0] == lumpy_sum[i + 1][2].split("_")[0] and lumpy_sum[i][0] == lumpy_sum[i + 1][0] and abs(int(str(lumpy_sum[i][1]))-int(str(lumpy_sum[i+1][1]))) > 50000 and str(lumpy_sum[i][0]) != "NC_001416.1" and 2 < abs(int(str(lumpy_sum[i][7].split(";")[1].split(":")[1].split(",")[0]))) < 20:
                    lumpy_intra.append(lumpy_sum[i] + lumpy_sum[i + 1])
        elif lumpy_sum[i][7].split(";")[0] != "SVTYPE=BND" and abs(int(str(lumpy_sum[i][7].split(";")[2].split("=")[1]))) > 50000 and 3 < abs(int(str(lumpy_sum[i][7].split(";")[1].split(":")[1].split(",")[0]))) < 20:
            del_1 = [lumpy_sum[i][0], lumpy_sum[i][1], lumpy_sum[i][2]+"_1", lumpy_sum[i][3], lumpy_sum[i][4], lumpy_sum[i][5], lumpy_sum[i][6], lumpy_sum[i][7], lumpy_sum[i][8],lumpy_sum[i][9]]
            del_2 = [lumpy_sum[i][0], lumpy_sum[i][7].split(";")[3].split("=")[1], lumpy_sum[i][2]+"_2", lumpy_sum[i][3], lumpy_sum[i][4], lumpy_sum[i][5], lumpy_sum[i][6], lumpy_sum[i][7], lumpy_sum[i][8],lumpy_sum[i][9]]
            lumpy_intra.append(del_1 + del_2)      
    return lumpy_intra

def svdetect_lumpy_intracompare(lumpy_intra, filter_intra):
    common_lumpy = []
    for j in range(len(filter_intra)):
        for i in range(len(lumpy_intra)):
            if lumpy_intra[i][0] == filter_intra[j][0] and (
                        int(filter_intra[j][1]) - 500) < int(lumpy_intra[i][1]) < (
                        int(filter_intra[j][2]) + 500) and (
                        int(filter_intra[j][4]) - 500) < int(lumpy_intra[i][11]) < (int(filter_intra[j][5]) + 500) and \
                            lumpy_intra[i] not in common_lumpy:
                common_lumpy.append(lumpy_intra[i])
            elif lumpy_intra[i][0] == filter_intra[j][0] and (
                        int(filter_intra[j][4]) - 500) < int(lumpy_intra[i][1]) < (
                        int(filter_intra[j][5]) + 500) and (
                        int(filter_intra[j][1]) - 500) < int(lumpy_intra[i][11]) < (int(filter_intra[j][2]) + 500) and \
                            lumpy_intra[i] not in common_lumpy:
                common_lumpy.append(lumpy_intra[i])
    return common_lumpy
    
def event_filter_wt(common_event, common_wt):
    uniq = []
    common = []
    for j in range(len(common_wt)):
        for i in range(len(common_event)):
            if common_event[i][0] == common_wt[j][0] and (
                        int(common_wt[j][1]) - 50000) < int(common_event[i][1]) < (
                        int(common_wt[j][1]) + 50000) and (
                        int(common_wt[j][11]) - 50000) < int(common_event[i][11]) < (int(common_wt[j][11]) + 50000) and \
                            common_event[i] not in common:
                common.append(common_event[i])
            elif common_event[i][0] == common_wt[j][0] and (
                        int(common_wt[j][11]) - 50000) < int(common_event[i][1]) < (
                        int(common_wt[j][11]) + 50000) and (
                        int(common_wt[j][1]) - 50000) < int(common_event[i][11]) < (int(common_wt[j][1]) + 50000) and \
                            common_event[i] not in common:
                common.append(common_event[i])
    for x in common_event:
    	if x not in common and x not in uniq:
    		uniq.append(x)
    return uniq
    
def output(common, file):
    with open(file, "w") as circos:
        for i in range(len(common)):
            if 2 < abs(int(str(common[i][7].split(";")[1].split(":")[1].split(",")[0]))) < 10:
                if common[i][7].split(";")[0] == "SVTYPE=BND": 
                    circos.write(str(int(i + 1)) + "\t" + "zm" + common[i][0] + "\t" + common[i][1] + "\t" + str(int(common[i][1]) + 1) + "\n" + str(int(i + 1)) + "\t" + "zm" + common[i][10] + "\t" +common[i][11] + "\t" + str(int(common[i][11]) + 1) + "\n")
                else: 
                    circos.write(str(int(i + 1)) + "\t" + "zm" + common[i][0] + "\t" + common[i][1] + "\t" + str(int(common[i][1]) + 1) + "\t" + str(common[i][7].split(";")[0].split("=")[1]) + "\n" + str(int(i + 1)) + "\t" + "zm" + common[i][10] + "\t" +common[i][11] + "\t" + str(int(common[i][11]) + 1) + "\t" + str(common[i][7].split(";")[0].split("=")[1]) + "\n")

def main():
    lumpy_file = sys.argv[1]
    filter_file = sys.argv[2]
    wt_lumpy_file = sys.argv[3]
    wt_filter_file = sys.argv[4]
    circos_file= sys.argv[5]
    filter_intra = filter_input(filter_file)
    lumpy_intra = lumpy_input(lumpy_file)
    common = svdetect_lumpy_intracompare(lumpy_intra, filter_intra)
    #print(common[5])
    #wt_filter_intra = filter_input(wt_filter_file)
    wt_lumpy_intra = lumpy_input(wt_lumpy_file)
    #wt_common = svdetect_lumpy_intracompare(wt_lumpy_intra, wt_filter_intra)
    uniq = event_filter_wt(common, wt_lumpy_intra)
    output(uniq, circos_file)
if __name__ == "__main__":
    main()