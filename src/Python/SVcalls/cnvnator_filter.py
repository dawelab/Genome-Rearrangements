#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  8 11:23:04 2018

@author: JianingLiu
"""

import csv
import sys
import math
import numpy as np


def load_genotype(file):
    cnv = []
    with open(file) as cnvator:
        c = csv.reader(cnvator, delimiter=" ")
        for lines in c:
            cnv.append(lines)
    return cnv


def load(file):
    cnv = []
    with open(file) as cnvator:
        c = csv.reader(cnvator, delimiter="\t")
        for lines in c:
            cnv.append(lines)
    return cnv


# Compare event with wt
def compare(wt_genotype, event_genotype, event):
    copy_number = []
    num_wt_genotype = []
    for item in wt_genotype:
        temp = [item[1].split(":")[0], item[1].split(":")[1].split("-")[0], item[1].split(":")[1].split("-")[1],
                item[3]]
        tempnum = [float(xx) for xx in temp]
        num_wt_genotype.append(tempnum)

    num_event_genotype = []
    for i, item in enumerate(event_genotype):
        temp = [item[1].split(":")[0], item[1].split(":")[1].split("-")[0], item[1].split(":")[1].split("-")[1],
                item[3]]

        tempnum = [float(xx) for xx in temp] + [1.0 if event[i][0] == "deletion" else 0.0, float(event[i][4]),
                                                float(event[i][8])]
        num_event_genotype.append(tempnum)

    arr_wt_genotype = np.array(num_wt_genotype)
    arr_event_genotype = np.array(num_event_genotype)

    for xx in sorted(np.unique(arr_wt_genotype[..., 0])):
        sub_wt = arr_wt_genotype[np.where(arr_wt_genotype[..., 0] == xx)]
        sub_event = arr_event_genotype[np.where(arr_event_genotype[..., 0] == xx)]
        sub_event = sub_event[np.where(sub_event[..., -2] < 0.01)]
        sub_event = sub_event[np.where(sub_event[..., -1] >= 0)]
        unique_found = [False for xxx in range(len(sub_event[..., 0]))]
        for j, evt in enumerate(sub_event):
            for i, wt in enumerate(sub_wt):
                cutoff = max([abs(wt[1] - wt[2]),abs(evt[1] - evt[2])])/3
                if abs(wt[1] - evt[1]) <= cutoff and abs(evt[2] - wt[2]) <= cutoff: 
                    unique_found[j] = True
                    if abs(evt[3] - wt[3]) > 0.5:
                        tmp = ["deletion" if evt[-3] == 1.0 else "duplication", int(xx), int(evt[1]), int(evt[2]),
                               int(evt[2] - evt[1]), evt[3], 'n']
                        if 10 > wt[3] > 1.5 and 0.8 < evt[3] < 1.2 and evt[-3] == 1.0 and tmp not in copy_number:
                            copy_number.append(tmp)
                        if 0.5 < wt[3] < 2.5 and 10 > evt[3] > 2.8 and evt[-3] == 0.0 and tmp not in copy_number:
                            copy_number.append(tmp)
#                elif len(np.intersect1d(np.arange(evt[1], evt[2]), np.arange(wt[1], wt[2]))) > abs(evt[1]-evt[2])/2: 
#                    unique_found[j] = True
#                    if abs(evt[3] - wt[3]) > 0.5:
#                        tmp = ["deletion" if evt[-3] == 1.0 else "duplication", int(xx), int(evt[1]), int(evt[2]),
#                               int(evt[2] - evt[1]), evt[3], 'n']
#                        if 10 > wt[3] > 1.5 and 0.8 < evt[3] < 1.2 and evt[-3] == 1.0 and tmp not in copy_number:
#                            copy_number.append(tmp)
#                        if 0.5 < wt[3] < 2.5 and 10 > evt[3] > 2.8 and evt[-3] == 0.0 and tmp not in copy_number:
#                            copy_number.append(tmp)

        for j, evt in enumerate(sub_event):
            if not unique_found[j]:
                if 0.8 < evt[3] < 1.2 and evt[-3] == 1.0 and abs(evt[1] - evt[2]) > 100000:
                    copy_number.append(
                        ["deletion" if evt[-3] == 1.0 else "duplication", int(xx), int(evt[1]), int(evt[2]),
                         int(evt[2] - evt[1]), evt[3], 'y'])
                if 2.8 < evt[3] < 10 and evt[-3] == 0.0 and abs(evt[1] - evt[2]) > 100000:
                    copy_number.append(
                        ["deletion" if evt[-3] == 1.0 else "duplication", int(xx), int(evt[1]), int(evt[2]),
                         int(evt[2] - evt[1]), evt[3], 'y'])
    return copy_number
    #
    #
    # for i in range(4200, 4300):
    #     x += 1
    #     for j in range(4200, 4300):
    #         y += 1
    #         chromosome_event = event_genotype[j][1].split(":")[0]
    #         chromosome_wt = wt_genotype[i][1].split(":")[0]
    #         start_coordinates_event = int(event_genotype[j][1].split(":")[1].split("-")[0])
    #         start_coordinates_wt = int(wt_genotype[i][1].split(":")[1].split("-")[0])
    #         end_coordinates_event = int(event_genotype[j][1].split(":")[1].split("-")[1])
    #         end_coordinates_wt = int(wt_genotype[i][1].split(":")[1].split("-")[1])
    #         genovalue_wt = float(wt_genotype[i][3])
    #         genovalue_event = float(event_genotype[j][3])
    #
    #         cnv_call_event = event[j][0]
    #         evalue_1 = float(event[j][4])
    #         q0 = int(event[j][8])
    #         if chromosome_wt == chromosome_event and evalue_1 < 0.01 and q0 >= 0:
    #             # Filter 1 - Quality ( evalue_1 < 0.01 and q0 != -1)
    #
    #                 z += 1
    #
    #                 if abs(genovalue_event - genovalue_wt) <= 0.5:
    #                     common.append(event[j])
    #                 else:
    #                     aa += 1
    #                     if 10 > genovalue_wt > 1.5 and 0.5 < genovalue_event < 1.5 and cnv_call_event == "deletion" and \
    #                                     event[j] not in cnv:
    #                         cnv.append(event[j])
    #                         # aa += 1
    #                         copy_number.append(
    #                             [cnv_call_event, chromosome_event, start_coordinates_event, end_coordinates_event,
    #                              end_coordinates_event - start_coordinates_event, genovalue_event,
    #                              abs(genovalue_event - genovalue_wt), start_coordinates_wt, start_coordinates_event,
    #                              end_coordinates_wt, end_coordinates_event])
    #                     if 0.5 < genovalue_wt < 2.5 and 10 > genovalue_event > 2.5 and cnv_call_event == "duplication" and \
    #                                     event[
    #                                         j] not in cnv:
    #                         cnv.append(event[j])
    #                         # aa += 1
    #                         copy_number.append(
    #                             [cnv_call_event, chromosome_event, start_coordinates_event, end_coordinates_event,
    #                              end_coordinates_event - start_coordinates_event, genovalue_event,
    #                              abs(genovalue_event - genovalue_wt), start_coordinates_wt, start_coordinates_event,
    #                              end_coordinates_wt, end_coordinates_event])
    #                     elif event[j] not in common:
    #                         common.append(event[j])
    #             else:
    #                 if 0.5 < genovalue_event < 1.5 and cnv_call_event == "deletion" and event[j] not in cnv:
    #                     cnv.append(event[j])
    #                     copy_number.append(
    #                         [cnv_call_event, chromosome_event, start_coordinates_event, end_coordinates_event,
    #                          end_coordinates_event - start_coordinates_event, genovalue_event, 0, j, i,
    #                          end_coordinates_wt, end_coordinates_event])
    #                 if 2.5 < genovalue_event < 10 and cnv_call_event == "duplication" and event[j] not in cnv:
    #                     cnv.append(event[j])
    #                     copy_number.append(
    #                         [cnv_call_event, chromosome_event, start_coordinates_event, end_coordinates_event,
    #                          end_coordinates_event - start_coordinates_event, genovalue_event, 0, j, i,
    #                          end_coordinates_wt, end_coordinates_event])
    #                 elif event[j] not in common:
    #                     common.append(event[j])
    #
    # print(x)
    # print(y)
    # print(z, aa)
    # print(len(common))


# def cov_change(cnv):
#     cnv_final = []
#     for i in range(len(cnv)):
#         chromosome = cnv[i][1].split(":")[0]
#         rice_chromosome = "os" + str(chromosome)
#         start_coordinate = cnv[i][1].split(":")[1].split("-")[0]
#         end_coordinate = cnv[i][1].split(":")[1].split("-")[1]
#         if cnv[i][3] == "inf":
#             cnv_final.append([rice_chromosome, start_coordinate, end_coordinate, 100])
#         else:
#             cnv_final.append([rice_chromosome, start_coordinate, end_coordinate, math.log(float(cnv[i][3])*2,2)])
#     return cnv_final


def main():
    sample = sys.argv[1]
    control_genotype = sys.argv[2]
    sample_genotype = sys.argv[3]
    unique_file = sys.argv[4]
    event = load(sample)
    wt_genotype = load_genotype(control_genotype)
    event_genotype = load_genotype(sample_genotype)
    unique = compare(wt_genotype, event_genotype, event)
    print(len(event_genotype))
    print(len(wt_genotype))
    print(len(unique))
    with open(unique_file, "w", newline='') as norm:
        sw = csv.writer(norm, delimiter='\t')
        for item in unique:
            temp = ["zm" + str(item[1])] + item[2:4] + [item[5]]
            sw.writerow(temp)
            # for x in range(len(unique)):
            #     norm.write("os" + str(unique[x][1]) + "\t" + str(unique[x][2]) + "\t" + str(unique[x][3]) + "\t" + str(
            #         unique[x][5]) + "\t" + str(unique[x][6]) + "\t" + str(unique[x][7]) + "\t" + str(
            #         unique[x][8]) + "\t" + str(unique[x][9]) + "\t" + str(unique[x][10]) + "\n")


if __name__ == "__main__":
    main()
