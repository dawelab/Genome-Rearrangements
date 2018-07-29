
import matplotlib.pyplot as plt
import numpy as np
import csv
import sys
import math
import os
import matplotlib.ticker as ticker

def load_genotype(file):
    cnv = []
    with open(file) as cnvator:
        c = csv.reader(cnvator, delimiter="\t")
        for lines in c:
            cnv.append(lines)
    return cnv

def cluster(cnv, start, end):
    cnv_final = []
    cnv = sorted(cnv, key=lambda x : int(x[2]))

    if int(start) < int(cnv[0][1]) - 1:
        cnv_final.append([cnv[0][0], int(start), int(cnv[0][1]), 2])
    if int(end) > int(cnv[-1][2]) + 1:
        cnv_final.append([cnv[-1][0], int(cnv[-1][2]) + 1, int(end) + 1, 2])
    for i in range(len(cnv)):
        cnv_final.append([cnv[i][0], int(cnv[i][1]), int(cnv[i][2]) + 1, round(float(cnv[i][3]))])
    for i in range(len(cnv)-1):
        if int(cnv[i+1][1]) - int(cnv[i][2]) != 1:
            cnv_final.append([cnv[i][0], int(cnv[i][2])+1, int(cnv[i+1][1]), 2])
    return cnv_final
def load_inter(file, chromosome):
    coord = []
    with open(file) as trans:
        t = csv.reader(trans, delimiter="\t")
        for lines in t:
            if lines[1] == str(chromosome):
                coord.append(lines[2])
    return coord

def load_intra(file):
    dups = []
    dels = []
    invs = []
    intras = []
    with open(file) as trans:
        t = csv.reader(trans, delimiter="\t")
        for lines in t:
            if len(lines) == 4:
                intras.append(lines[2])
            elif len(lines) == 5:
                if lines[4] == "DUP":
                    dups.append(lines[2])
                elif lines[4] == "INV":
                    invs.append(lines[2])
                elif lines[4] == "DEL":
                    dels.append(lines[2])
              
    return dups,dels,invs,intras


def plot():
    inter_chrom = load_inter('1419-12-inter.txt', "os1")
    inter_lambda = load_inter('1419-12-inter-lambda.txt', "os1")
    intras = load_intra('1419-12-intra.txt')[3]
    dups = load_intra('1419-12-intra.txt')[0]
    dels = load_intra('1419-12-intra.txt')[1]
    invs = load_intra('1419-12-intra.txt')[2]
    cnv = load_genotype('1419-12-chr1.txt')
    cnv_final = cluster(cnv, 18000000, 22000000)
    fig = plt.figure(figsize=(12,3),facecolor='whitesmoke')

    ax = fig.add_subplot(1, 1, 1)
    ax.set_yticks(range(5))
    ax.set_ylim(0, 4)
    ax.set_xlim(18000000,22000000)
    #plt.axvline(x=32,ymin=0.25, ymax=0.5)
    plt.yticks(range(5))
    for i in range(len(intras)):
        plt.axvline(x=int(intras[i]),ymin=0.0, ymax=0.12,linewidth=1,color='blue',zorder=10)
    for i in range(len(dups)):
        plt.axvline(x=int(dups[i]),ymin=0.0, ymax=0.12,linewidth=1,color='orangered',zorder=10)
    for i in range(len(dels)):
        plt.axvline(x=int(dels[i]),ymin=0.0, ymax=0.12,linewidth=1,color='grey',zorder=10)
    for i in range(len(invs)):
        plt.axvline(x=int(invs[i]),ymin=0.0, ymax=0.12,linewidth=1,color='darkgreen',zorder=10)
    for i in range(len(inter_chrom)):
        plt.axvline(x=int(inter_chrom[i]),ymin=0.0, ymax=0.24,linewidth=1,color='dimgrey',zorder=5)
    for i in range(len(inter_lambda)):
        plt.axvline(x=int(inter_lambda[i]),ymin=0.0, ymax=0.24,linewidth=1,color='violet',zorder=5)


    for i in range(len(cnv_final)):
        x = np.array(range(int(cnv_final[i][1]), int(cnv_final[i][2])))
        y = [int(cnv_final[i][3])] * len(x)
        if int(cnv_final[i][3]) == 1:
            plt.plot(x, y, color='grey', linewidth=1.5)
        elif int(cnv_final[i][3]) == 2:
            plt.plot(x, y, color='blue', linewidth=1.2)
        elif int(cnv_final[i][3]) == 3:
            plt.plot(x, y, color='orange', linewidth=1.5)
        elif int(cnv_final[i][3]) >= 4:
            plt.plot(x, y, color='red', linewidth=1.5)
#    for i in range(len(cnv_final)):
#        x = np.array(range(int(cnv_final[i][1]), int(cnv_final[i][2])))
#        if int(cnv_final[i][3]) == 1:
#            plt.axhline(y=1,  xmin = int(cnv_final[i][1]), xmax= int(cnv_final[i][2]), color='grey', linewidth=1.5)
#        elif int(cnv_final[i][3]) == 2:
#            plt.axhline(y=2,  xmin = int(cnv_final[i][1]), xmax= int(cnv_final[i][2]), color='royalblue', linewidth=1.5)
#        elif int(cnv_final[i][3]) == 3:
#            plt.axhline(y=3,  xmin = int(cnv_final[i][1]), xmax= int(cnv_final[i][2]), color='orange', linewidth=1.5)
#        elif int(cnv_final[i][3]) >= 4:
#            plt.axhline(y=4,  xmin = int(cnv_final[i][1]), xmax= int(cnv_final[i][2]), color='red', linewidth=1.5)


# set ticks and tick labels
    ax.set_xticks(np.arange(18000000,22500000,500000))
    ax.set_xticklabels(np.arange(18,22.5,0.5))

    #ax.set_xlabel('Chromosomal position (Mb)', fontsize=12, fontweight='roman')
    ax.set_ylabel('Copy number', fontsize=12, fontweight='roman')
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax.yaxis.set_major_locator(ticker.MultipleLocator(1.00))
#ax.annotate('Chr11', xy=(.925, 0.95), xycoords='figure fraction', horizontalalignment='left', verticalalignment='top', fontsize=11, fontweight='roman')

# Hide the right and top spines
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
# Move left and bottom spines outward
    ax.spines['left'].set_position(('outward', 10))
    ax.spines['bottom'].set_position(('outward', 1))

    plt.tight_layout()
    plt.savefig('chr1-2.pdf', format ='pdf', transparent=True)
