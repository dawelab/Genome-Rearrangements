#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 29 21:13:28 2018

@author: JianingLiu
"""

import pandas as pd 
import numpy as np
import sys
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

def load_file(file):
    blast = pd.read_csv(file, delimiter="\t",header=None)
    return blast

def filter(blast):
    tmp = []
    blast = blast[blast[3]>25] # filter the alignment size < 30bp
    # blast = blast[blast[4]<1e-5] # filter the evalue higher than 1e-5
    blast = blast.sort_values(by=[3],ascending=False).sort_values(by=[0]) #Sort by the node name and alignment length
    for i in range(len(blast.index)):
        for j in range(i+1,len(blast.index)):
            if blast.iloc[i,0] == blast.iloc[j,0]:
                # if row i and row j have >20bp overlap in the case of coordinates, then remove row j (in this case, the alignment len of j is smaller than that of i)
                if len(set(range(blast.iloc[i,6], blast.iloc[i,7])).intersection(range(blast.iloc[j,6], blast.iloc[j,7]))) > 18: 
                    tmp.append(j)
    blast = blast.drop(blast.index[tmp]) 
    return blast.sort_values(by=[0,6])

def extract(filtered):
    orientation = []
    relative_orientation = []
    overlap = []
    for i in range(len(filtered.index)):
        if filtered.iloc[i,8] < filtered.iloc[i,9]:
            orientation.append("HT")
        else:
            orientation.append("TH")
    for i in range(len(filtered.index)-1):
        if filtered.iloc[i,0] == filtered.iloc[i+1,0]:
            relative_orientation.append(str(orientation[i][1]) + str(orientation[i+1][0]))
            overlap.append(filtered.iloc[i+1,6]-filtered.iloc[i,7]-1)
    size = filtered[3].tolist()
    return relative_orientation, overlap, size

#def plot_orientation(orientation): 
#    x = ["HT", "TH", "HH", "TT"]
#    y = [orientation.count("HT")/len(orientation), orientation.count("TH")/len(orientation), orientation.count("HH")/len(orientation), orientation.count("TT")/len(orientation)]
#    ax = fig.add_subplot(131)
#    plt.bar(x,y, color="darkblue",width=0.5,alpha=0.7)
#    plt.rc('xtick',labelsize=13)
#    plt.rc('ytick',labelsize=12)
#    ax.set_title('Distribution of relative orientations between junction fragments', fontsize=14)
#    ax.set_ylabel('Frequency', fontsize=14)
#    ax.spines['right'].set_visible(False)
#    ax.spines['top'].set_visible(False)
#    ax.spines['left'].set_visible(False)
#    ax.set_axisbelow(True)
#    ax.xaxis.set_ticks_position('none') 
#    ax.yaxis.set_ticks_position('none') 
#    ax.yaxis.grid(linestyle='-', linewidth=0.3)
#    
#def plot_size(size):
#    ax1 = fig.add_subplot(132)
#    plt.hist(size, normed=True, bins=30, color="teal",alpha=0.7)
#    plt.rc('xtick',labelsize=13)
#    plt.rc('ytick',labelsize=12)
#    ax1.set_title('Distribution of the sizes of junction fragments', fontsize=14)
#    ax1.set_ylabel('Frequency', fontsize=14)
#    ax1.spines['right'].set_visible(False)
#    ax1.spines['top'].set_visible(False)
#    ax1.spines['left'].set_visible(False)
#    ax1.set_axisbelow(True)
#    ax1.xaxis.set_ticks_position('none') 
#    ax1.yaxis.set_ticks_position('none') 
#    ax1.yaxis.grid(linestyle='-', linewidth=0.3)
#    
#    
#def plot_overlap(overlap):  
#    x = np.array(overlap)
#    blunt = len(x[x==0])
#    insert = len(x[x>0])
#    nhej = sum(np.logical_and(x<0, x>-5))
#    mmej = sum(np.logical_and(x>=-25, x<=-5)) # 5-25
#    ssa = len(x[x<-25])
#    junction = ['', "Blunt ends", "Insertion of novel sequence", "1-4bp microhomology", "5-25bp microhomology",  ">25bp microhomology"]
#    ax2 = fig.add_subplot(133)
#    plt.bar(np.arange(5),np.array([blunt,insert,nhej,mmej,ssa])/len(x), color="darkblue",width=0.5,alpha=0.7)
#    plt.rc('xtick',labelsize=14)
#    plt.rc('ytick',labelsize=12)
#    ax2.set_title('Distribution of relative orientations between junction fragments', fontsize=14)
#    ax2.set_xticks(np.arange(len(junction)))
#    ax2.set_xticklabels(junction, rotation=-25, ha="left")
#    ax2.set_ylabel('Frequency', fontsize=13)
#    ax2.spines['right'].set_visible(False)
#    ax2.spines['top'].set_visible(False)
#    ax2.spines['left'].set_visible(False)
#    ax2.set_axisbelow(True)
#    ax2.xaxis.set_major_locator(MaxNLocator(integer=True))
#    ax2.xaxis.set_ticks_position('none') 
#    ax2.yaxis.set_ticks_position('none') 
#    ax2.yaxis.grid(linestyle='-', linewidth=0.3)
#blast = load_file("/Users/JianingLiu/Desktop/RiceE24-spades-pe.txt")
#filtered = filter(blast)
#orientation = extract(filtered)[0]
#plot_orientation(orientation)
#overlap = extract(filtered)[1]
#plot_overlap(overlap)
#size = extract(filtered)[2] 
#plot_orientation(size)
def main():
    file = sys.argv[1]
    output = sys.argv[2]
    csv_file = sys.argv[3]
    blast = load_file(file)
    filtered = filter(blast)
    filtered.to_csv(csv_file)
    orientation = extract(filtered)[0]
    overlap = extract(filtered)[1] 
    size = extract(filtered)[2] 
    fig = plt.figure(figsize=(24,8),facecolor='w')
    ax = fig.add_subplot(131)
    x = ["HT", "TH", "HH", "TT"]
    y = [orientation.count("HT")/len(orientation), orientation.count("TH")/len(orientation), orientation.count("HH")/len(orientation), orientation.count("TT")/len(orientation)]
    plt.bar(x,y, color="darkblue",width=0.5,alpha=0.7)
    plt.rc('xtick',labelsize=13)
    plt.rc('ytick',labelsize=12)
    ax.set_title('Distribution of relative orientations between junction fragments', fontsize=14)
    ax.set_ylabel('Frequency', fontsize=14)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.set_axisbelow(True)
    ax.xaxis.set_ticks_position('none') 
    ax.yaxis.set_ticks_position('none') 
    ax.yaxis.grid(linestyle='-', linewidth=0.3)
    ax1 = plt.subplot(132)
    plt.hist(size, normed=True, bins=30, color="teal",alpha=0.7)
    plt.rc('xtick',labelsize=14)
    plt.rc('ytick',labelsize=12)
    ax1.set_title('Distribution of the sizes of junction fragments', fontsize=14)
    ax1.set_ylabel('Frequency', fontsize=13)
    ax1.spines['right'].set_visible(False)
    ax1.spines['top'].set_visible(False)
    ax1.spines['left'].set_visible(False)
    ax1.set_axisbelow(True)
    ax1.xaxis.set_ticks_position('none') 
    ax1.yaxis.set_ticks_position('none') 
    ax1.yaxis.grid(linestyle='-', linewidth=0.3)
    x = np.array(overlap)
    blunt = len(x[x==0])
    insert = len(x[x>0])
    nhej = sum(np.logical_and(x<0, x>-5))
    mmej = sum(np.logical_and(x>=-25, x<=-5)) # 5-25
    ssa = len(x[x<-25])
    junction = ['', "Blunt ends", "Insertion of novel sequence", "1-4bp microhomology", "5-25bp microhomology",  ">25bp microhomology"]
    ax2 = fig.add_subplot(133)
    plt.bar(np.arange(5),np.array([blunt,insert,nhej,mmej,ssa])/len(x), color="darkblue",width=0.5,alpha=0.7)
    plt.rc('xtick',labelsize=14)
    plt.rc('ytick',labelsize=12)
    ax2.set_title('Distribution of relative orientations between junction fragments', fontsize=14)
    ax2.set_xticks(np.arange(len(junction)))
    ax2.set_xticklabels(junction, rotation=-25, ha="left")
    ax2.set_ylabel('Frequency', fontsize=13)
    ax2.spines['right'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax2.spines['left'].set_visible(False)
    ax2.set_axisbelow(True)
    ax2.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax2.xaxis.set_ticks_position('none') 
    ax2.yaxis.set_ticks_position('none') 
    ax2.yaxis.grid(linestyle='-', linewidth=0.3)
    plt.tight_layout()
    plt.savefig(output)
if __name__ == "__main__":
    main()
