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
    blast = blast[blast[3]>30] # filter the alignment size < 30bp
    #blast = blast[blast[10]<1e-1] # filter the evalue higher than 1e-5
    blast = blast.groupby([1]).apply(lambda x: x.sort_values(by=[3],ascending=False)).reset_index(drop=True) #Sort by the node name and alignment length
    for i in range(len(blast.index)):
        for j in range(i+1,len(blast.index)):
            if blast.iloc[i,0] == blast.iloc[j,0]:
                # if row i and row j have >20bp overlap in the case of coordinates, then remove row j (in this case, the alignment len of j is smaller than that of i)
                if len(set(range(blast.iloc[i,8], blast.iloc[i,9])).intersection(range(blast.iloc[j,8], blast.iloc[j,9]))) > 20: 
                    tmp.append(j)
    blast = blast.drop(blast.index[tmp]) 
    blast = blast.sort_values(by=[1,8,9])
    return blast

#def extract(filtered):
#    orientation = []
#    relative_orientation = []
#    overlap = []
#    for i in range(len(filtered.index)):
#        if filtered.iloc[i,8] < filtered.iloc[i,9]:
#            orientation.append("HT")
#        else:
#            orientation.append("TH")
#    for i in range(len(filtered.index)-1):
#        if filtered.iloc[i,0] == filtered.iloc[i+1,0]:
#            relative_orientation.append(str(orientation[i][1]) + str(orientation[i+1][0]))
#            overlap.append(filtered.iloc[i+1,6]-filtered.iloc[i,7]-1)
#    size = filtered[3].tolist()
#    return relative_orientation, overlap, size
def main():    
    file = sys.argv[1]
    output = sys.argv[2]
    x = load_file(file)
    filtered = filter(x)
    filtered.to_csv(output,header=False)
if __name__ == '__main__':
  main()   
