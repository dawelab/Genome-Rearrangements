
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 26 17:00:40 2018
@author: JianingLiu
"""

import numpy as np
import random
from functools import reduce
import pandas as pd

# Read fasta file 
def open_fasta(file):
    dna = ""
    with open(file) as seq:
        for lines in seq:
            if lines[0] != ">":
                dna += lines.strip()
    return dna
# Output modified fasta file 
def write_sim(seq, file, title):
    with open(file, "w") as out:
        out.write(">Chr" + str(title) + "\n" + seq + "\n")
# Invert DNA sequence        
def invert(DNA_sequence):
    inverted_seq = ""        
    for x in DNA_sequence.upper():
        if x == "A":
            inverted_seq += "T"
        elif x == "C":
            inverted_seq += "G"
        elif x == "G":
            inverted_seq += "C"
        elif x == "T":
            inverted_seq += "A"
        else:
            inverted_seq += "N"
    return inverted_seq[::-1]
# Probability to define whether sequence should be inverted
def prob(p):
    x = random.uniform(0, 1)
    if x > float(p):
        return False
    else:
        return True

#def deletion(seq, start,end): 
#        
#    return seq.replace(seq[start:end],'')
#
#def duplication(seq,num,start,end): 
#    
#    return seq[:start]+ seq[start:end]*int(num) + seq[end:]

    
# Break lambda DNA into pieces and invert half of the pieces 
def shatter(seq,num): # number is the breakpont number  
#    seq = open_fasta('/Volumes/WD2T/RSVsim/NCBI_lambda_genome.fa')
#    num = 5
    coordinates = []
    isolates = np.sort(np.random.randint(0, len(seq), size=int(num)))
    coordinates = [0] + list(isolates) + [len(seq)]
    coordinate_pairs = []
    for i in range(len(coordinates)-1):
        coordinate_pairs.append([coordinates[i],coordinates[i+1]])
    random.shuffle(coordinate_pairs, random.random)
    coordinate_pairs_reversed = []
    for x in coordinate_pairs:
        if prob(0.5):
            x.reverse()
            coordinate_pairs_reversed.append(x)
        else: 
            coordinate_pairs_reversed.append(x)
    random.shuffle(coordinate_pairs_reversed, random.random) 
    return coordinate_pairs_reversed 

# Ligate the 100 broken pieces generatged by shatter function above to generate 50 pieces
def ligate(seq,lambda_pieces):
#    seq = open_fasta('/Volumes/WD2T/RSVsim/NCBI_lambda_genome.fa')
    coordinate_pairs = shatter(seq,lambda_pieces-1)# lambda_pieces=49+1
#    coordinate_pairs_2 = shatter(seq,3)
#    coordinate_pairs_3 = shatter(seq,2)
#    coordinate_pairs = coordinate_pairs_1 
#    coordinate_pairs = coordinate_pairs_1 + coordinate_pairs_2 + coordinate_pairs_3 # 103 pieces
    coordinate_pairs_reversed_sequence = []
    for i in range(0,len(coordinate_pairs)):
        if coordinate_pairs[i][0] < coordinate_pairs[i][1]:
            coordinate_pairs_reversed_sequence.append(seq[coordinate_pairs[i][0]:coordinate_pairs[i][1]])
        else:
            coordinate_pairs_reversed_sequence.append(invert(seq[coordinate_pairs[i][1]:coordinate_pairs[i][0]]))  
    x = list(range(1,lambda_pieces))
    random.shuffle(x)
    cut = [0] + list(np.sort(list(x[0:int(lambda_pieces/2)-1]))) + [lambda_pieces+1] # 49, 103
    ligated_coordinate_pairs = []
    sequence_final = []
    trans_intra = []
    for i in range(0, len(cut)-1):
        individual = coordinate_pairs[cut[i]:cut[i+1]]
        ligated_coordinate_pairs.append(individual)
        sub = sum(individual, [])[1:-1]   
        trans_intra += [sub[x:x+2] for x in range(0, len(sub), 2)]
        sequence_final.append(''.join(str(r) for r in coordinate_pairs_reversed_sequence[cut[i]:cut[i+1]]))
    print(trans_intra)
    return ligated_coordinate_pairs, sequence_final, trans_intra

# Generate coordinates in genome to place del, dup and ins
def location(seq,dup_num, del_num, ins_num):
    num = int(dup_num) + int(del_num) 
    x = list(range(1,int(dup_num) + int(del_num) + int(ins_num))) # 103
    random.shuffle(x)
    coord = np.sort((np.random.beta(2,2,num*2 + int(ins_num))*len(seq)).astype(int)) # generate a list of coordinates of beta distribution 
    a = np.array(range(0,len(coord)-1))
    random.shuffle(a)
    coord2 = np.sort(a[:num])
    dup_del = []
    for i in coord2:
        dup_del.append([coord[i], coord[i+1]])
    random.shuffle(dup_del)
    dupco = np.sort(dup_del[:int(dup_num)],axis=0)
    delco = np.sort(dup_del[int(dup_num):],axis=0)
    s = dupco.tolist() + delco.tolist()
    if len(s)>0:
        s2 = reduce(lambda x,y: x+y,s)
    else:
        s2 = []
    insco = []
    for y in coord:
        if y not in s2:
            insco.append([y])
    dupco = dupco.tolist()
    delco = delco.tolist()
    return dupco, delco, insco

# Create modified sequnece by inserting del, dup and ins at the same time
def del_dup_ins(seq_genome, dupco, delco, insco,ligated,ligated_coordinate_pairs): 
    trans = []
    total = delco + dupco + insco
    total.sort(key=lambda x: x[0])
    final = seq_genome[:total[0][0]]
    c = 0 # counter for counting the number of double pairs
    for i in range(0, len(total)-1):      
        if len(total[i]) == 2:
            c += 1
            if total[i] in delco:
                final += seq_genome[total[i][1]:total[i+1][0]]
                
            else:
                final += 2*seq_genome[total[i][0]:total[i][1]] + seq_genome[total[i][1]:total[i+1][0]]
        if len(total[i]) == 1: # Insertion
            final += str(ligated[i-c]) + seq_genome[total[i][0]:total[i+1][0]]
            trans.append([ligated_coordinate_pairs[i-c][0][0],int(total[i][0])])
            trans.append([ligated_coordinate_pairs[i-c][-1][1], int(total[i][0])])
    if len(total[-1]) == 2: 
        final += seq_genome[total[-1][1]:]
    else:
        trans.append([ligated_coordinate_pairs[-1][0][0],int(total[-1][0])])
        trans.append([ligated_coordinate_pairs[-1][-1][1],int(total[-1][0])])
        final += str(ligated[-1]) + seq_genome[total[-1][0]:]
    return final, trans



def sv(delco, dupco, trans, trans_intra):
    chr_dupco = pd.DataFrame({'Chromosome':[1]*len(dupco)})
    type_dupco = pd.DataFrame({'SV':['DUP']*len(dupco)})
    dupco = pd.DataFrame(dupco, columns = ['coordinate_1', 'coordinate_2'])
    combined_dup = pd.concat([chr_dupco, dupco, type_dupco], axis=1)
    chr_delco = pd.DataFrame({'Chromosome':[1]*len(delco)})
    type_delco = pd.DataFrame({'SV':['DEL']*len(delco)})
    delco = pd.DataFrame(delco, columns = ['coordinate_1', 'coordinate_2'])
    combined_delco = pd.concat([chr_delco, delco, type_delco], axis=1)
    chr1_trans = pd.DataFrame({'Chromosome':[1]*len(trans)})
    chr2_trans = pd.DataFrame({'Chromosome':['Chr3']*len(trans)})
    type_trans = pd.DataFrame({'SV':['Trans_inter']*len(trans)})
    trans = pd.DataFrame(trans, columns = ['Chr_2', 'Chr_1'])
    trans1 = trans['Chr_1']
    trans2 = trans['Chr_2']
    combined_trans = pd.concat([chr1_trans, trans1, chr2_trans, trans2, type_trans], axis=1)    
    chr_trans_intra = pd.DataFrame({'Chromosome':['Chr3']*len(trans_intra)})
    type_trans_intra = pd.DataFrame({'SV':['Trans_intra']*len(trans_intra)})
    trans_intra = pd.DataFrame(trans_intra, columns = ['Chr_1', 'Chr_2'])
    trans_intra1 = trans_intra['Chr_1']
    trans_intra2 = trans_intra['Chr_2']
    combined_trans_intra = pd.concat([chr_trans_intra, trans_intra1, chr_trans_intra, trans_intra2, type_trans_intra], axis=1)
    trans_combined = pd.concat([combined_trans_intra, combined_trans])
    indel_combined = pd.concat([combined_dup, combined_delco])
    trans_combined.index = trans_combined.index + 1
    indel_combined.index = indel_combined.index + 1
    return trans_combined, indel_combined
def main():
	lambda_seq = open_fasta(file1)
	copy = 1
	breaks_per_copy = 100
	seq = lambda_seq * copy
	lambda_insertion_fragments = breaks_per_copy * copy
	genome_seq = open_fasta(file2)
	dupco = location(genome_seq, 0,0,lambda_insertion_fragments/2)[0]
	delco = location(genome_seq, 0,0,lambda_insertion_fragments/2)[1]
	insco = location(genome_seq, 0,0,lambda_insertion_fragments/2)[2]
	ligated = ligate(seq,lambda_insertion_fragments)[1]
	ligated_coordinate_pairs = ligate(seq,lambda_insertion_fragments)[0]
	trans_intra = ligate(seq,lambda_insertion_fragments)[2]
	final = del_dup_ins(genome_seq, dupco, delco, insco,ligated,ligated_coordinate_pairs)[0]
	trans = del_dup_ins(genome_seq, dupco, delco, insco,ligated,ligated_coordinate_pairs)[1]
	trans_combined = sv(delco, dupco, trans, trans_intra)[0]
	indel_combined = sv(delco, dupco, trans, trans_intra)[1]
	#indel_combined.to_csv('/Volumes/WD2T/RSVsim/Indel.csv', header=True, index = True) 
	trans_combined.to_csv(trans_file, header=True, index = True)  
	write_sim(final, path+str(copy)+filename, 'derivative')
	with open(path+str(copy)+filename, "w") as out:
    	for x in range(len(ligated)):
        	out.write(">Fragment" + str(x+1) + "\n" + ligated[x]+ "\n")


'''
Pseudo code:
    
1. read lambda file
2. shatter lambda and generate lambda piece arrays with function shatter
3. ligate lambda and generate lambda piece arrays with function ligate
4. read genome file
5. generate a list of coordinates as candidates for insertions, duplications/deletions 
6. Ligate lambda and genome together
7. Output sv file and final derivative genome
'''
