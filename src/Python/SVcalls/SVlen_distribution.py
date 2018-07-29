#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 10 23:26:48 2018

@author: JianingLiu
"""
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import statsmodels.api as sm
from statsmodels.formula.api import ols
from scipy import stats
circular = pd.read_csv("circular_cnv.csv", sep=',')

x1 = circular[circular['Events'].str.contains('10A-1')].Size
x2 = circular[circular['Events'].str.contains('10A-2')].Size
x3 = circular[circular['Events'].str.contains('10A-3')].Size
x4 = circular[circular['Events'].str.contains('10A-4')].Size
x5 = circular[circular['Events'].str.contains('10A-5')].Size
x6 = circular[circular['Events'].str.contains('10A-6')].Size
x7 = circular[circular['Events'].str.contains('12A-1')].Size
x8 = circular[circular['Events'].str.contains('12A-2')].Size
x9 = circular[circular['Events'].str.contains('12A-3')].Size
x10 = circular[circular['Events'].str.contains('12A-4')].Size
x11 = circular[circular['Events'].str.contains('12A-5')].Size
x12 = circular[circular['Events'].str.contains('12A-6')].Size
simple = x1.append(x2).append(x3)
complexe =  x4.append(x5).append(x6)
pval1 = stats.f_oneway(simple,complexe)[1]
simple2 = x7.append(x8).append(x9)
complexe2 =  x10.append(x11).append(x12)
pval2 = stats.f_oneway(simple2,complexe2)[1]
pval_simple = stats.f_oneway(x1,x2,x3)[1]
pval_simple2 = stats.f_oneway(x7,x8,x9)[1]
pval_complexe = stats.f_oneway(x4,x5,x6)[1]
pval_complexe2 = stats.f_oneway(x10,x11,x12)[1]

#moore_lm = ols('Size ~ C(Events, Sum)*C(CNV, Sum)', data=circular).fit()
#table = sm.stats.anova_lm(moore_lm, typ=2) # Type 2 ANOVA DataFrame
circular['Size'] = circular['Size'] /1000000

'''
PLOT
'''
sns.set_style("whitegrid")
plt.figure(figsize=(30,10))
#ax = sns.violinplot(x="Events", y="Size", hue="CNV",data=circular, palette="Set3")
ax = sns.violinplot(x="Events", y="Size", hue="CNV",data=circular, palette="Set2")

#ax.set_yscale('log2')
ax = sns.swarmplot(x="Events", y="Size", hue="CNV",data=circular,split=True)
plt.xticks(fontsize=22,rotation=-30)
plt.yticks(fontsize=22)
plt.ylabel("Size (Mb)",fontsize=26, rotation=90)
plt.xlabel("Events",fontsize=26)
plt.legend(fontsize=18)
#x1, x2,x3 = -0.4, 1, 2.0   # columns 'Sat' and 'Sun' (first column: 0, see plt.xticks())
#y, h, col = circular['Size'].max()-4, 1, 'k'
#plt.plot([x1, x1,x1,x2,x2,x2,x2,x2,x2,x3,x3,x3], [y,y+h, y+h, y+h, y+h,y+h, y+h,y+h, y+h,y+h, y+h, y], lw=1.5, c=col)
#x4, x5,x6 = 2.6, 4, 5.15   # columns 'Sat' and 'Sun' (first column: 0, see plt.xticks())
#y, h, col = circular['Size'].max()-0.2, 1, 'k'
#plt.plot([x4, x4, x5,x5, x5, x5,x5, x5, x5, x5,x6,x6], [y,y+h, y+h, y+h, y+h,y+h, y+h,y+h, y+h,y+h, y+h, y], lw=1.5, c=col)
#plt.text((x4+x5+x6)*.35, y+h+0.1, "p=0.041927", ha='center', va='bottom', color=col,fontsize=22)
#x1, x2,x3 = 5.6, 7, 8.0   # columns 'Sat' and 'Sun' (first column: 0, see plt.xticks())
#y, h, col = circular['Size'].max()-4, 1, 'k'
#plt.plot([x1, x1,x1,x2,x2,x2,x2,x2,x2,x3,x3,x3], [y,y+h, y+h, y+h, y+h,y+h, y+h,y+h, y+h,y+h, y+h, y], lw=1.5, c=col)
#x4, x5,x6 = 8.6, 9, 11.15   # columns 'Sat' and 'Sun' (first column: 0, see plt.xticks())
#y, h, col = circular['Size'].max()-0.2, 1, 'k'
#plt.plot([x4, x4, x5,x5, x5, x5,x5, x5, x5, x5,x6,x6], [y,y+h, y+h, y+h, y+h,y+h, y+h,y+h, y+h,y+h, y+h, y], lw=1.5, c=col)
#plt.text((x4+x5+x6)*.35, y+h+0.1, "p=0.019771", ha='center', va='bottom', color=col,fontsize=22)
x1 = 6
y, h, col = circular['Size'].max()-4, 1, 'k'
plt.plot([x1], [y+h], lw=1.5, c=col)

