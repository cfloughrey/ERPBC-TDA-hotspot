# -*- coding: utf-8 -*-
"""
Created on Thu Mar 24 16:15:42 2022

@author: ciara
"""


import seaborn as sns
import matplotlib.pyplot as plt



#---------------------------------read in files----------------------------#
#read in the dataframe generated from hotspot search describing the results for each parameter combination
h_df = "metabric/hotspot_dataframe.csv" 



#-----------------Run through survival for different parameters ------------#
#list of intervals and list of overlap
interval_list = np.arange(10,32,2)
overlap_list = np.arange(10,50,5)

param_list = [np.zeros(len(interval_list)) ]* len(overlap_list)
#dataframe with na is more accurate than 0 
param_df = pd.DataFrame( columns = interval_list, index = overlap_list)


#for each interval and overlap, assign the survival p-value to the matrix 
interval = h_df["interval"]
overlap = h_df["overlap"]
prog = h_df["logrank"]

for i in enumerate(interval):
    o = int(overlap.iloc[i[0]] * 100)

    #if multiple hotspots to same parameters, keep hotspot with lower p-value
    if pd.isnull(param_df.loc[45,10]):
        param_df.loc[o,i[1]] = prog.iloc[i[0]]
    elif prog.iloc[i[0]] < param_df.loc[o,i[1]]:
        param_df.loc[o,i[1]] = prog.iloc[i[0]]

param_df = round(param_df,3)
mask=param_df.isnull()
param_df = param_df.fillna(0)



#---------------------------plot results--------------------------------------#
#The smaller the p-value the greate the significance 
# plot the heatmap
fig, ax = plt.subplots(figsize=(20, 10))
sns.heatmap(param_df,xticklabels=param_df.columns, yticklabels=param_df.index,
            cmap = "Blues_r",annot=True, cbar=True, linewidths=.3, mask = mask,
            cbar_kws={'label': 'log-rank p-value'})
#plt.title("Survival significance of Hotspot", fontsize=40)
plt.xlabel('Intervals', fontsize = 25) # x-axis label with fontsize 15
plt.ylabel('Overlap', fontsize = 25) # y-axis label with fontsize 15
ax.set_facecolor("lightgrey")
plt.show()