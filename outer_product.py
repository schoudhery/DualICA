from operator import methodcaller
import pandas as pd 
import sys
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import statsmodels.api as sm
import statsmodels.stats
import scipy.stats as scistat
from sklearn import linear_model
from sklearn.metrics import mean_squared_error

# python3 linear_reg.py G_matrix C_Matrix lfcs_data

# G_matrix : The signal matrix of the ICA where genes are the rows

# C_matrix : The signal matrix of the ICA where conditions  are the rows

# lfcs data : genes x conditions and each value is the LFC of the condition compared to DMSO.

# categories data : a metadata file with atelast two colums- Name and Classification, where name is the condition names in lfcs data and classification is the class of the condition name. 
# NOTE: categories data file is used ONLY to label the output file for the user and not used in the calculation of the blocks. 


fig, (ax1, ax2) = plt.subplots(1, 2, figsize =(40,40),tight_layout = True)
G_matrix = pd.read_csv(sys.argv[1], index_col=0, low_memory=False)
C_matrix = pd.read_csv(sys.argv[2], index_col=0, low_memory=False)
split_index = C_matrix.index.str.split(" ", expand = True)

lfcs_data = pd.read_csv(sys.argv[3], sep="\t", index_col=0,low_memory=False)
categories_data = pd.read_csv(sys.argv[4], sep=",")
categories_data = categories_data.loc[categories_data["Name"].isin(lfcs_data.columns), :]
categories_data =categories_data.sort_values(by=['Classification'])

lfcs_data= lfcs_data[categories_data["Name"]]

LIM=6
new_cols =[]
for col in lfcs_data.columns:
    classification = categories_data[categories_data["Name"]==col]["Classification"].iloc[0]
    new_cols.append(col+" ("+classification+")")

    lfcs_data[col] = (lfcs_data[col] - lfcs_data[col].mean())/lfcs_data[col].std()
    lfcs_data[col][lfcs_data[col] >= LIM] = LIM
    lfcs_data[col][lfcs_data[col] <= -LIM] = -LIM 
lfcs_data.columns = new_cols

lfcs_data = lfcs_data.replace([np.inf,-np.inf],np.nan)
lfcs_data = lfcs_data.dropna()


#### ------------------------------ Outer Produxt
sums_matrix= []
for i in G_matrix.columns:
    sum_row=[]
    G_col = G_matrix[i]
    for j in C_matrix.columns:
        C_col = C_matrix[j]
        outer_result = np.outer(G_col, C_col)
        #print(outer_result[0])
        inner_result= np.multiply(outer_result, lfcs_data.to_numpy())
        result = inner_result.sum()
        sum_row.append(result)
    sums_matrix.append(sum_row)

Sums_df = pd.DataFrame(sums_matrix, index = G_matrix.columns, columns=C_matrix.columns)
s1=sns.heatmap(Sums_df, ax=ax1, cmap=sns.color_palette("vlag", as_cmap=True), center=0, linewidth=0.5, cbar_kws={"aspect":"100"},xticklabels=True, annot=False, fmt='.2f')
s1.set(xlabel='Condition Components', ylabel='Gene Components', title= " ")
Sums_df.to_csv("outer_prods.csv")