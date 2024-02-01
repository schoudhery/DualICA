from sklearn.decomposition import FastICA
import pandas as pd
import numpy as np
import sys
import scipy
import seaborn as sns
import matplotlib.pyplot as plt


# python3 dual_ICA.py lfcs_data categories_data NCIC NGIC
# lfcs data : genes x conditions and each value is the LFC of the condition compared to DMSO.

# categories data : a metadata file with atelast two colums- Name and Classification, where name is the condition names in lfcs data and classification is the class of the condition name. 
# NOTE: categories data file is used ONLY to label the output file for the user and not used in the calculation of the blocks. 

#NCIC : number of compoenents to extract from conditions -> signal matrix of this has genes on the rows

#NGIC : number of compoenents to extract from genes -> signal matrix of this has conditions on the rows

lfcs_data = pd.read_csv(sys.argv[1], sep="\t", index_col=0)
lfcs_data = lfcs_data.loc[:, ~lfcs_data.columns.isin(["sentinel", "nTA"])]

categories_data = pd.read_csv(sys.argv[2], sep=",")
categories_data = categories_data.loc[categories_data["Name"].isin(lfcs_data.columns), :]
categories_data =categories_data.sort_values(by=['Classification'])

lfcs_data= lfcs_data[categories_data["Name"]]


LIM = 6 # limit lfcs beyond +/-6

new_cols =[]
for col in lfcs_data.columns:
    print(col,len(lfcs_data[col][lfcs_data[col] >= LIM]), len(lfcs_data[col][lfcs_data[col] <= -LIM]))
    classification = categories_data[categories_data["Name"]==col]["Classification"].iloc[0]
    new_cols.append(col+" ("+classification+")")

    lfcs_data[col] = (lfcs_data[col] - lfcs_data[col].mean())/lfcs_data[col].std()
    lfcs_data[col][lfcs_data[col] >= LIM] = LIM
    lfcs_data[col][lfcs_data[col] <= -LIM] = -LIM 
lfcs_data.columns = new_cols

lfcs_data = lfcs_data.replace([np.inf,-np.inf],np.nan)
lfcs_data = lfcs_data.dropna()

NGENES = len(lfcs_data.index)
NCONDS = len(lfcs_data.columns)

NCICs = int(sys.argv[3]) #100
NGICs = int(sys.argv[4]) #70

ica1 = FastICA(n_components=NCICs, max_iter=5000)
G = pd.DataFrame(ica1.fit_transform(lfcs_data))  # Reconstruct signals
G.index = lfcs_data.index

A1_ = pd.DataFrame(ica1.mixing_)  # Get estimated mixin g matrix
A1_.index = lfcs_data.columns


ica2 = FastICA(n_components=NGICs)
C = pd.DataFrame(ica2.fit_transform(lfcs_data.T))  # Reconstruct signals
C.index = lfcs_data.columns

A2_ = pd.DataFrame(ica2.mixing_)  # Get estimated mixing matrix
A2_.index = lfcs_data.index


# clustering conditions
cond_set_assignments ={}
for col in lfcs_data.columns:
    cond_set_assignments[col]=[]

cond_modules={}
cond_modules["9999"]=[]

for i in C.columns:
    cond_modules[str(i)+"-"]=[]
    cond_modules[str(i)+"+"]=[]
    current_list = C[i]
    stat,pval = scipy.stats.normaltest(current_list, axis=0, nan_policy='propagate')
    cond_set= []
    while pval<0.05:
        max_idx = current_list.abs().idxmax()
        if current_list[max_idx]<0:
            curr_i = str(i)+"-"
        else:
            curr_i = str(i)+"+"

        cond_set.append(max_idx)
        current_list = current_list.drop(labels = [max_idx]) 
        cond_modules[curr_i] =  cond_modules[curr_i] +[max_idx]
        cond_set_assignments[max_idx] = cond_set_assignments[max_idx]+[curr_i]
        stat,pval = scipy.stats.normaltest(current_list, axis=0, nan_policy='propagate')

#assign remaining conditions
for c in cond_set_assignments:
    if cond_set_assignments[c]==[]: 
        max_idx = C.loc[c,:].abs().idxmax()
        if C.loc[c,:][max_idx] <0:
            curr_i =str(max_idx)+ "-"
        else:
            curr_i =str(max_idx)+ "+"
        cond_set_assignments[c].append(curr_i)
        cond_modules[curr_i] =  cond_modules[curr_i] +[c]

#clustering genes
gene_set_assignments ={}
for gene in lfcs_data.index:
    gene_set_assignments[gene]=[]

gene_modules = {}
gene_modules["9999"]=[]

for i in G.columns:
    gene_modules[str(i)+"-"]=[]
    gene_modules[str(i)+"+"]=[]
    #print(i)
    current_list = G[i]
    stat,pval = scipy.stats.normaltest(current_list, axis=0, nan_policy='propagate')
    gene_set= []
    main_stat = stat
    stat_list=[]
    while pval<0.05:
        stat_list.append(stat)
        if len(stat_list)>=10:
            test_list = stat_list[-10:]
            res = all(i < j for i, j in zip(test_list, test_list[1:]))
            if res:
                gene_set = gene_set[: len(gene_set) - 10]
                break
        max_idx = current_list.abs().idxmax()
        if current_list[max_idx]<0:
            curr_i = str(i)+"-"
        else:
            curr_i = str(i)+"+"
        gene_set.append(max_idx)
        current_list = current_list.drop(labels = [max_idx])

        gene_modules[curr_i] =  gene_modules[curr_i] +[max_idx]
        gene_set_assignments[max_idx] = gene_set_assignments[max_idx]+[curr_i]
        stat,pval = scipy.stats.normaltest(current_list, axis=0, nan_policy='propagate')

no_assigned= []
multiple_assigned=[]
for gene in gene_set_assignments:
    if len(gene_set_assignments[gene])==0:
        no_assigned.append(gene)
        #assign remaining conditions
        gene_set_assignments[gene]=["9999"]
        gene_modules["9999"] = gene_modules["9999"]+["gene"]
    if len(gene_set_assignments[gene])>1:multiple_assigned.append(gene)

df_list= []
for i in gene_modules:
    print(i, len(gene_modules[i]))
    sub_data = lfcs_data[lfcs_data.index.isin(gene_modules[i])]
    sub_data.insert(0, 'Module', i)
    df_list.append(sub_data)

out_df = pd.concat(df_list)
df_list = []
for i in cond_modules:
    sub_data = out_df[cond_modules[i]]
    cond_idx = ([i]*len(cond_modules[i]))
    cond_row = pd.DataFrame([cond_idx], index=['Condition Modules'], columns=cond_modules[i])
    sub_data = pd.concat([cond_row,sub_data])
    df_list.append(sub_data)

out2_df = pd.concat(df_list, axis=1)
out2_df.insert(0, 'Gene Modules', ["NA"]+out_df["Module"].values.tolist())
out2_df.to_csv(sys.argv[5], sep="\t")

#print condition clusters
for c in cond_modules:
    print(c)
    print("----------------------------------------------")
    for i in cond_modules[c]:
        print(i)
    print("==============================================")

