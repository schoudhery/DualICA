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
#C_matrix.index = split_index.get_level_values(0)

lfcs_data = pd.read_csv(sys.argv[3], sep="\t", index_col=0,low_memory=False)
categories_data = pd.read_csv(sys.argv[4], sep=",")
categories_data = categories_data.loc[categories_data["Name"].isin(lfcs_data.columns), :]
categories_data =categories_data.sort_values(by=['Classification'])

lfcs_data= lfcs_data[categories_data["Name"]]

lfcs_data = lfcs_data.loc[:, ~lfcs_data.columns.isin(["sentinel", "nTA"])]
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


out_list = []
y_out=[]
X_out_list=[]
lfcs_data = lfcs_data.reset_index()
melted_lfcs = pd.melt(lfcs_data, id_vars="index")
melted_lfcs.columns = ["gene","condition","LFC"]
for idx,row in melted_lfcs.iterrows():
    g = row["gene"]
    c= row["condition"]
    lfc = row["LFC"]
    interaction_terms = []
    IC_gene = G_matrix.loc[g,:]
    IC_cond = C_matrix.loc[c,:]
    interaction_terms = np.outer(IC_gene, IC_cond, out=None)
    flattened_interactions=[]
    for i in interaction_terms:
         flattened_interactions.extend(i) 
    res_list = [lfc,g,c]+flattened_interactions
    out_list.append(res_list)
    y_out.append(lfc)
    X_out_list.append(flattened_interactions)

print("Fitting Linear Regression")
X= sm.add_constant(X_out_list)
model = sm.OLS(y_out, X)
results = model.fit(method="qr",use_t=True)
pvals = results.pvalues
adjusted_pvals = statsmodels.stats.multitest.fdrcorrection(pvals,alpha=0.05)[1]
coeff = results.params
zscores = scistat.zscore(coeff)

results_df = pd.DataFrame({"pvals":pvals,"adj_pvals":adjusted_pvals,"coeff":coeff, "z-score":zscores})
index=["cnst"]
for i in G_matrix: 
    for j in C_matrix:
        index.append(str(i)+"_"+str(j))
results_df.index=index
results_df.to_csv("./linreg_results.csv")

results_df = results_df.iloc[1:,:]
results_df["Gene Set"] = results_df.index.str.split("_").str[0]
results_df["Gene Set"] = results_df["Gene Set"].astype(int)
results_df["Condition Set"] = results_df.index.str.split("_").str[1]
results_df["Condition Set"] = results_df["Condition Set"].astype(int)

coeff_matrix = results_df.pivot(index='Gene Set', columns='Condition Set', values='coeff')
coeff_matrix=coeff_matrix.sort_index()
coeff_matrix=coeff_matrix.sort_index(axis=1)
coeff_matrix.to_csv("./linear_coeffs.tsv", sep="\t")