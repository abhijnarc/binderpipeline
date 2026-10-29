# -*- coding: utf-8 -*-
"""
Created on Mon Oct 12 16:12:35 2026

merge all dataframes and calculate correlation with expt data

@author: sujan
"""

import pandas as pd
from functools import reduce
from scipy.stats import spearmanr, kendalltau

#%%

# file names
expt_file = 'expt_kd.csv'
intf_file = 'interaction_scores_entropy.csv'
ifm_file = 'ifmetric_summary.csv'
prod_file = 'bench20_prod.csv'

# get the data
expt_df = pd.read_csv(expt_file)
intf_df = pd.read_csv(intf_file)
ifm_df = pd.read_csv(ifm_file)
prod_df = pd.read_csv(prod_file)

#%%

# nake a list of dataframes
df_list = [expt_df, intf_df, ifm_df, prod_df]

master_df = reduce(lambda left, right: 
                   pd.merge(left, right, on = 'complex', how = 'inner'),
                   df_list)

#%%    
# Now we calculate the correlations
expt_logKd = master_df.iloc[:, 1]

other_cols = master_df.iloc[:, 2:]

spearman_corr = other_cols.corrwith(expt_logKd, method = 'spearman')


print(spearman_corr)

master_df.to_csv("bench20_allparam.csv", encoding = 'utf-8', index = False)
