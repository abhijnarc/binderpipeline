import pandas as pd
import numpy as np

# Load input file
df = pd.read_csv("bench20_ifmetric.csv")

#%%

# calculate additional parameters

# DG normalized by dSASA
df['dG_norm'] = np.where(
    df['dSASA'] > 0, 
    (df['dG'] / df['dSASA'])*100,
    0)

# Unsaturation density
df['unsat_density'] = np.where(
    df['dSASA'] > 0, 
    (df['unsatHB'] / df['dSASA'])*1000,
    0)

#%%
# agregate and calculate mean values for each complex

# mean of 3rd column onwards

cols_to_mean = df.columns[2:]

# now calculate mean score and entropy for each complex

df_complex= df.groupby(df.columns[0])[cols_to_mean].mean().reset_index()

#%%

# Save output
df_complex.to_csv("ifmetric_summary.csv", index=False)

print(df_complex)
