import pandas as pd
import numpy as np
import os

PREFIX = os.environ.get("PREFIX", "pilot10")

# Load input file
df = pd.read_csv(f"{PREFIX}_ifmetric.csv")

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

# The two derived columns are added above, so columns [2:] includes both
# dG_norm and unsat_density as well as the original InterfaceAnalyzer metrics.

cols_to_mean = df.columns[2:]

# now calculate mean score and entropy for each complex

df_complex= df.groupby(df.columns[0])[cols_to_mean].mean().reset_index()

#%%

# Save output
df_complex.to_csv(f"{PREFIX}_ifmetric_summary.csv", index=False)

print(df_complex)
