# Merge pilot feature tables and calculate pairwise-complete Spearman/Kendall
# correlations against experimental log10(Kd in M).

import os
from functools import reduce

import numpy as np
import pandas as pd
from scipy.stats import kendalltau, spearmanr

PREFIX = os.environ.get("PREFIX", "pilot10")
EXPT_FILE = os.environ.get("EXPT_FILE", "expt_kd.csv")
INTF_FILE = os.environ.get("INTF_FILE", f"{PREFIX}_interaction_scores_entropy.csv")
IFM_FILE = os.environ.get("IFM_FILE", f"{PREFIX}_ifmetric_summary.csv")
PROD_FILE = os.environ.get("PROD_FILE", f"{PREFIX}_prodigy.csv")

expt_df = pd.read_csv(EXPT_FILE)
intf_df = pd.read_csv(INTF_FILE)
ifm_df = pd.read_csv(IFM_FILE)
prod_df = pd.read_csv(PROD_FILE)

if len(expt_df.columns) < 2 or expt_df.columns[1] != "logKd":
    raise ValueError("Experimental CSV must have columns 'complex,logKd' in that order")
if "logKd" in prod_df.columns:
    prod_df = prod_df.rename(columns={"logKd": "prodigy_logKd"})

dataframes = [expt_df, intf_df, ifm_df, prod_df]
master_df = reduce(
    lambda left, right: pd.merge(left, right, on="complex", how="inner"),
    dataframes,
)
master_df["logKd"] = pd.to_numeric(master_df["logKd"], errors="coerce")

rows = []
for column in master_df.columns:
    if column in ("complex", "logKd"):
        continue
    values = pd.to_numeric(master_df[column], errors="coerce")
    paired = pd.DataFrame({"experimental": master_df["logKd"], "feature": values}).dropna()
    paired = paired[np.isfinite(paired.to_numpy()).all(axis=1)]
    n = len(paired)
    if n >= 2:
        rho, rho_p = spearmanr(paired["experimental"], paired["feature"])
        tau, tau_p = kendalltau(paired["experimental"], paired["feature"])
    else:
        rho = rho_p = tau = tau_p = float("nan")
    rows.append({
        "parameter": column,
        "n": n,
        "spearman_rho": rho,
        "spearman_pvalue": rho_p,
        "kendall_tau": tau,
        "kendall_pvalue": tau_p,
        "passes_abs_rho_0.3": bool(pd.notna(rho) and abs(rho) >= 0.3),
    })

correlations = pd.DataFrame(rows)
master_path = f"{PREFIX}_allparam.csv"
corr_path = f"{PREFIX}_correlations.csv"
master_df.to_csv(master_path, index=False)
correlations.to_csv(corr_path, index=False)
print(correlations.to_string(index=False))
print(f"Wrote {master_path} and {corr_path}")
