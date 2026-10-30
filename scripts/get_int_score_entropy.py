import pandas as pd
import numpy as np
import os

PREFIX = os.environ.get("PREFIX", "pilot10")

# Load input file
df = pd.read_csv(f"{PREFIX}_intsummary.csv")

# Define interaction weights
weights = {
    "SaltBridge": 0.4,
    "PiStack": 0.3,
    "HydrogenBond": 0.2,
    "Hydrophobic": 0.1
}

# Calculate weighted interaction score
df["weighted_score"] = df["interaction"].map(weights) * df["count"]

# Group by complex and pose
group_cols = ["complex", "nconf"]

result = (
    df.groupby(group_cols, as_index=False)
      .agg(score=("weighted_score", "sum"))
)

# Calculate total interactions for each complex-pose
df["total_count"] = df.groupby(group_cols)["count"].transform("sum")

# Calculate fraction of each interaction type
df["p"] = df["count"] / df["total_count"].where(df["total_count"] > 0)

# Calculate -p * ln(p), treating p = 0 as zero
df["entropy_term"] = 0.0
positive = df["p"] > 0
df.loc[positive, "entropy_term"] = -df.loc[positive, "p"] * np.log(df.loc[positive, "p"])

# Sum entropy terms for each complex-pose
entropy = (
    df.groupby(group_cols, as_index=False)
      .agg(shannon_entropy=("entropy_term", "sum"))
)

# Merge score and entropy
result = result.merge(entropy, on=group_cols)

# now calculate mean score and entropy for each complex
result_complex = result.groupby('complex', as_index=False).agg(
            {'score': 'mean', 'shannon_entropy': 'mean'})

# Save output
result_complex.to_csv(f"{PREFIX}_interaction_scores_entropy.csv", index=False)

print(result_complex)
