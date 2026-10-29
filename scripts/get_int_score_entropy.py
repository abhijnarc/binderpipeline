import pandas as pd
import numpy as np

# Load input file
df = pd.read_csv("bench20_intsummary.csv")

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
df["p"] = df["count"] / df["total_count"]

# Calculate -p * ln(p), treating p = 0 as zero
df["entropy_term"] = np.where(
    df["p"] > 0,
    -df["p"] * np.log(df["p"]),
    0
)

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
result_complex.to_csv("interaction_scores_entropy.csv", index=False)

print(result_complex)
