import streamlit as st
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import ttest_ind, mannwhitneyu, shapiro

# Simulate RNAseq data
np.random.seed(42)
n_samples = 30

data = pd.DataFrame(
    {
        "sample": [f"S{i+1}" for i in range(n_samples * 2)],
        "condition": ["cancer"] * n_samples + ["normal"] * n_samples,
        "BRCA1": np.concatenate(
            [
                np.random.normal(loc=10, scale=2, size=n_samples),  # cancer
                np.random.normal(loc=8, scale=2, size=n_samples),  # normal
            ]
        ),
        "BRCA2": np.concatenate(
            [
                np.random.normal(loc=7, scale=1.5, size=n_samples),  # cancer
                np.random.normal(loc=6.5, scale=1.2, size=n_samples),  # normal
            ]
        ),
    }
)

# UI
st.title("🎻 Boxplot Visualizer for Gene Expression")

gene = st.selectbox("Select Gene", ["BRCA1", "BRCA2"])

# Prepare data
df = data[["condition", gene]].rename(columns={gene: "expression"})

# Subset by groups
group1 = df[df["condition"] == "cancer"]["expression"]
group2 = df[df["condition"] == "normal"]["expression"]

# Normality test (Shapiro)
p1 = shapiro(group1)[1]
p2 = shapiro(group2)[1]

# Choose test
if p1 > 0.05 and p2 > 0.05:
    test_name = "t-test"
    stat, p_value = ttest_ind(group1, group2)
else:
    test_name = "Wilcoxon rank-sum"
    stat, p_value = mannwhitneyu(group1, group2)

# Plot
fig, ax = plt.subplots()
sns.boxplot(data=df, x="condition", y="expression", ax=ax, palette="Set2")
sns.stripplot(
    data=df, x="condition", y="expression", ax=ax, color="black", jitter=True, size=5
)

# Display p-value
ax.set_title(f"{gene} expression ({test_name}, p = {p_value:.4f})")

# Render
st.pyplot(fig)
