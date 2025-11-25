import pandas as pd

data = pd.read_csv("all_data.csv")
data.target.value_counts().reset_index().to_csv("label_distribution.csv", index=False)
