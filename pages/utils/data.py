import pandas as pd


def transform_data_for_raincloud(data: pd.DataFrame, gene_col: str, group_col: str) -> pd.DataFrame:
    """Transform data from wide to long format suitable for raincloud plotting."""
    plot_data = pd.DataFrame({"value": data[gene_col].values, "condition": data[group_col].values, "gene": gene_col})
    plot_data = plot_data.dropna()
    return plot_data
