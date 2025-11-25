import streamlit as st
import pandas as pd
import plotly.express as px

summary_label = pd.read_csv("data/label_distribution.csv")
summary_category = pd.read_csv("data/category_distribution.csv")
st.set_page_config(
    page_title="RIVER-VIS - miRNA Cancer Atlas", page_icon="🧬", layout="wide", initial_sidebar_state="expanded"
)

st.title("🧬 RIVER-VIS")
st.subheader("Rapid Interactive Visualization Environment for miRNA Expression Analysis")
st.markdown("## 📊 Dataset Summary")
metric_col1, metric_col2, metric_col3, metric_col4 = st.columns(4)

with metric_col1:
    st.metric("Total Samples", "16,190")
with metric_col2:
    st.metric("Cancer Types", "13")
with metric_col3:
    st.metric("Study ID", "GSE211692")
with metric_col4:
    st.metric("Sample Type", "Serum miRNA")


st.markdown(
    """
Welcome to **RIVER-VIS** - your comprehensive toolkit for visualizing and analyzing miRNA expression 
data across multiple cancer types.

---

## 🎯 Project Overview

This application provides rapid visualization and statistical analysis of serum miRNA profiles from 
**GSE211692**: a large-scale study of 9,921 patients with 13 types of human solid cancers. 
Built with Streamlit, it enables interactive exploration, group comparisons, and generation of 
publication-ready visualizations.

## 📊 Dataset Information
"""
)


viz_col1, viz_col2 = st.columns(2)


with viz_col1:
    st.markdown("#### 🎯 Cancer Types Distribution")
    fig_label = px.bar(
        summary_label,
        x="label",
        y="count",
        color="label",
        labels={"label": "Cancer Type", "count": "Number of Samples"},
        color_discrete_sequence=px.colors.qualitative.Set3,
    )
    fig_label.update_layout(
        showlegend=False,
        xaxis_tickangle=-45,
        plot_bgcolor="white",
        font=dict(size=11),
        xaxis=dict(title_font=dict(weight="bold")),
        yaxis=dict(title_font=dict(weight="bold")),
    )
    st.plotly_chart(fig_label, width="stretch")

with viz_col2:
    st.markdown("#### 📋 Sample Categories Distribution")
    fig_category = px.bar(
        summary_category,
        x="category",
        y="count",
        color="category",
        labels={"category": "Sample Category", "count": "Number of Samples"},
        color_discrete_sequence=px.colors.qualitative.Pastel,
    )
    fig_category.update_layout(
        showlegend=False,
        xaxis_tickangle=-45,
        plot_bgcolor="white",
        font=dict(size=11),
        xaxis=dict(title_font=dict(weight="bold")),
        yaxis=dict(title_font=dict(weight="bold")),
    )
    st.plotly_chart(fig_category, width="stretch")


col1, col2, col3 = st.columns(3)

with col1:
    st.info(
        """
    ### 📊 Raincloud Plots
    
    Visualize miRNA distributions with:
    - Box plots & violin plots
    - Jittered data points
    - Statistical annotations
    - Multiple cancer comparisons
    """
    )

with col2:
    st.info(
        """
    ### 📈 Statistical Analysis
    
    Perform comprehensive tests:
    - T-tests & ANOVA
    - Non-parametric tests
    - Multiple testing corrections
    - Effect size calculations
    """
    )

with col3:
    st.info(
        """
    ### 🔬 Multi-Cancer Exploration
    
    Compare across cohorts:
    - 13 cancer types
    - Benign conditions
    - Healthy controls
    - miRNA biomarker discovery
    """
    )

st.markdown("---")

st.markdown(
    """
## 📖 How to Use

1. **Navigate** using the sidebar menu (👈)
2. **Select miRNAs** - Choose specific miRNAs for analysis
3. **Define Groups** - Filter by cancer type or condition
4. **Visualize** - Generate raincloud plots and distributions
5. **Analyze** - Perform statistical tests between groups
6. **Explore** - Generate heatmaps and comparative analyses

### 💡 Quick Tips
- Start with the **miRNA Analysis** page to explore expression patterns
- Use **Raincloud Plots** for beautiful distribution visualizations
- Compare multiple cancer types simultaneously
- Statistical significance is automatically assessed (α = 0.05)
- Identify potential miRNA biomarkers across cancer types

### 🔍 Research Applications
- **Cancer Biomarker Discovery** - Identify miRNAs that differentiate cancer types
- **Clinical Validation** - Compare cancer vs. benign vs. control samples
- **Comparative Analysis** - Explore expression patterns across 13 cancer types
- **Statistical Power** - Large sample sizes enable robust analysis

---

*RIVER-VIS is designed for cancer researchers, bioinformaticians, and data scientists working with 
miRNA expression data and cancer genomics.*
"""
)
