import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
import io
from plotting.raincloud import make_raincloud_plot
from pages.utils.data import transform_data_for_raincloud
from pages.utils.statistics import get_statistical_test
from pages.utils.mirna_ui import sidebar_controls, create_bar_plot

st.set_page_config(page_title="Raincloud Plot", layout="wide")

# Main function for the raincloud plot page
st.title("📊 miRNAGene Expression Analysis")

# Show data info
col1, col2, col3 = st.columns(3)
with col1:
    st.metric("Rows", "16,910 samples")
with col2:
    st.metric("Columns", "2,554 miRNAs")
with col3:
    st.metric("Metadata", "2 groups")

# Data transformation section
st.markdown("---")
st.subheader("🔧 Data Selection and Transformation")

# Identify numeric and categorical columns
numeric_cols = [col for col in pd.read_csv("data/all_data.csv", nrows=0).columns.tolist() if col.startswith("hsa-")]
categorical_cols = [
    col for col in pd.read_csv("data/all_data.csv", nrows=0).columns.tolist() if not col.startswith("hsa-")
]
# Selection controls
col1, col2 = st.columns(2)

with col1:
    default_gene = "hsa-miR-1246"
    selected_gene = st.selectbox(
        "Select gene/feature to visualize:",
        numeric_cols,
        index=numeric_cols.index(default_gene) if default_gene in numeric_cols else 0,
        help="Choose the numeric column representing expression values",
    )

    with col2:
        default_grouping = "label"
        grouping_var = st.selectbox(
            "Select grouping variable:",
            categorical_cols,
            index=categorical_cols.index(default_grouping) if default_grouping in categorical_cols else 0,
            help="Variable to group and compare samples",
            key="grouping_var_select",
        )

# Load data
st.session_state.data = pd.read_csv("data/all_data.csv", usecols=[grouping_var, selected_gene])
data = st.session_state["data"]

# Get all unique groups - force reload when grouping_var changes
if "prev_grouping_var" not in st.session_state:
    st.session_state.prev_grouping_var = None

# Check if grouping variable has changed
if st.session_state.prev_grouping_var != grouping_var:
    st.session_state.prev_grouping_var = grouping_var
    # Clear the multiselect key to force reload
    if f"multiselect_{grouping_var}" in st.session_state:
        del st.session_state[f"multiselect_{grouping_var}"]

all_groups = sorted(data[grouping_var].unique().tolist())

# Get default selection (first 2 groups or previous selection)
selected_groups = st.multiselect(
    f"Select groups from '{grouping_var}' to compare:",
    options=all_groups,
    default=(
        ["Breast (BR)", "Healthy Controls"]
        if "Breast (BR)" in all_groups and "Healthy Controls" in all_groups
        else all_groups
    ),
    help=f"Select 2 or more groups. Found {len(all_groups)} unique groups.",
    key=f"multiselect_{grouping_var}",
)


# Filter data
if len(selected_groups) < 2:
    st.error("❌ Please select at least 2 groups to compare.")
    st.stop()
else:
    viz_data = data[[grouping_var, selected_gene]]
    viz_data = viz_data[viz_data[grouping_var].isin(selected_groups)]
    viz_data[grouping_var] = pd.Categorical(viz_data[grouping_var], categories=selected_groups, ordered=True)
    viz_data = viz_data.sort_values(grouping_var)
    # transform data
    try:
        transformed_data = transform_data_for_raincloud(viz_data, selected_gene, grouping_var)

        with st.expander("📊 Preview Transformed Data"):
            st.markdown(f"**Data shape:** {transformed_data.shape[0]} rows × {transformed_data.shape[1]} columns")

            group_summary = (
                transformed_data.groupby("condition")["value"]
                .agg([("Count", "count"), ("Mean", "mean"), ("Median", "median"), ("Std", "std")])
                .round(4)
            )

            st.markdown(f"**Summary by {grouping_var}:**")
            st.dataframe(group_summary, width="stretch")

            st.markdown("**First 10 rows:**")
            st.dataframe(transformed_data.head(10), width="stretch", hide_index=True)

        st.success(f"✅ Data ready! {len(selected_groups)} groups selected")

    except Exception as e:
        st.error(f"❌ Error: {str(e)}")
        st.exception(e)
        st.stop()

st.markdown("---")
# Get controls from sidebar
unique_conditions = transformed_data["condition"].unique()
controls = sidebar_controls(
    unique_conditions=unique_conditions,
    default_controls={
        "plot_title": f"{selected_gene} Expression",
        "x_label": grouping_var,
        "y_label": "Expression Level",
    },
)

# Prepare data for plotting
df = transformed_data[["condition", "value"]].rename(columns={"value": "expression"})
df["condition"] = df["condition"].map(lambda x: controls["group_labels"].get(x, str(x)))
groups_data = [transformed_data[transformed_data["condition"] == group]["value"].values for group in unique_conditions]
group_names = [controls["group_labels"].get(g, str(g)) for g in unique_conditions]
# Choose and perform appropriate statistical test
test_results = get_statistical_test(groups_data, group_names)
# Create plot based on selection
if controls["plot_type"] == "Bar Plot with Points":
    fig = create_bar_plot(df, controls, test_results, controls["group_labels"])
else:
    plot = make_raincloud_plot(df, controls, selected_gene, test_results["pvalue"], test_results["test_name"])
    fig = plot.draw()

# Display plot
col1, col2, col3 = st.columns([1, 3, 1])
with col2:
    st.pyplot(fig)
# Statistical analysis results
st.subheader("📊 Statistical Analysis Results")
col1, col2, col3, col4 = st.columns(4)
with col1:
    st.metric("Test", test_results["test_name"].split("(")[0].strip())
with col2:
    st.metric("Statistic", f"{test_results['statistic']:.4f}")
with col3:
    st.metric("P-value", f"{test_results['pvalue']:.4e}")
with col4:
    pval = test_results["pvalue"]
    if pval < 0.001:
        sig = "*** p < 0.001"
    elif pval < 0.01:
        sig = "** p < 0.01"
    elif pval < 0.05:
        sig = "* p < 0.05"
    else:
        sig = "ns (p ≥ 0.05)"
    st.metric("Significance", sig)

# Effect size metrics
col1, col2, col3 = st.columns(3)
with col1:
    st.metric(f"Effect Size ({test_results['effect_size_type']})", f"{test_results['effect_size']:.4f}")
with col2:
    st.metric("Sample Type", "Large sample (n > 5000)" if test_results["large_sample"] else "Standard")
with col3:
    st.metric("Total Sample Size", f"{test_results['total_sample_size']:,}")

# Detailed test information
with st.expander("ℹ️ Test Selection Details"):
    col1, col2, col3 = st.columns(3)

    with col1:
        st.markdown("**Test Characteristics:**")
        st.write(f"- Number of groups: {test_results['n_groups']}")
        st.write(f"- Min sample size: {test_results['min_sample_size']:,}")
        st.write(f"- Max sample size: {test_results['max_sample_size']:,}")
        st.write(f"- Total samples: {test_results['total_sample_size']:,}")
        st.write(f"- Test type: {test_results['test_type'].upper()}")

    with col2:
        st.markdown("**Assumptions Check:**")
        st.write(f"- All groups normal: {'✅ Yes' if test_results['all_normal'] else '❌ No'}")
        st.write(f"- Equal variances: {'✅ Yes' if test_results['equal_variances'] else '❌ No'}")
        st.write(f"- Large sample mode: {'✅ Yes' if test_results['large_sample'] else '❌ No'}")
        if test_results["n_groups"] >= 2:
            st.write(f"- Levene's test p-value: {test_results['levene_pval']:.4f}")

    with col3:
        st.markdown("**Normality Tests (Shapiro-Wilk):**")
        for i, (group_name, pval) in enumerate(zip(group_names, test_results["normality_pvals"])):
            status = "✅ Normal" if pval > 0.05 else "❌ Non-normal"
            st.write(f"- {group_name}: p = {pval:.4f} {status}")

    st.markdown("---")
    st.markdown(
        f"""
    **Test Selection Logic:**

    **Sample Size:** {test_results['total_sample_size']:,} samples

    **For Large Samples (n > 5000):**
    - Normality tests are overly sensitive, so Central Limit Theorem applies
    - Parametric tests (t-test/ANOVA) are preferred due to robustness
    - Levene's test used for variance equality

    **Test Decision Tree:**
    - **2 Groups:** Welch's t-test (unequal variance assumption) → Mann-Whitney U (if non-parametric)
    - **3+ Groups:** Welch's ANOVA (unequal variance assumption) → Kruskal-Wallis (if non-parametric)

    **Effect Size Interpretation:**
    - {test_results['effect_size_type']}: {test_results['effect_size']:.4f}
    - {"Small" if abs(test_results['effect_size']) < 0.2 else "Medium" if abs(test_results['effect_size']) < 0.5 else "Large"} effect
    """
    )

# Summary statistics
with st.expander("📈 Summary Statistics by Group"):
    summary_stats = (
        df.groupby("condition")["expression"]
        .agg(
            [
                "count",
                "mean",
                "median",
                "std",
                "min",
                "max",
                ("Q1", lambda x: x.quantile(0.25)),
                ("Q3", lambda x: x.quantile(0.75)),
            ]
        )
        .round(4)
    )
    st.dataframe(summary_stats, width="stretch")

# Raw data
with st.expander("📋 View Raw Data"):
    st.dataframe(df, width="stretch", hide_index=True)

# Download options
csv_data = df.to_csv(index=False).encode("utf-8")
st.sidebar.download_button(
    "📄 Download Data (CSV)", data=csv_data, file_name=f"{selected_gene}_expression_data.csv", mime="text/csv"
)

# Plot download
if controls["generate_plot_download"]:
    buf = io.BytesIO()
    fig.savefig(buf, format=controls["plot_format"].lower(), dpi=controls["dpi"], bbox_inches="tight")
    buf.seek(0)
    mimes = {"PNG": "image/png", "PDF": "application/pdf", "SVG": "image/svg+xml", "JPEG": "image/jpeg"}
    st.sidebar.download_button(
        f"📥 Download {controls['plot_format']}",
        data=buf,
        file_name=f"{selected_gene}_{controls['plot_type'].replace(' ', '_').lower()}.{controls['plot_format'].lower()}",
        mime=mimes[controls["plot_format"]],
    )
    st.sidebar.success(f"✅ Plot ready at {controls['dpi']} DPI")
    plt.close(fig)
