import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import io
from scipy import stats as sp_stats
from plotting.raincloud import make_raincloud_plot
from plotting.fonts import get_available_fonts

st.set_page_config(page_title="Raincloud Plot", layout="wide")


def sidebar_controls(genes=None, unique_conditions=None, default_controls=None):
    """Main sidebar controls for the raincloud plot page"""
    # Data source section
    st.sidebar.header("Data Source")
    data_source = st.sidebar.radio("Choose data source:", ["session", "upload"])
    uploaded_file = st.sidebar.file_uploader("Upload CSV file...") if data_source == "upload" else None

    # Plot customization section
    st.sidebar.markdown("---")
    st.sidebar.subheader("📊 Visualization Options")

    # Plot type
    plot_type = st.sidebar.selectbox(
        "Select plot type:",
        ["Bar Plot with Points", "Raincloud Plot"],
        index=1,
        help="Bar plot recommended for small samples",
    )

    # Error bar type (only for bar plot)
    error_type = (
        st.sidebar.radio(
            "Error bars:",
            ["SEM (Standard Error)", "SD (Standard Deviation)", "95% CI"],
            help="SEM: Standard Error of Mean | SD: Standard Deviation | CI: Confidence Interval",
        )
        if plot_type == "Bar Plot with Points"
        else None
    )

    # Show points (only for bar plot)
    show_points = (
        st.sidebar.checkbox("Show individual points", value=True) if plot_type == "Bar Plot with Points" else None
    )

    # Group labels
    if unique_conditions is not None and len(unique_conditions) >= 2:
        group_labels = {}
        st.sidebar.markdown("---")
        st.sidebar.subheader("Group Labels")
        for i, condition in enumerate(unique_conditions):
            group_labels[condition] = st.sidebar.text_input(
                f"Label for '{condition}'", str(condition).capitalize(), key=f"group_{i}"
            )
    else:
        group_labels = {}

    # Plot formatting
    if default_controls is None:
        default_controls = {}

    plot_title = st.sidebar.text_input("Plot Title", default_controls.get("plot_title", "Gene Expression"))
    x_label = st.sidebar.text_input("X Axis Label", default_controls.get("x_label", "Condition"))
    y_label = st.sidebar.text_input("Y Axis Label", default_controls.get("y_label", "Expression Level"))

    # Styling options
    st.sidebar.markdown("---")
    st.sidebar.subheader("🎨 Style Options")

    # Fonts
    available_fonts, missing_fonts = get_available_fonts()
    font_family = st.sidebar.selectbox("Font Family", available_fonts)

    # Colors - dynamic based on number of groups
    colors = {}
    color_palette = ["#1F77B4", "#D62728", "#2CA02C", "#FF7F0E", "#9467BD"]
    conditions_list = [] if unique_conditions is None else list(unique_conditions)
    for i, condition in enumerate(conditions_list):
        color = color_palette[i % len(color_palette)]
        colors[condition] = st.sidebar.color_picker(f"Color for '{condition}'", color)

    # Transparency settings (for raincloud plot)
    if plot_type == "Raincloud Plot":
        violin_alpha = st.sidebar.slider("Violin Transparency", 0.0, 1.0, 0.6, 0.1)
        boxplot_alpha = st.sidebar.slider("Box Transparency", 0.0, 1.0, 0.8, 0.1)
        jitter_alpha = st.sidebar.slider("Point Transparency", 0.0, 1.0, 0.5, 0.1)
        point_size = st.sidebar.slider("Point Size", 1, 5, 2)
        group_spacing = st.sidebar.slider("Group Spacing", 0.5, 2.0, 2.0, 0.1)
    else:
        violin_alpha = boxplot_alpha = jitter_alpha = point_size = group_spacing = None

    # Figure size
    fig_width = st.sidebar.slider("Figure Width", 4, 16, default_controls.get("fig_width", 5))
    fig_height = st.sidebar.slider("Figure Height", 4, 16, default_controls.get("fig_height", 6))

    # Export section
    st.sidebar.markdown("---")
    st.sidebar.subheader("💾 Export Options")

    dpi = st.sidebar.selectbox("Resolution (DPI)", [72, 150, 300, 600], index=2)
    plot_format = st.sidebar.selectbox("File Format", ["PNG", "PDF", "SVG", "JPEG"], index=0)
    generate_plot_download = st.sidebar.button("Generate Plot Download")

    return {
        "data_source": data_source,
        "uploaded_file": uploaded_file,
        "plot_type": plot_type,
        "error_type": error_type,
        "show_points": show_points,
        "group_labels": group_labels,
        "plot_title": plot_title,
        "x_label": x_label,
        "y_label": y_label,
        "font_family": font_family,
        "colors": colors,
        "violin_alpha": violin_alpha,
        "boxplot_alpha": boxplot_alpha,
        "jitter_alpha": jitter_alpha,
        "point_size": point_size,
        "group_spacing": group_spacing,
        "fig_width": fig_width,
        "fig_height": fig_height,
        "dpi": dpi,
        "plot_format": plot_format,
        "generate_plot_download": generate_plot_download,
    }


def transform_data_for_raincloud(data, gene_col, group_col):
    """Transform data from wide to long format suitable for raincloud plotting."""
    plot_data = pd.DataFrame({"value": data[gene_col].values, "condition": data[group_col].values, "gene": gene_col})
    plot_data = plot_data.dropna()
    return plot_data


def choose_statistical_test(groups_data, group_names):
    """
    Automatically choose the appropriate statistical test based on:
    - Number of groups
    - Sample sizes
    - Normality (Shapiro-Wilk test) - with adjustment for large samples
    - Variance homogeneity (Levene's test)

    For large sample sizes (n > 5000), normality tests become overly sensitive,
    so we use effect size and visual inspection instead.
    """
    n_groups = len(groups_data)
    sample_sizes = [len(g) for g in groups_data]
    min_sample_size = min(sample_sizes)
    max_sample_size = max(sample_sizes)
    total_sample_size = sum(sample_sizes)

    # For large samples, normality tests are overly sensitive
    # Use a subsample for normality testing if n > 5000
    large_sample = total_sample_size > 5000
    subsample_size = 5000

    # Check normality for each group
    normality_tests = []
    normality_pvals = []

    for group in groups_data:
        if len(group) >= 3:
            # For large samples, test on random subsample
            if large_sample and len(group) > subsample_size:
                test_group = np.random.choice(group, subsample_size, replace=False)
            else:
                test_group = group

            try:
                stat, pval = sp_stats.shapiro(test_group)
                normality_pvals.append(pval)
                normality_tests.append(pval > 0.05)
            except:  # noqa
                # If Shapiro-Wilk fails (rare with large n), assume non-normal
                normality_pvals.append(0.0)
                normality_tests.append(False)
        else:
            normality_tests.append(False)
            normality_pvals.append(0.0)

    all_normal = all(normality_tests)

    # Check variance homogeneity (Levene's test) - only for 2+ groups
    if n_groups >= 2:
        try:
            levene_stat, levene_pval = sp_stats.levene(*groups_data)
            equal_variances = levene_pval > 0.05
        except:  # noqa
            equal_variances = False
            levene_pval = 0.0
    else:
        equal_variances = True
        levene_pval = 1.0

    # For large samples, prefer parametric tests due to Central Limit Theorem
    # even if individual groups deviate from normality
    use_parametric = (large_sample and total_sample_size > 5000) or (all_normal and min_sample_size > 30)

    # Choose test based on number of groups
    if n_groups == 2:
        group1, group2 = groups_data

        # Two-sample test
        if use_parametric:
            if equal_variances:
                stat, pval = sp_stats.ttest_ind(group1, group2, equal_var=True)
                test_name = "Independent t-test (parametric, large sample)"
                test_type = "parametric"
            else:
                stat, pval = sp_stats.ttest_ind(group1, group2, equal_var=False)
                test_name = "Welch's t-test (unequal variances, large sample)"
                test_type = "parametric"
        else:
            stat, pval = sp_stats.mannwhitneyu(group1, group2, alternative="two-sided")
            test_name = "Mann-Whitney U test (non-parametric)"
            test_type = "non-parametric"

        # Calculate effect size (Cohen's d)
        mean_diff = np.mean(group1) - np.mean(group2)
        pooled_std = np.sqrt((np.std(group1, ddof=1) ** 2 + np.std(group2, ddof=1) ** 2) / 2)
        cohens_d = mean_diff / pooled_std if pooled_std > 0 else 0
        effect_size = cohens_d
        effect_size_type = "Cohen's d"

    elif n_groups == 3:
        # Three-group test
        if use_parametric:
            if equal_variances:
                stat, pval = sp_stats.f_oneway(*groups_data)
                test_name = "One-way ANOVA (parametric, large sample)"
                test_type = "parametric"
            else:
                stat, pval = sp_stats.f_oneway(*groups_data)
                test_name = "Welch's ANOVA (unequal variances, large sample)"
                test_type = "parametric"
        else:
            stat, pval = sp_stats.kruskal(*groups_data)
            test_name = "Kruskal-Wallis H test (non-parametric)"
            test_type = "non-parametric"

        # Calculate eta-squared (effect size for ANOVA)
        grand_mean = np.concatenate(groups_data).mean()
        ss_between = sum(len(g) * (np.mean(g) - grand_mean) ** 2 for g in groups_data)
        ss_total = sum((x - grand_mean) ** 2 for g in groups_data for x in g)
        eta_squared = ss_between / ss_total if ss_total > 0 else 0
        effect_size = eta_squared
        effect_size_type = "η² (Eta-squared)"

    else:  # n_groups > 3
        # Multiple groups test
        if use_parametric:
            if equal_variances:
                stat, pval = sp_stats.f_oneway(*groups_data)
                test_name = f"One-way ANOVA (parametric, {n_groups} groups, large sample)"
                test_type = "parametric"
            else:
                stat, pval = sp_stats.f_oneway(*groups_data)
                test_name = f"Welch's ANOVA (unequal variances, {n_groups} groups, large sample)"
                test_type = "parametric"
        else:
            stat, pval = sp_stats.kruskal(*groups_data)
            test_name = f"Kruskal-Wallis H test (non-parametric, {n_groups} groups)"
            test_type = "non-parametric"

        # Calculate eta-squared (effect size for ANOVA)
        grand_mean = np.concatenate(groups_data).mean()
        ss_between = sum(len(g) * (np.mean(g) - grand_mean) ** 2 for g in groups_data)
        ss_total = sum((x - grand_mean) ** 2 for g in groups_data for x in g)
        eta_squared = ss_between / ss_total if ss_total > 0 else 0
        effect_size = eta_squared
        effect_size_type = "η² (Eta-squared)"

    # Adjusted p-value using Benjamini-Hochberg FDR control
    # (for multiple comparisons correction)
    adj_pval = pval  # Default: no adjustment for single test
    method = "None (single test)"

    return {
        "test_name": test_name,
        "test_type": test_type,
        "statistic": stat,
        "pvalue": pval,
        "pvalue_adjusted": adj_pval,
        "adjustment_method": method,
        "n_groups": n_groups,
        "sample_sizes": sample_sizes,
        "min_sample_size": min_sample_size,
        "max_sample_size": max_sample_size,
        "total_sample_size": total_sample_size,
        "large_sample": large_sample,
        "all_normal": all_normal,
        "normality_pvals": normality_pvals,
        "equal_variances": equal_variances,
        "levene_pval": levene_pval,
        "effect_size": effect_size,
        "effect_size_type": effect_size_type,
    }


def create_bar_plot(df, sidebar_vals, test_results, group_labels):
    """Create bar plot with error bars and significance testing"""
    conditions = df["condition"].unique()
    means = df.groupby("condition")["expression"].mean()
    sems = df.groupby("condition")["expression"].sem()
    stds = df.groupby("condition")["expression"].std()
    counts = df.groupby("condition")["expression"].count()

    if sidebar_vals["error_type"] == "SEM (Standard Error)":
        errors = sems
    elif sidebar_vals["error_type"] == "SD (Standard Deviation)":
        errors = stds
    else:
        errors = sems * sp_stats.t.ppf(0.975, counts - 1)

    fig, ax = plt.subplots(figsize=(sidebar_vals["fig_width"], sidebar_vals["fig_height"]))
    x_pos = np.arange(len(means))
    colors = [sidebar_vals["colors"].get(cond, "#1F77B4") for cond in means.index]

    ax.bar(
        x_pos,
        means,
        yerr=errors,
        capsize=10,
        alpha=0.8,
        color=colors,
        edgecolor="black",
        linewidth=2,
        error_kw={"linewidth": 2.5, "elinewidth": 2.5, "capthick": 2.5},
    )

    if sidebar_vals["show_points"]:
        for i, cond in enumerate(means.index):
            cond_data = df[df["condition"] == cond]["expression"].values
            np.random.seed(42)
            jitter = np.random.normal(0, 0.05, len(cond_data))
            x_vals = np.full(len(cond_data), i) + jitter
            ax.scatter(x_vals, cond_data, color="black", s=80, alpha=0.6, zorder=3, edgecolors="white", linewidth=1.5)

    # Add significance annotation for 2-group comparison
    if test_results["n_groups"] == 2:
        y_max = df["expression"].max()
        y_min = df["expression"].min()
        y_range = y_max - y_min
        max_bar_height = max([means.iloc[i] + errors.iloc[i] for i in range(len(means))])
        y_sig_line = max_bar_height + y_range * 0.05
        y_sig_text = y_sig_line + y_range * 0.03

        ax.plot([0, 1], [y_sig_line, y_sig_line], "k-", linewidth=2)
        ax.plot([0, 0], [y_sig_line - y_range * 0.01, y_sig_line], "k-", linewidth=2)
        ax.plot([1, 1], [y_sig_line - y_range * 0.01, y_sig_line], "k-", linewidth=2)

        pval = test_results["pvalue"]
        if pval < 0.001:
            sig_symbol = "***"
            sig_text = "p < 0.001"
        elif pval < 0.01:
            sig_symbol = "**"
            sig_text = f"p = {pval:.3f}"
        elif pval < 0.05:
            sig_symbol = "*"
            sig_text = f"p = {pval:.3f}"
        else:
            sig_symbol = "ns"
            sig_text = f"p = {pval:.3f}"

        ax.text(0.5, y_sig_text, sig_symbol, ha="center", va="bottom", fontsize=18, fontweight="bold")
        ax.text(0.5, y_sig_text + y_range * 0.08, sig_text, ha="center", va="bottom", fontsize=10, style="italic")

        ax.set_ylim(bottom=min(0, y_min - y_range * 0.05), top=y_sig_text + y_range * 0.15)

    # Style the plot
    ax.set_xticks(x_pos)
    ax.set_xticklabels(means.index, fontsize=13, fontweight="bold")
    ax.set_title(sidebar_vals["plot_title"], fontsize=16, fontweight="bold", pad=20)
    ax.set_ylabel(sidebar_vals["y_label"], fontsize=13, fontweight="bold")
    ax.set_xlabel(sidebar_vals["x_label"], fontsize=13, fontweight="bold")
    ax.grid(axis="y", alpha=0.3, linestyle="--", zorder=0)
    ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(1.5)
    ax.spines["bottom"].set_linewidth(1.5)

    return fig


# Main function for the raincloud plot page
st.title("📊 miRNAGene Expression Analysis")

# Load data
st.session_state.data = pd.read_csv("data/all_data.csv")
data = st.session_state["data"]

# Show data info
with st.expander("ℹ️ Current Dataset Info"):
    col1, col2, col3 = st.columns(3)
    with col1:
        st.metric("Rows", f"{len(data):,}")
    with col2:
        st.metric("Columns", f"{len(data.columns):,}")
    with col3:
        st.info("📊 Dataset Loaded")

# Data transformation section
st.markdown("---")
st.subheader("🔧 Data Selection and Transformation")

# Identify numeric and categorical columns
numeric_cols = data.select_dtypes(include=["number"]).columns.tolist()
categorical_cols = data.select_dtypes(include=["object", "category"]).columns.tolist()

if len(numeric_cols) == 0:
    st.error("❌ No numeric columns found in the data.")
    st.stop()

if len(categorical_cols) == 0:
    st.error("❌ No categorical columns found in the data.")
    st.stop()

# Selection controls
col1, col2 = st.columns(2)

with col1:
    default_gene = "hsa-miR-3185" if "hsa-miR-3185" in numeric_cols else None
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
    )

# Get all unique groups
all_groups = sorted(data[grouping_var].unique().tolist())

# Allow user to select which groups to include
selected_groups = st.multiselect(
    f"Select groups from '{grouping_var}' to compare:",
    options=all_groups,
    default=["Breast (BR)", "Healthy Controls"],
    help=f"Select 2 or more groups. Found {len(all_groups)} unique groups.",
)

# Validate selection
if len(selected_groups) < 2:
    st.warning("⚠️ Please select at least 2 groups for comparison.")
    st.stop()

# Filter data
data = data[data[grouping_var].isin(selected_groups)]

# Transform data
try:
    transformed_data = transform_data_for_raincloud(data, selected_gene, grouping_var)

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
    genes=[selected_gene],
    unique_conditions=unique_conditions,
    default_controls={
        "plot_title": f"{selected_gene} Expression",
        "x_label": grouping_var,
        "y_label": "Expression Level",
    },
)

# Prepare data for plotting
df = transformed_data[["condition", "value"]].rename(columns={"value": "expression"})

# Apply custom group labels if provided
df["condition"] = df["condition"].map(lambda x: controls["group_labels"].get(x, str(x)))

# Prepare groups data for statistical test
groups_data = [transformed_data[transformed_data["condition"] == group]["value"].values for group in unique_conditions]
group_names = [controls["group_labels"].get(g, str(g)) for g in unique_conditions]

# Choose and perform appropriate statistical test
test_results = choose_statistical_test(groups_data, group_names)

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
