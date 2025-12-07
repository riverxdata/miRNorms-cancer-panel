import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
from plotting.fonts import get_available_fonts
from scipy import stats as sp_stats


def sidebar_controls(unique_conditions=None, default_controls=None):
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
        # print(unique_conditions)
        group_labels = {}
        st.sidebar.markdown("---")
        st.sidebar.subheader("Group Labels")
        for condition in unique_conditions:
            group_labels[condition] = st.sidebar.text_input(
                f"Label for '{condition}'", str(condition).capitalize(), key=f"group_label_{condition}"
            )
        print(group_labels)
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
    available_fonts, _ = get_available_fonts()
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


def create_bar_plot(df, sidebar_vals, test_results, group_labels):
    """Create bar plot with error bars and significance testing"""
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
