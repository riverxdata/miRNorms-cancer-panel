from plotnine import (
    ggplot,
    aes,
    geom_violin,
    geom_boxplot,
    geom_jitter,
    theme_minimal,
    labs,
    scale_fill_manual,
    theme,
    element_text,
    element_line,
    element_blank,
    position_jitter,
)


def make_raincloud_plot(df, controls, gene, pvalue, test_name):
    return (
        ggplot(df, aes(x="condition", y="expression", fill="condition"))
        + geom_violin(
            aes(x="condition"),
            alpha=controls["violin_alpha"],
            width=0.4 * controls["group_spacing"],
            show_legend=False,
            style="left",
        )
        + geom_boxplot(
            aes(x="condition"),
            width=0.15 * controls["group_spacing"],
            alpha=controls["boxplot_alpha"],
            show_legend=False,
            outlier_alpha=0,
        )
        + geom_jitter(
            alpha=controls["jitter_alpha"],
            size=controls["point_size"],
            show_legend=False,
            position=position_jitter(width=0.05 * controls["group_spacing"], height=0),
        )
        + labs(title=f"{gene} expression ({test_name}, p = {pvalue:.4f})", x="Condition", y="Expression Level")
        + theme_minimal()
        + theme(
            figure_size=(controls["fig_width"], controls["fig_height"]),
            plot_title=element_text(size=14, weight="bold", family=controls["font_family"]),
            axis_title_x=element_text(size=12, family=controls["font_family"]),
            axis_title_y=element_text(size=12, family=controls["font_family"]),
            axis_text_x=element_text(size=10, family=controls["font_family"]),
            axis_text_y=element_text(size=10, family=controls["font_family"]),
            panel_border=element_blank(),
            axis_line=element_line(color="black", size=1),
            panel_grid_major=element_blank(),
            panel_grid_minor=element_blank(),
            axis_ticks_major_y=element_line(color="black", size=0.5),
            axis_ticks_length_major=8,
        )
    )
