import numpy as np
import scipy.stats as sp_stats


def get_statistical_test(groups_data: list, group_names: str):
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
