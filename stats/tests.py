from scipy.stats import ttest_ind, mannwhitneyu, shapiro

def test_normality(group):
    p = shapiro(group)[1]
    return p

def choose_test(group1, group2):
    p1, p2 = test_normality(group1), test_normality(group2)
    if p1 > 0.05 and p2 > 0.05:
        stat, pval = ttest_ind(group1, group2)
        return "t-test", stat, pval, p1, p2
    else:
        stat, pval = mannwhitneyu(group1, group2)
        return "Wilcoxon rank-sum", stat, pval, p1, p2