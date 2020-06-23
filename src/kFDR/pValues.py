from scipy import stats
#---------------------------------------------------------------------------#
#Student's t-test assumes that the two populations have normal distributions with equal variances.
#Welch's t-test is designed for unequal variances, but the assumption of normality is maintained.
#Welch's t-test is an approximate solution to the Behrens–Fisher problem.
# https://en.wikipedia.org/wiki/Welch%27s_t-test
# https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.ttest_ind.html
#---------------------------------------------------------------------------#
def t_test(X1, X2):
    X1 = X1
    X2 = X2
    t_value, p_value = stats.ttest_ind(X1, X2, axis=1, equal_var=False, nan_policy='omit')
    return t_value, p_value
