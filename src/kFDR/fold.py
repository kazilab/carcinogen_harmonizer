import numpy as np

#---------------------------------------------------------------------------#
#calculate fold #
#---------------------------------------------------------------------------#

def cal_fold(X1, X2):
    X1 = X1
    X2 = X2
    X1_av = np.average(X1, axis=1, weights=None, returned=False)
    X2_av = np.average(X2, axis=1, weights=None, returned=False)
    pre_fold = X2_av/X1_av
    for i, x in enumerate(pre_fold): # enumerate: mention (a number of things) one by one.
        if x >= 1:
            pre_fold[i] = x
        else:
            pre_fold[i] = -1/x
    return X1_av, X2_av, pre_fold
