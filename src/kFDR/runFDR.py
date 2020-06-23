#import pandas as pd
import kFDR
from kFDR.data_in import data_in
from kFDR.smartGuess import guessLog
from kFDR.pValues import t_test
from kFDR.fold import cal_fold
from kFDR.join import t_test_fold, t_test_fold_plus
from kFDR.calFDR import fdr_values


import pandas as pd
#-----------------------------------------------#
#-----------------------------------------------#
#-----------------------------------------------#

def runFDR(data, nLabels, nX1, nX2, data_table):
    data = data
    nLabels = nLabels
    nX1 = nX1
    nX2 = nX2
    data_table = data_table
    labels, X1, X2 = data_in(data=data, nLabels=nLabels, nX1=nX1, nX2=nX2)
    X1 , X2 = guessLog(X1, X2)
    t_value, p_value = t_test(X1, X2)
    X1_av, X2_av, fold = cal_fold(X1, X2)
    if data_table == False:
        labels_t_p_value_fold = t_test_fold(t_value, p_value, X1_av, X2_av, fold, labels)
    else:
        labels_t_p_value_fold = t_test_fold_plus(t_value, p_value, X1_av, X2_av, fold, data)
    allFDR = fdr_values(p_value)
    out=pd.concat([labels_t_p_value_fold, allFDR], axis=1, sort=False)
    return out
