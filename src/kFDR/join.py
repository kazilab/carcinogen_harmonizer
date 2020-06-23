import pandas as pd

#---------------------------------------------------------------------------#
#---------------------------------------------------------------------------#
#---------------------------------------------------------------------------#

def t_test_fold(t_value, p_value, X1_av, X2_av, fold, labels):
    t_valueD = pd.DataFrame(t_value)
    t_valueD.columns = ['t-value']
    p_valueD = pd.DataFrame(p_value)
    p_valueD.columns = ['p-value']
    X1_avD = pd.DataFrame(X1_av)
    X1_avD.columns = ['Group 1 average']
    X2_avD = pd.DataFrame(X1_av)
    X2_avD.columns = ['Group 2 average']
    foldD = pd.DataFrame(fold)
    foldD.columns = ['fold']
    labelsD = pd.DataFrame(labels)
    labels_t_p_value_fold=pd.concat([labelsD, X1_avD, X2_avD, t_valueD, p_valueD, foldD], axis=1, sort=False)
    return labels_t_p_value_fold

def t_test_fold_plus(t_value, p_value, X1_av, X2_av, fold, data):
    t_valueD = pd.DataFrame(t_value)
    t_valueD.columns = ['t-value']
    p_valueD = pd.DataFrame(p_value)
    p_valueD.columns = ['p-value']
    X1_avD = pd.DataFrame(X1_av)
    X1_avD.columns = ['Group 1 average']
    X2_avD = pd.DataFrame(X1_av)
    X2_avD.columns = ['Group 2 average']
    foldD = pd.DataFrame(fold)
    foldD.columns = ['fold']
    dataF = pd.DataFrame(data)
    labels_t_p_value_fold_plus=pd.concat([dataF, X1_avD, X2_avD, t_valueD, p_valueD, foldD], axis=1, sort=False)
    return labels_t_p_value_fold_plus
