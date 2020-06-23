# import data from excel file
import pandas as pd

#----------------------------------------------#
#----------------------------------------------#
#----------------------------------------------#

def data_in(data, nLabels, nX1, nX2):
    data = data
    nLabels = nLabels
    nX1 = nX1+nLabels
    nX2 = nX2+nX1
    labels = data.iloc[:,0:nLabels]
    X1 = data.iloc[:,nLabels:nX1].values
    X2 = data.iloc[:,nX1:nX2].values
    return labels, X1, X2
