#import pandas as pd
import kFDR
from kFDR.runFDR import runFDR


#-----------------------------------------------#
#-----------------------------------------------#
#-----------------------------------------------#

def runkFDR(data, nLabels, nX1, nX2, data_table):
    data = data
    nLabels = nLabels
    nX1 = nX1
    nX2 = nX2
    data_table = data_table
    out = runFDR(data, nLabels, nX1, nX2, data_table)
    return out
