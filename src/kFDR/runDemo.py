import kFDR
from kFDR.runkFDR import runkFDR
import pandas as pd
import os
#---------------------------------------------------------------------------#
#---------------------------------------------------------------------------#
#---------------------------------------------------------------------------#

def runDemoT():
    data = pd.read_excel(os.path.join(os.path.dirname(kFDR.__file__), 'data/demo_data1'))
    out1 = runkFDR(data = data, nLabels=2, nX1=4, nX2=4, data_table=True)
    return out1

def runDemoF():
    data = pd.read_excel(os.path.join(os.path.dirname(kFDR.__file__), 'data/demo_data2'))
    out2 = runkFDR(data = data, nLabels=2, nX1=4, nX2=4, data_table=False)
    return out2
