Installation: (tested on Mac OS X only)

On Windows:
python setup.py install

On Linux:
sudo python setup.py install 

On Mac OS X:
sudo python ./setup.py install 


To run demo data

Jupyter Notebook script

import kFDR
Example1:
kFDR.runDemo.runDemoT()
Example2:
kFDR.runDemo.runDemoF()


To run kFDR:

data file containing at least one column with data labels (Gene symbol, etc.), two data groups with replicates in column.
Example:
--------------------------------------------------------------------
|Symbol|Gr-1 rep1|Gr-1 rep2|Gr-1 rep3|Gr-2 rep1|Gr-2 rep1|Gr-2 rep1|
--------------------------------------------------------------------
Jupiter Notebook script

import pandas as pd
from kFDR.runkFDR import runkFDR
data = pd.read_csv('csv data file including path') or
data = pd.read_excel('csv data file including path')
runkFDR(data=data, nLabels=2, nX1=4, nX2=4, data_table=True)

nLabels = number of columns with data labels,
nX1 = number of replicates for group 1 (columns),
nX2 = number of replicates for group 2 (columns)
data_table = True or False. True, if data table is included in the output otherwise False.

returns pandas data frame with following parameters:
Group 1 average	
Group 2 average	
t-value	
p-value	
fold	
Bonferroni	
Sidak	
Holm	
Holm-Sidak	
Simes-Hochberg	
Hommel	
FDR Benjamini-Hochberg	
FDR Benjamini-Yekutieli	
FDR 2-stage Benjamini-Hochberg	
FDR 2-stage Benjamini-Krieger-Yekutieli	
FDR adaptive Gavrilov-Benjamini-Sarkar
-------------------------------------
T.TEST: https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.ttest_ind.html
Adjusted p-values: https://www.statsmodels.org/dev/generated/statsmodels.stats.multitest.multipletests.html
------------------------------------