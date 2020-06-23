import numpy as np
#---------------------------------------#
#---------------------------------------#
#---------------------------------------#

def guessLog(X1, X2):
    X1 = X1
    X2 = X2
    X12 = np.column_stack((X1,X2))
    max_X12 = np.amax(X12)
    if max_X12 < 20:
        X1 = 2**X1
        X2 = 2**X2
    else:
        X1 = X1
        X2 = X2
    if max_X12 < 20:
        print ('Guess: log2 values identified :)!')
    else: print ('Guess: Data is not log transformed :)!')
    return X1, X2
