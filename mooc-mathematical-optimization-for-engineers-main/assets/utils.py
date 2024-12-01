##
#  @file itils.py
#
#  @brief Utils for training of artificial neural network in Keras and export to file that is readable by MeLOn.
#
# ==============================================================================\n
#   Aachener Verfahrenstechnik-Systemverfahrenstechnik, RWTH Aachen University  \n
# ==============================================================================\n
#
#  @author Artur M. Schweidtmann, Friedrich von Bülow, Jing Cui, Laurens Lueg, and Alexander Mitsos
#  @date 20. January 2020
#  
#  @author Johannes M. M. Faust, Ashutosh Manchanda
#  @date 11. Septemeber 2021
#  @change: reduce to used functions only
##

import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

def scale(X, scaleInput):
    # scale Input values to range [-1,1] in each dimension
    if (scaleInput):
        nom = (X -  X.min(axis=0))*2
        denom = X.max(axis=0) - X.min(axis=0)
        denom[denom==0] = 1
        return -1 + nom/denom
    else:
        return X

def normalize(y, normalizeOutput):
    # normalize output to z-score
    if(normalizeOutput):
        y_norm = (y - np.mean(y, axis=0))/np.std(y, axis=0);
        return y_norm
    else:
        return y
