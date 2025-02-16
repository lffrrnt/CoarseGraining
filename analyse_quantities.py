import numpy as np
from scipy import stats


def silence_proba(data):
    """
    Calculate the average probability of silence at each coarse
    graining step.  
    --------------------------------------
    Inputs: Coarse-grained variables matrix in shape (N, T) where 
    T is the length of the time series, and N is the number 
    of new variables obtained at step k.
    --------------------------------------
    Return: Real value.
    """
    P = np.zeros(len(data))
    for i in range(len(data)):
        P[i] = (data.shape[1] - np.count_nonzero(data[i,:])) / data.shape[1]
    return np.mean(P), stats.sem(P)


def compute_variance(data):
    """
    Calculate the variance of all clusters at each coarse graining 
    step k.
    --------------------------------------
    Inputs: Non-normalized coarse-grained variables matrix in shape (T, N)
    where T is the length of the time series, and N is the number 
    of new variables obtained at step k.
    --------------------------------------
    Return: Variance array of all steps.
    """
    var = np.zeros(len(data))
    for i in range(len(data)):
        var[i] = np.var(data[i])

    return var


def acrl(x,t):
   xs = x - np.mean(x)
   nrm = np.sqrt(np.mean(xs**2)*np.mean(xs**2))
   tot = len(x)
   if nrm != 0:
       autocorr = np.array([np.mean(xs[i:]*xs[:tot-i])/nrm for i in t])
   else:
       autocorr = np.zeros(len(t))
   return autocorr

def find_auto_correlation_variables(data, t):
    auto_corr = {}
    for idx, i in enumerate(data):
        auto_corr[idx] = np.zeros((i.shape[0],t))
        for j in range(i.shape[0]):
            auto_corr[idx][j,:] = acrl(i[j,:], range(t))
    return auto_corr


    

























