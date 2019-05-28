# -*- coding: utf-8 -*-
"""
Created on Sun Oct 29 15:55:19 2017

@author: ciullo
"""
import random
import numpy as np

# werklijn function: step-wise distribution of high discharges.
def werklijn_cdf(Xlist, A):

    X = np.asarray(Xlist)
    # number of connecting points of piece-wise relation
    nl = np.shape(A)[0]
    a = A['a'].values
    b = A['b'].values
    A['Q'].loc[nl + 1] = np.inf
    XL = A['Q'].values     # limits of piece-wise relation

    # derive P-values
    P = np.repeat(np.nan, np.size(X))
    for j in range(0, nl):
        indexlow = X >= XL[j]
        indexup = X < XL[j + 1]
        index = np.where((indexlow * indexup) == True)[0]  # AND condition
        P[index] = np.exp(-np.exp(-(X[index] - b[j]) / a[j]))
    return P


def werklijn_inv(Plist, A):
    # inverse probability distribution function
    # probability is translated to frequency.
    # X is a piece-wise linear function of log(frequency)

    # input
    #   P:    probability of non-exceedance
    #   A:  parameters of the werklijn
    #
    # output
    #   X:    x-value, asociated with P

    # read input parameters
    P = np.asarray(Plist)
    # number of connecting points of piece-wise relation
    nl = np.shape(A)[0]
    a = A['a'].values
    b = A['b'].values
    # this is where you actually extrapolate way beyond the third connecting
    # point
    A['RP'].loc[nl + 1] = np.inf
    RPL = A['RP'].values    # limits of piece-wise relation
                            # transform probability of non-exceedance
                            # frequency of exceedance
    Fe = -np.log(P)
    RP = 1 / Fe    # return period

# derive X-values
    X = np.repeat(np.nan, np.size(P))
    for j in range(0, nl):
        indexlow = RP >= RPL[j]
        indexup = RP < RPL[j + 1]
        index = np.where((indexlow * indexup) == True)[0]  # AND condition
        X[index] = a[j] * np.log(RP[index]) + b[j]

    return X

def werklijn_pdf(Xlist, A):
    # pdf according to "werklijn"
    # probability is translated to frequency.
    # X is a piece-wise linear function of log(frequency)

    # input
    #   X:    x-value
    #   A:  parameters of the werklijn
    #
    # output
    #   P:    probability density

    # read input parameters
    X = np.array(Xlist)

    # number of connecting points of piece-wise relation
    nl = np.shape(A)[0]
    a = A['a'].values
    b = A['b'].values
    A['Q'].loc[nl + 1] = np.inf
    XL = A['Q'].values     # limits of piece-wise relation

    # derive P-values
    P = np.repeat(np.nan, np.size(X))
    for j in range(0, nl):
        indexlow = X >= XL[j]
        indexup = X < XL[j + 1]
        index = np.where((indexlow * indexup) == True)[0]  # AND condition
        P[index] = werklijn_cdf(X[index], A) * \
            np.exp(-(X[index] - b[j]) / a[j]) * (1 / a[j])
    return P

# randomly sample from werklijn
def rand_werklijn(A):  # compute x = inverseFx(Xu)
    u = random.random()     # compute Xu
    return werklijn_inv([u], A)  # return inverseFx(Xu)





















