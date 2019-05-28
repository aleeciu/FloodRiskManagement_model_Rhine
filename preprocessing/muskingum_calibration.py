# -*- coding: utf-8 -*-
"""
Created on Mon Nov 20 13:16:06 2017

@author: ciullo
"""
import pandas as pd
import scipy as sp
import numpy as np
import matplotlib.pyplot as plt
from functions.funs_dikes import Lookuplin

def coefficients(deltaT, K, X):
    den = (2 * K * (1 - X) + deltaT)

    C1 = (deltaT - 2 * K * X) / den
    C2 = (deltaT + 2 * K * X) / den
    C3 = (2 * K * (1 - X) - deltaT) / den

    C = (C1, C2, C3)

    return C


def Muskingum(params, upstream, deltaT=1):
    ' It simulates Muskingum between two locations, used for multiple '
    '  Muskingum over the stretch                                     '
    k, X = params
    C1, C2, C3 = coefficients(deltaT, k, X)

    Qsim = np.repeat(upstream[0], len(upstream))

    for t in range(1, len(upstream)):
        Qsim[t] = upstream[t] * C1 + upstream[t - 1] * C2 + Qsim[t - 1] * C3

    return Qsim


data_net = pd.read_excel('./data/Rhine_network.xlsx', 
                         dtype=object, index_col=0)

''' BENCHMARKING '''

files = list('0ABCDE')  # BovnRijn, pnKan, Waal, Lek, IJssel

branches = [91, 101, 201, 301, 401, 501]

link = dict(zip(branches, files))

k0 = 1
X0 = 0.47
init_guesses = np.array([k0, X0])
params = pd.DataFrame()

#plt.figure(figsize = (15,15))
for b in branches:
    fl = link[b]
    data = pd.read_excel(
        './data/hydraulics/SOBEK_discharges/{}_branch_Q.xlsx'.format(fl)).iloc[:, 1:]
    # from hours to days, neglect last number == rest of division
    # .to_dict('list')
    data = data.drop(
        data.index[-1]).groupby((data.index[:-1] / 24).astype(int)).mean()

    data.columns = data_net['NodeName'][data_net['BRANCH'] == b]

#        # Only use most downstream for a two location calibration
#        data = data.iloc[:, [0, -1]]

    params_branch = pd.DataFrame(index=data.columns[:-1], columns=['K','X'])
    # from upstream to penultime location
    for col in data.columns[:-1]:

        ind = data.columns.values.tolist().index(col)
        upstream = data[col]
        downstream = data[data.columns[ind + 1]]

        def errfunc(prms): return Muskingum(
            prms, upstream=upstream) - downstream

        result = sp.optimize.least_squares(
            errfunc, init_guesses, bounds=([0, 0], [np.inf, 0.5])).x
#                result = sp.optimize.least_squares(errfunc, init_guesses).x

        sim = Muskingum(result, upstream=upstream)

        K, X = result

        params_branch.loc[col, 'K'], params_branch.loc[col, 'X'] = K, X
        params_branch.loc[col, 'C1'] = (1 - 2 * K * X) / (2 * K * (1 - X) + 1)
        params_branch.loc[col, 'C2'] = (1 + 2 * K * X) / (2 * K * (1 - X) + 1)
        params_branch.loc[col, 'C3'] = (
            2 * K * (1 - X) - 1) / (2 * K * (1 - X) + 1)

        plt.plot(sim, '-.', label='Muskingum')
        plt.plot(downstream, label='SOBEK')
#       plt.plot(upstream, label = 'upstream')
        plt.title(col)
        plt.legend()
        plt.savefig('./Muskingmum_fit_pics/{}.png'.format(col))
        plt.show()

    params = pd.concat([params, params_branch], axis=0)

#params.to_excel('./data/hydraulics/musk_calib_factor.xlsx')

#''' TEST the found factors on a wave different than the one used in calibration'''
#from functools import partial
#params = pd.read_excel('./data/hydraulics/musk_calib_factor.xlsx', dtype = object)
#
### TEST DATA:
#discharges = pd.read_excel('./data/hydraulics/SOBEK_discharges/4checking/hydrographs.xlsx').iloc[:, 1:]
#discharges = discharges.drop(discharges.index[-1]).groupby((discharges.index[:-1] / 24).astype(int)).mean()
#
### CALIBRATION DATA:
#discharges = pd.DataFrame()
#for fl in ['A', 'C', 'D', 'E']:
#if fl == 'A':
##               data = pd.read_excel('./data/hydraulics/SOBEK_discharges/{}_branch.xlsx'.format(fl)).iloc[:, [1]]
#else:
##               data = pd.read_excel('./data/hydraulics/SOBEK_discharges/{}_branch.xlsx'.format(fl)).iloc[:, [-1]]
##
##        data = data.drop(data.index[-1]).groupby((data.index[:-1] / 24).astype(int)).mean()
##        discharges = pd.concat([discharges, data], axis = 1)
#
## Name of upstream location per branch
#discharges.columns = ['101.860', '201.960', '401.977', '501.977']
#Q_distr = pd.read_excel('./data/hydraulics/distr_factors.xlsx')
#nodes_to_skip = data_net['NodeName'][(data_net['type'] == 'downstream') | (data_net['type'] == 'bifurcation')]
#
#seqs = [(101, 201), (101, 301, 401), (101, 301, 501)]
#
# for seq in seqs[0:1]:
#        upstream = discharges['101.860']
#
##        f, s, t = seq
#
#        f, s = seq
#
#        t = 20000
#
#        new_data = data_net[(data_net['BRANCH'] == f) | (data_net['BRANCH'] == s) | (data_net['BRANCH'] == t)]
#
#        for node_name in new_data['NodeName']:
#
#                if node_name not in nodes_to_skip.values:
#
#                        index = params['NodeName'] == node_name
#
#                        prms = params.loc[index, ['K', 'X']].values[0]
#                        downstream = Muskingum(prms, upstream, deltaT = 1)
#                        upstream = downstream
#
#                        if node_name == '{}.0'.format(s):
#
#                                 cols = ['{}.1000_Q'.format(f), '{}_ratio'.format(node_name)]
#
#                                 upstream *= map(partial(Lookuplin, Q_distr[cols].values, 0, 1),
#                                                 upstream)
#
#                        elif node_name == '{}.0'.format(t):
#
#                                 cols = ['{}.1000_Q'.format(s), '{}_ratio'.format(node_name)]
#
#                                 upstream *= map(partial(Lookuplin, Q_distr[cols].values, 0, 1),
#                                                 upstream)
#
#                if node_name in ['201.960', '401.977', '501.977']:
#                        plt.plot(range(len(discharges[node_name])), discharges[node_name])
#                        plt.plot(range(len(downstream)), downstream, '--')
#                        plt.title('{}_{}'.format(seq[0], seq[-1]))
#        plt.show()
