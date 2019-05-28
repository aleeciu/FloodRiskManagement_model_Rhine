# -*- coding: utf-8 -*-
"""
Created on Fri Mar 16 19:37:55 2018

@author: ciullo
"""
from copy import deepcopy
import numpy as np
import pandas as pd
from scipy.stats import norm
from functools import partial
import matplotlib.pyplot as plt
from functions.funs_hydrostat import werklijn_pdf, rand_werklijn, werklijn_inv
import networkx as nx
from dike_model_function import DikeNetwork

''' Generate water levels '''

# Monte Carlo approach
nsample = 3000
A = pd.read_excel('./data/hydraulics/werklijn_params.xlsx')

load = True                           # load files of Q, strenghts, w levels
load_Q = False
IS = True                             # important sampling
calibration = True                     # calibration or testing
policy = False                         # new policy or target failure probability

if load:
    samples = pd.read_excel(
    './data/fragility_curves/calibration/Qpeaks_{}.xlsx'.format(nsample),
    index_col=0).T.values[0]
#    random_numbers = pd.read_excel('./data/fragility_curves/calibration/strengths_{}.xlsx'.format(nsample), index_col=0).values

#     outcomes =  pd.read_excel('./data/fragility_curves/calibration/maxWL_{}_divpolicy.xlsx'.format(nsample)).to_dict()
    outcomes = pd.read_excel(
    './data/fragility_curves/calibration/maxWL_{}.xlsx'.format(nsample), index_col=0)

    corr_f = werklijn_pdf(samples, A) / (1 / float(np.max(samples) - np.min(samples)))
    # Prepare input data:
    dike_nodes = DikeNetwork(samples).dikenodes
    G = DikeNetwork(samples).G

else:
    # if you only want to generate new water levels:
    if load_Q:
        samples = pd.read_excel(
        './data/fragility_curves/calibration/Qpeaks_{}.xlsx'.format(nsample)).T.values[0]  # .values

    else:
        # important sampling:
        if IS:
            # lowest VNK pf is 81, VNK discharges have been simulated from one year RT
            # up to events with a RT of 12500. IS helps speeding up convergence!
            low, high = werklijn_inv([1 - 1 / 1.001, 1 - 1 / 12500], A)
            samples = np.random.uniform(low, high, nsample)
            corr_f = werklijn_pdf(samples, A) / (1 / float(high - low))
        else:
            # normal sampling:
            samples = [rand_werklijn(A) for _ in range(0, nsample)]

        # Save samples:
#        pd.DataFrame(samples).to_excel(
#            './data/fragility_curves/calibration/Qpeaks_{}.xlsx'.format(nsample))
    # Save strengths:
    random_numbers = np.random.uniform(size = nsample)
#    pd.DataFrame(random_numbers).to_excel('./data/fragility_curves/calibration/strengths_{}.xlsx'.format(nsample))

    inputs = {}
    factors = {'Bmax': 200, 'pfail': 0.5, 'Brate': 0.5, 'DikeIncrease':0}

    diversion = {'Diversion_401.0': -2.0, 'Diversion_201.0': 0.0}

    # Prepare input data - !modify dikenetwork in order for it to get Qpeak in input!
    network = DikeNetwork(samples)
    dike_nodes = network.dikenodes
    G = network.G

    for node in dike_nodes:
        inputs.update({'{}_{}'.format(node, fac)
                      : factors[fac] for fac in factors})

    inputs.update(diversion)

    # Compute high wate levels - specify them as model output in the DikeNet function
    outcomes = network.__call__(**inputs)
    new_outcome_keys = [key for key in outcomes.keys() if key.split('_')[0] == 'wl']
    outcomes = {key: outcomes[key] for key in new_outcome_keys}
    
#    pd.DataFrame(outcomes).to_excel('./data/fragility_curves/calibration/maxWL_{}_divpolicy.xlsx'.format(nsample))
#    pd.DataFrame(outcomes).to_excel(
#        './data/fragility_curves/calibration/maxWL_{}_{}_{}.xlsx'.format(nsample, int(low), int(high)))

pfailure = {n: 0 for n in dike_nodes}

if calibration:
#    outcomes = pd.read_excel(
#    './data/fragility_curves/calibration/maxWL_{}.xlsx'.format(nsample)).to_dict()
    # initialize the adjstment factors matrix
    guesses_df = pd.DataFrame(np.zeros((len(dike_nodes), 3)), index=dike_nodes,
                              columns=['low_guess', 'mean', 'up_guess'])

    for node in dike_nodes:
        guesses_df.loc[node]['low_guess'] = -100
        guesses_df.loc[node]['up_guess'] = 100

    pf = {'different': 'VNK', 'equal': 1250.0}
    pf_type = 'different'

    # target probability of failure
    pf_target = {}
    for node in dike_nodes:
        if pf_type == 'different':
            pf_target.update({node: 1 / round(G.node[node]['VNK_pf'], 1)})

        elif pf_type == 'equal':
            pf_target.update({node: 1 / 1250.0})

    #pf_target = nx.get_node_attributes(G, 'pfailure')
    print('Calibrate Fragility Curves')
    
    _pfailure = np.asarray(list(pfailure.values()))
    _ptarget = np.asarray(list(pf_target.values()))
    
    while max(abs(guesses_df['low_guess'] - guesses_df['up_guess'])) > 0.0001 and any(
            abs(_pfailure - _ptarget) / _ptarget > 0.0001):

        guesses_df['mean'] = guesses_df[['low_guess', 'up_guess']].mean(axis=1)
        frag_curve = deepcopy(nx.get_node_attributes(G, 'f'))

        for node in dike_nodes:
            p_failure = []
            frag_curve[node]['mu'] += guesses_df['mean'].loc[node]

            strength = norm.ppf(
                0.5,
                frag_curve[node]['mu'],
                frag_curve[node]['sd'])

#            o = 'wl_{}'.format(node)
            o = '{}_maxwl'.format(node)
            z = outcomes.loc[:,o] > np.repeat(strength, nsample)

            if IS:
                #   p_failure.append(np.sum(z*corr_f)/float(nsample))
                p_failure = np.sum(z * corr_f) / float(nsample)
            else:
                #                   p_failure.append(np.sum(z)/float(nsample))
                p_failure = np.sum(z) / float(nsample)

#                    pfailure[node] = np.mean(p_failure)
            pfailure[node] = p_failure

            if abs(pfailure[node] - pf_target[node])/pf_target[node] < 0.0001:
                   print('{}_converged'.format(node))
            print(node, pfailure[node])

            if pfailure[node] > pf_target[node]:  # here put a tolerance level
                guesses_df['low_guess'].loc[node] = guesses_df['mean'].loc[node]

            elif pfailure[node] < pf_target[node]:  # here put a tolerance level
                guesses_df['up_guess'].loc[node] = guesses_df['mean'].loc[node]

        _pfailure = np.asarray(list(pfailure.values()))
        _ptarget = np.asarray(list(pf_target.values()))

    correct = pd.DataFrame(nx.get_node_attributes(
        G, 'VNK_pf'), index=['correct']).T.dropna()
    estimated = pd.DataFrame({key: round(1 / pfailure[key]) for key in correct.index},
                             index=['estimated']).T
    p_estimation = pd.concat([correct, estimated], axis=1)

    guesses_df['NodeName'] = guesses_df.index
    guesses_df.index = range(len(dike_nodes))
#    guesses_df[['mean', 'NodeName']].to_excel(
#        './data/fragility_curves/calibration/FC_calfactors_{}_{}_IS{}.xlsx'.format(pf[pf_type], nsample, IS))

elif policy:
    # mind that it works only for dike raising where you just need to shift the
    # fragility curve, with rfr you may need to run the whole model

    mean_FC = pd.read_excel(
        './data/fragility_curves/calibration/FC_calfactors_VNK_{}_IS{}.xlsx'.format(
            nsample, IS), dtype=object)
    mean_FC = mean_FC.set_index('NodeName')

    pol = pd.read_excel(
        './results/RWS/results_optimization_RWS_100000nfe.xlsx').iloc[:, :-3]
    for dh in pol.columns:
        dikecode = dh.split('_')[0]
        mean_FC.loc[dikecode] += pol[dh].values / 10.0

    frag_curve = deepcopy(nx.get_node_attributes(G, 'f'))

    for node in dike_nodes:
        p_failure = []
        frag_curve[node]['mu'] += mean_FC.loc[node].values[0]

        strength = norm.ppf(0.5, frag_curve[node]['mu'],
                            frag_curve[node]['sd'])

        o = '{}_maxwl'.format(node)
        z = outcomes[o].values() > np.repeat(strength, nsample)

        if IS:
            p_failure.append(np.sum(z * corr_f) / float(nsample))
        else:
            p_failure.append(np.sum(z) / float(nsample))

        pfailure[node] = np.mean(p_failure)

    correct = pd.DataFrame(nx.get_node_attributes(
        G, 'VNK_pf'), index=['correct']).T.dropna()

    estimated = pd.DataFrame({key: '{:f}'.format(1.0 / pfailure[key]) for key in pfailure.keys()},
                             index=['estimated']).T

    p_estimation = pd.concat([correct, estimated], axis=1)
    p_estimation['dike_increase'] = pol.T.values / 10.0
    p_estimation
