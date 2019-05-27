# -*- coding: utf-8 -*-
"""
Created on Fri Mar 16 19:37:55 2018

@author: ciullo
"""

from __future__ import division
from copy import deepcopy
from ema_workbench import ema_logging

from funs_generate_network import get_network
from funs_dikes import Lookuplin, dikefailure, init_node
import numpy as np
import pandas as pd
from scipy.stats import norm
from functools import partial
import matplotlib.pyplot as plt
from funs_hydrostat import werklijn_pdf, rand_werklijn, werklijn_inv, werklijn_cdf
import networkx as nx


def Muskingum(C1, C2, C3, Qn0_t1, Qn0_t0, Qn1_t0):
    Qn1_t1 = C1 * Qn0_t1 + C2 * Qn0_t0 + C3 * Qn1_t0

    return Qn1_t1


class DikeNetwork(object):
    def __init__(self):
        # load network
        network = get_network(Musk_params=True, fragcurves=False, Losses=True)

        G = network[0]
        branches = network[1]
        upstream_node = network[2]
        transboundary_nodes = network[3]
        national_nodes = network[4]

        self.G = G
        self.transboundary_nodes = transboundary_nodes
        self.national_nodes = national_nodes
        self.dikenodes = transboundary_nodes + national_nodes
        self.branches = branches
        self.upstream_node = upstream_node

        self.sb = False
        ema_logging.info('model initialized')

    def _initialize_hydroloads(self, node, time, Q_0):
        node['cumVol'], node['wl'], node['Qpol'], node['hbas'] = (
            init_node(0.0, time) for _ in range(4))
        node['Qin'], node['Qout'] = (init_node(Q_0, time) for _ in range(2))
        node['status'] = init_node(False, time)
        node['tbreach'] = np.nan
        return node

    def _initialize_rfr(self, G, dikenodes):
        # Initialize room for the river
        G.node['rfr']['rfr_cost'] = 0
        # Create a copy of the rating curve that will be used in the
        # simulation:
        for node in dikenodes:
            G.node[node]['rnew'] = deepcopy(G.node[node]['r'])
        return G

    def __call__(self, timestep=1, div_frac=None,
                 Qpeaks=[20000], q=4, **kwargs):

        G = self.G
        dikenodes = self.dikenodes
        upst_node = self.upstream_node

        self._initialize_rfr(G, dikenodes)

        # load all kwargs into network
        for item in kwargs:
            string1, string2 = item.split('_')
            G.node[string1][string2] = kwargs[item]

            # string1: usually dikename
            # string2: usually name of uncertainty or lever
            # exception only for rfr strategies where: 'project number id_rfr'

        # Dictionary storing outputs:
        data = {}
        # Outputs of interest:
#        Qpol = {dike: [] for dike in dikenodes}
#            Qout = {dike: [] for dike in dikenodes}
#        Qin = {dike: [] for dike in dikenodes}
        maxwl = {dike: [] for dike in dikenodes}
#        status = {dike: [] for dike in dikenodes}

        for Qp in Qpeaks:
            #           THE upstream event:
            time = np.arange(
                0, G.node[upst_node]['Qevents_shape'].loc[q].shape[0], timestep)
            G.node[upst_node]['Qin'] = Qp * \
                G.node[upst_node]['Qevents_shape'].loc[q]

            # Select a branch:
            for j in self.branches.keys():
                branchnodes = self.branches[j]

                # Initialize upstream hydrpgraph
                G.node[branchnodes[0]]['Qout'] = G.node[G.node[branchnodes[0]]
                                                        ['prec_node']]['Qin']
                Q_0 = int(G.node[branchnodes[0]]['Qout'][0])

                # Impose initial conditions at each node:
                for n in branchnodes[1:]:
                    node = G.node[n]

                    if node['type'] == 'dike':
                        self._initialize_hydroloads(node, time, Q_0)
                    # this could go before the loop, but we are not sure if pfail
                    # is iterated before or after dikeincrease:
                        node['critWL'] = norm.ppf(node['pfail'], node['f']['mu'],
                                                  node['f']['sd'])

                    elif node['type'] == 'bifurcation':

                        # Assign the diversion measure to the corresponding
                        # branches:
                        Qbif = '{}_Q'.format(n)
                        r1 = '{}_ratio'.format(node['subs_node1'])
                        r2 = '{}_ratio'.format(node['subs_node2'])

                        G.node[node['subs_node2']]['div_measure'] = pd.DataFrame(
                            {k: G.node['Diversion'][k].values()
                             for k in [Qbif, r2]},
                            columns=[r2, Qbif]).values

                        # Apply the diversion policy to the even branches:
                        G.node[node['subs_node2']]['div_measure'][:, 0] *= 1 + \
                            G.node['Diversion']['{}'.format(
                                node['subs_node2'])] / 10.0

                        G.node[node['subs_node1']]['div_measure'] = pd.DataFrame(
                            {k: G.node['Diversion'][k].values()
                             for k in [Qbif, r1]},
                            columns=[r1, Qbif]).values

                        G.node[node['subs_node1']]['div_measure'][:, 0] = 1 - \
                            G.node[node['subs_node2']]['div_measure'][:, 0]

                        # Qin for the three nodes of the bifurcation:
                        node['Qin'] = (init_node(Q_0, time))

                        G.node[node['subs_node1']]['Qin'] = map(partial(Lookuplin,
                                                                        G.node[node['subs_node1']]['div_measure'], 1, 0),
                                                                node['Qin']) * node['Qin']

                        G.node[node['subs_node2']]['Qin'] = map(partial(Lookuplin,
                                                                        G.node[node['subs_node2']]['div_measure'], 1, 0),
                                                                node['Qin']) * node['Qin']

                # Propagate the discharge wave:
                for t in range(1, len(time)):

                    # Run over each node of the branch:
                    for n in range(0, len(branchnodes)):
                        node = G.node[branchnodes[n]]

                        if node['type'] in ['upstream',
                                            'downstream', 'conjunction']:
                            continue

                        # Muskingum parameters:
                        C1 = node['C1']
                        C2 = node['C2']
                        C3 = node['C3']

                        prec_node = G.node[node['prec_node']]

                        node['Qin'][t] = Muskingum(C1, C2, C3,
                                                   prec_node['Qout'][t],
                                                   prec_node['Qout'][t - 1],
                                                   node['Qin'][t - 1])

                        if node['type'] == 'dike':

                            node['wl'][t] = Lookuplin(
                                node['rnew'], 0, 1, node['Qin'][t])

                            res = dikefailure(self.sb,
                                              node['Qin'][t], node['wl'][t],
                                              node['hbas'][t], node['hground'],
                                              node['status'][t -
                                                             1], node['Bmax'],
                                              node['Brate'], time[t],
                                              node['tbreach'], node['critWL'])

                            node['Qout'][t] = res[0]
                            node['Qpol'][t] = res[1]
                            node['status'][t] = res[2]
                            node['tbreach'] = res[3]

                            # TODO: discretization time, make it automatic:
                            node['cumVol'][t] = np.trapz(node['Qpol']) * 24

                            # the 1/125.0 area is considered the minimum one:
                            Area = Lookuplin(
                                node['table'], 6, 0, node['wl'][t])
                            node['hbas'][t] = node['cumVol'][t] / float(Area)

                        elif node['type'] == 'bifurcation':

                            # Distribute Qin into the two branches. Note:
                            # Qins and Qouts identical, no breach can happen.

                            G.node[node['subs_node1']]['Qin'][t] = Lookuplin(
                                G.node[node['subs_node1']]['div_measure'],
                                1, 0, node['Qin'][t]) * node['Qin'][t]

                            G.node[node['subs_node2']]['Qin'][t] = Lookuplin(
                                G.node[node['subs_node2']]['div_measure'],
                                1, 0, node['Qin'][t]) * node['Qin'][t]

                for n in branchnodes:
                    if G.node[n]['type'] == 'dike':
                        node = G.node[n]

#                        Qout[n].append(node['Qout'])
#                        Qin[n].append(node['Qin'])
#
#                        Qpol[n].append(node['Qpol'])
                        maxwl[n].append(np.max(node['wl']))
#                        status[n].append(node['status'][-1])

        for n in dikenodes:
            data.update({
                #                             '{}_Qout'.format(n): Qout[n],
                #                             '{}_Qin'.format(n): Qin[n],
                '{}_maxwl'.format(n): maxwl[n],
                #                             '{}_Qpol'.format(n): Qpol[n],
                #                             '{}_status'.format(n): status[n]
            })

        return data, G


''' Generate water levels '''

# Monte Carlo approach
nsample = 3000
A = pd.read_excel('./data/hydraulics/werklijn_params.xlsx')

load = True                            # load files of Q, strenghts, w levels
load_Q = False
IS = True                              # important sampling
calibration = True                    # calibration or testing
policy = False                          # new policy or target failure probability

if load:
    samples = pd.read_excel(
        './data/fragility_curves/calibration/Qpeaks_{}.xlsx'.format(nsample)).T.values[0]  # .values
#     random_numbers = pd.read_excel('./data/fragility_curves/calibration/strengths_{}.xlsx'.format(nsample, IS)).values

#     outcomes =  pd.read_excel('./data/fragility_curves/calibration/maxWL_{}_divpolicy.xlsx'.format(nsample)).to_dict()
    outcomes = pd.read_excel(
        './data/fragility_curves/calibration/maxWL_{}.xlsx'.format(nsample)).to_dict()

    corr_f = werklijn_pdf(samples, A) / \
        (1 / float(np.max(samples) - np.min(samples)))
    # Prepare input data:
    dike_nodes = DikeNetwork().dikenodes
    G = DikeNetwork().G

else:
    # if you only want to generate new water levels:
    if load_Q:
        samples = pd.read_excel(
            './data/fragility_curves/calibration/Qpeaks_{}.xlsx'.format(nsample)).T.values[0]  # .values

    else:
        # important sampling:
        if IS:
            # lowest VNK pf is 81, VNK discharges have been simualted up to
            # events of 125000
            low, high = werklijn_inv([1 - 1 / 1.002, 1 - 1 / 13000.0], A)
            samples = np.random.uniform(low, high, nsample)
            corr_f = werklijn_pdf(samples, A) / (1 / float(high - low))

        else:
            # normal sampling:
            # divide by the fraction reaching the IJssel
            samples = [rand_werklijn(A) for _ in range(0, nsample)]

        # Save samples:
        pd.DataFrame(samples).to_excel(
            './data/fragility_curves/calibration/Qpeaks_{}.xlsx'.format(nsample))
    # Save strengths:
#     random_numbers = np.random.uniform(size = nsample)
#     pd.DataFrame(random_numbers).to_excel('./data/fragility_curves/calibration/strengths_{}.xlsx'.format(nsample))

    inputs = {}
    factors = {'Bmax': 200, 'pfail': 0.5, 'Brate': 0.5}

    diversion = {'Diversion_401.0': -2.0, 'Diversion_201.0': 0.0}

    # Prepare input data:
    dike_nodes = DikeNetwork().dikenodes
    G = DikeNetwork().G

    for node in dike_nodes:
        inputs.update({'{}_{}'.format(node, fac)
                      : factors[fac] for fac in factors})

    inputs.update(diversion)

    # Compute high wate levels
    outcomes, G = DikeNetwork().__call__(Qpeaks=samples, **inputs)

#     pd.DataFrame(outcomes).to_excel('./data/fragility_curves/calibration/maxWL_{}_divpolicy.xlsx'.format(nsample))
    pd.DataFrame(outcomes).to_excel(
        './data/fragility_curves/calibration/maxWL_{}_{}_{}.xlsx'.format(nsample, int(low), int(high)))

pfailure = {n: 0 for n in dike_nodes}

if calibration:
    outcomes = pd.read_excel(
        './data/fragility_curves/calibration/maxWL_{}.xlsx'.format(nsample)).to_dict()
    # initialize the adjstment factors matrix
    guesses_df = pd.DataFrame(np.zeros((len(dike_nodes), 3)), index=dike_nodes,
                              columns=['low_guess', 'mean', 'up_guess'])

    for node in dike_nodes:
        guesses_df.loc[node]['low_guess'] = -10
        guesses_df.loc[node]['up_guess'] = 10

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
    print 'Calibrate Fragility Curves'
    while max(abs(guesses_df['low_guess'] - guesses_df['up_guess'])) > 0.0001 and any(
            abs(np.asarray(pfailure.values()) - np.asarray(pf_target.values())) / pf_target.values() > 0.0001):

        guesses_df['mean'] = guesses_df[['low_guess', 'up_guess']].mean(axis=1)
        frag_curve = deepcopy(nx.get_node_attributes(G, 'f'))

        for node in dike_nodes:
            p_failure = []
            frag_curve[node]['mu'] += guesses_df['mean'].loc[node]

            strength = norm.ppf(
                0.5,
                frag_curve[node]['mu'],
                frag_curve[node]['sd'])

            o = '{}_maxwl'.format(node)
            z = outcomes[o].values() > np.repeat(strength, nsample)

            if IS:
                #                   p_failure.append(np.sum(z*corr_f)/float(nsample))
                p_failure = np.sum(z * corr_f) / float(nsample)
            else:
                #                   p_failure.append(np.sum(z)/float(nsample))
                p_failure = np.sum(z) / float(nsample)

#                    pfailure[node] = np.mean(p_failure)
            pfailure[node] = p_failure

#            print (node, pfailure[node], guesses_df['low_guess'].loc[node],
#                   guesses_df['mean'].loc[node], guesses_df['up_guess'].loc[node],
#                   (pfailure[node] > pf_target), (pfailure[node] < pf_target))
#
#            if abs(pfailure[node] - pf_target[node])/pf_target[node] < 0.0001:
#                   print '{}_converged'.format(node)
#                   continue

            if pfailure[node] > pf_target[node]:  # here put a tolerance level
                guesses_df['low_guess'].loc[node] = guesses_df['mean'].loc[node]

            elif pfailure[node] < pf_target[node]:  # here put a tolerance level
                guesses_df['up_guess'].loc[node] = guesses_df['mean'].loc[node]

    correct = pd.DataFrame(nx.get_node_attributes(
        G, 'VNK_pf'), index=['correct']).T.dropna()
    estimated = pd.DataFrame({key: round(1 / pfailure[key]) for key in pfailure.keys()},
                             index=['estimated']).T
    p_estimation = pd.concat([correct, estimated], axis=1)

    guesses_df['NodeName'] = guesses_df.index
    guesses_df.index = range(len(dike_nodes))
    guesses_df[['mean', 'NodeName']].to_excel(
        './data/fragility_curves/calibration/FC_calfactors_{}_{}_IS{}.xlsx'.format(pf[pf_type], nsample, IS))

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
