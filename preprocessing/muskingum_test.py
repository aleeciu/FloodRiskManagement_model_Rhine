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


def Muskingum(C1, C2, C3, Qn0_t1, Qn0_t0, Qn1_t0):
    Qn1_t1 = C1 * Qn0_t1 + C2 * Qn0_t0 + C3 * Qn1_t0

    return Qn1_t1


class DikeNetwork(object):
    def __init__(self):
        # load network
        G, branches, upstream_node, dike_nodes, bif_nodes = get_network(Musk_params=True,
                                                                        fragcurves=False,
                                                                        Losses=True)
        self.G = G
        self.bifnodes = bif_nodes
        self.dikenodes = dike_nodes
        self.branches = branches
        self.upstream_node = upstream_node

        self.sb = False
        ema_logging.info('model initialized')

    def _initialize_hydroloads(self, node, time, Q_0):
        node['cumVol'], node['wl'], node['Qpol'], node['hbas'] = (
            init_node(0, time) for _ in range(4))
        node['Qin'], node['Qout'] = (init_node(Q_0, time) for _ in range(2))
        node['status'] = init_node(False, time)
        node['tbreach'] = np.nan
        return node

    def _initialize_rfr(self, G, dikenodes):
        # Initialize room for the river
        G.node['rfr_projects']['rfr_cost'] = 0
        # Create a copy of the rating curve that will be used in the
        # simulation:
        for node in dikenodes:
            G.node[node]['rnew'] = deepcopy(G.node[node]['r'])
        return G

    def __call__(self, timestep=1, div_frac=None, q=4, Q=20000, **kwargs):

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
        Qout = {dike: [] for dike in dikenodes}
#        Qin = {dike: [] for dike in dikenodes}
#        maxwl = {dike: [] for dike in dikenodes}
#        status = {dike: [] for dike in dikenodes}

        # THE upstream event:
#            time = np.arange(0, G.node[upst_node]['Qevents_shape'].loc[q].shape[0], timestep)
#            G.node[upst_node]['Qin'] = Q*G.node[upst_node]['Qevents_shape'].loc[q]

        Qq = pd.read_excel(
            './data/hydraulics/SOBEK_discharges/A_branch.xlsx'.format(fl)).iloc[:, [1]]
        G.node[upst_node]['Qin'] = Qq.groupby(
            (Qq.index / 24).astype(int)).mean().values
        time = range(len(G.node[upst_node]['Qin']))

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

                    G.node[node['subs_node1']]['div_measure'] = pd.DataFrame(
                        {k: G.node['diversion'][k].values()
                         for k in [Qbif, r1]},
                        columns=[r1, Qbif]).values

                    G.node[node['subs_node2']]['div_measure'] = pd.DataFrame(
                        {k: G.node['diversion'][k].values()
                         for k in [Qbif, r2]},
                        columns=[r2, Qbif]).values

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

                    if node['type'] in ['upstream', 'downstream']:
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
                                          node['status'][t - 1], node['Bmax'],
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
                            node['table'].values, 4, 0, node['wl'][t])
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

            # Neglect first and last locations of the branch: not dike type
            for n in branchnodes[1:-1]:
                node = G.node[n]

                Qout[n].append(node['Qout'])
#                        Qin[n].append(node['Qin'])
#
#                        Qpol[n].append(node['Qpol'])
#                        maxwl[n].append(np.max(node['wl']))
#                        status[n].append(node['status'][-1])

        for n in dikenodes:
            data.update({
                '{}_Qout'.format(n): Qout[n],
                #                             '{}_Qin'.format(n): Qin[n],
                #                             '{}_maxwl'.format(n): maxwl[n],
                #                             '{}_Qpol'.format(n): Qpol[n],
                #                             '{}_status'.format(n): status[n]
            })

        return data, G


import matplotlib.pyplot as plt

''' See if dike model function properly simulates Muskingum as in muskingum_calibration.py'''

discharges = pd.DataFrame()
for fl in ['A', 'C', 'D', 'E']:
    if fl == 'A':
        data = pd.read_excel(
            './data/hydraulics/SOBEK_discharges/{}_branch.xlsx'.format(fl)).iloc[:, [1]]
    else:
        data = pd.read_excel(
            './data/hydraulics/SOBEK_discharges/{}_branch.xlsx'.format(fl)).iloc[:, [-2]]

    data = data.drop(
        data.index[-1]).groupby((data.index[:-1] / 24).astype(int)).mean()
    discharges = pd.concat([discharges, data], axis=1)
discharges.columns = ['101.860', '201.959', '401.973', '501.960']


# Prepare input data:
dike_nodes = DikeNetwork().dikenodes
inputs = {}
factors = {'Bmax': 200, 'pfail': 0.5, 'Brate': 0.5}

for dike in dike_nodes:
    inputs.update({'{}_{}'.format(dike, fac): factors[fac] for fac in factors})
#
#
outcomes, G = DikeNetwork().__call__(**inputs)
#
for l in ['201.959', '401.973', '501.960']:
    o = '{}_Qout'.format(l)
    plt.plot(range(len(outcomes[o][0])), outcomes[o][0], '--')
    plt.title('{}'.format(o.split('_')[0]))

    plt.plot(range(len(discharges[l])), discharges[l])

    plt.show()
