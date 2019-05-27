# -*- coding: utf-8 -*-
"""
Created on Mon Mar 26 10:23:02 2018

@author: ciullo
"""
from __future__ import division
from copy import deepcopy
from ema_workbench import ema_logging

import pandas as pd
from scipy.stats import norm

from funs_generate_network import get_network
from funs_dikes import Lookuplin, dikefailure, init_node
import numpy as np
from functools import partial
import matplotlib.pyplot as plt
from funs_hydrostat import rand_werklijn, werklijn_pdf, werklijn_inv
import networkx as nx


def Muskingum(C1, C2, C3, Qn0_t1, Qn0_t0, Qn1_t0):
    Qn1_t1 = C1 * Qn0_t1 + C2 * Qn0_t0 + C3 * Qn1_t0

    return Qn1_t1


class DikeNetwork(object):
    def __init__(self):
        # load network
        G, branches, upstream_node, dike_nodes, bif_nodes, corr = get_network(Musk_params=True,
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

                        G.node[node['subs_node1']]['div_measure'] = pd.DataFrame(
                            {k: G.node['Diversion'][k].values()
                             for k in [Qbif, r1]},
                            columns=[r1, Qbif]).values

                        G.node[node['subs_node2']]['div_measure'] = pd.DataFrame(
                            {k: G.node['Diversion'][k].values()
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
                                              node['status'][t -
                                                             1], node['Bmax'],
                                              node['Brate'], time[t],
                                              node['tbreach'], node['critWL'])

                            node['Qout'][t] = res[0]
                            node['Qpol'][t] = res[1]
                            node['status'][t] = res[2]
                            node['tbreach'] = res[3]

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


A = pd.read_excel('./data/hydraulics/werklijn_params.xlsx')
#FC = pd.read_excel('./data/fragility_curves/FC.xlsx')

# Prepare input data:
dike_nodes = DikeNetwork().dikenodes
inputs = {}
factors = {'Bmax': 200, 'pfail': 0.5, 'Brate': 0.5}

for node in dike_nodes:
    inputs.update({'{}_{}'.format(node, fac): factors[fac] for fac in factors})

output = {'{}_maxwl'.format(d): [] for d in dike_nodes}

for RT in [125.0, 1250.0, 12500.0, 125000.0]:
    Q = werklijn_inv([1 - 1 / RT], A)
    o, G = DikeNetwork().__call__(Qpeaks=Q, **inputs)

    [output[d].extend(o[d]) for d in output.keys()]

''' Build losses file'''


# TODO: Tentative areas, could not load the VNK scenario maps.
area = np.array([1.3 * 1e8, 1.5 * 1e8, 1.6 * 1e8, 1.75 * 1e8])
volume = np.array([2.8 * 1e7, 6.1 * 1e7, 9.3 * 1e7, 2 * 1e8])
h = volume / area

tp = ['-1', '0', '+1', '+2']
pos = [0, 1, 2, 3]

corr = dict(zip(pos, tp))


def back_damage(Din, Dtpcm, cm):
    Dout = Din - Dtpcm * cm
    return Dout


def forw_damage(Din, Dtpcm, cm):
    Dout = Din + Dtpcm * cm
    return Dout


damages = pd.read_excel('./data/damages/VNK_damages.xls', index_col=[0],
                        usecols="B:X", skiprows=[1])

# '0.1', '1.2', '2.3' eg 0.1 increase rate from pos 0 to pos 1 (TP-1 to TP)
damages.columns = ['{}_Schade'.format(TP) for TP in ['-1', '0', '+1', '+2']
                   ] + ['{}_Schade_per_cm'.format(TP) for TP in ['0.1', '1.2', '2.3']
                        ] + ['{}_Slachtoffers'.format(TP) for TP in ['-1', '0', '+1', '+2']
                             ] + ['{}_Slachtoffers_per_cm'.format(TP) for TP in ['0.1', '1.2', '2.3']
                                  ] + ['{}_Getroffenen'.format(TP) for TP in ['-1', '0', '+1', '+2']
                                       ] + ['{}_Getroffenen_per_cm'.format(TP) for TP in ['0.1', '1.2', '2.3']
                                            ] + ['Decimeringshoogte (cm)']

# TODO: make a function instead of doing three times the same thing in the loop

for node in dike_nodes:
    naam = G.node[node]['VNKNAAM']
    # deaths
    # convert all zeros into a reasonable number -> copy the decim nearby
    decim = damages.loc[naam, ['{}_Slachtoffers_per_cm'.format(TP) for TP in [
        '0.1', '1.2', '2.3']]]
    decim = decim.replace(0, np.nan)
    decim = decim.fillna(method='ffill')
    decim = decim.fillna(method='bfill')
    # identify where nan and non nan are:
    deaths = damages.loc[naam,
                         ['{}_Slachtoffers'.format(TP) for TP in tp]].values
    nans_bool = np.isnan(deaths)
    pos_nans = list(np.where(nans_bool)[0])            # position of nans
    non_nan = list(np.where(nans_bool == False)[0])      # position of non nans

#        print(damages.loc[naam, ['{}_Slachtoffers'.format(TP) for TP in tp]])
    if pos_nans == range(4):
        continue
    else:
        while pos_nans:  # while you have nans
            for pos_nan in pos_nans:
                # if nan has an upper valid damage estimation
                if (pos_nans[0] + 1) in non_nan:
                    damages.loc[naam, '{}_Slachtoffers'.format(corr[pos_nan])] = back_damage(
                        damages.loc[naam,
                                    '{}_Slachtoffers'.format(corr[pos_nan + 1])],
                        decim['{}.{}_Slachtoffers_per_cm'.format(
                            pos_nan, pos_nan + 1)],
                        damages.loc[naam, 'Decimeringshoogte (cm)'])

                    # drop it from the nan and append it to the non nan
                    pos_nans.remove(pos_nan)
                    non_nan.append(pos_nan)

                # if nan has a lower valid damage estimation
                elif (pos_nans[0] - 1) in non_nan:
                    damages.loc[naam, '{}_Slachtoffers'.format(corr[pos_nan])] = forw_damage(
                        damages.loc[naam,
                                    '{}_Slachtoffers'.format(corr[pos_nan - 1])],
                        decim['{}.{}_Slachtoffers_per_cm'.format(
                            pos_nan - 1, pos_nan)],
                        damages.loc[naam, 'Decimeringshoogte (cm)'])

                    # drop it from the nan and append it to the non nan
                    pos_nans.remove(pos_nan)
                    non_nan.append(pos_nan)

#        print(damages.loc[naam, ['{}_Slachtoffers'.format(TP) for TP in tp]])

    # losses
    decim = damages.loc[naam, ['{}_Schade_per_cm'.format(TP) for TP in [
        '0.1', '1.2', '2.3']]]
    decim = decim.replace(0, np.nan)
    decim = decim.fillna(method='ffill')
    decim = decim.fillna(method='bfill')

    deaths = damages.loc[naam, ['{}_Schade'.format(TP) for TP in tp]].values
    nans_bool = np.isnan(deaths)
    pos_nans = list(np.where(nans_bool)[0])            # position of nans
    non_nan = list(np.where(nans_bool == False)[0])      # position of non nans

#        print(damages.loc[naam, ['{}_Schade'.format(TP) for TP in tp]])
    if pos_nans == range(4):
        continue
    else:
        while pos_nans:  # while you have nans
            for pos_nan in pos_nans:
                if (pos_nans[0] + 1) in non_nan:
                    damages.loc[naam, '{}_Schade'.format(corr[pos_nan])] = back_damage(
                        damages.loc[naam,
                                    '{}_Schade'.format(corr[pos_nan + 1])],
                        decim['{}.{}_Schade_per_cm'.format(
                            pos_nan, pos_nan + 1)],
                        damages.loc[naam, 'Decimeringshoogte (cm)'])

                    pos_nans.remove(pos_nan)
                    non_nan.append(pos_nan)
                elif (pos_nans[0] - 1) in non_nan:
                    damages.loc[naam, '{}_Schade'.format(corr[pos_nan])] = forw_damage(
                        damages.loc[naam,
                                    '{}_Schade'.format(corr[pos_nan - 1])],
                        decim['{}.{}_Schade_per_cm'.format(
                            pos_nan - 1, pos_nan)],
                        damages.loc[naam, 'Decimeringshoogte (cm)'])

                    pos_nans.remove(pos_nan)
                    non_nan.append(pos_nan)

#        print(damages.loc[naam, ['{}_Schade'.format(TP) for TP in tp]])

    # affected
    decim = damages.loc[naam, ['{}_Getroffenen_per_cm'.format(TP) for TP in [
        '0.1', '1.2', '2.3']]]
    decim = decim.replace(0, np.nan)
    decim = decim.fillna(method='ffill')
    decim = decim.fillna(method='bfill')

    deaths = damages.loc[naam,
                         ['{}_Getroffenen'.format(TP) for TP in tp]].values
    nans_bool = np.isnan(deaths)
    pos_nans = list(np.where(nans_bool)[0])            # position of nans
    non_nan = list(np.where(nans_bool == False)[0])      # position of non nans

#        print(damages.loc[naam, ['{}_Getroffenen'.format(TP) for TP in tp]])
    if pos_nans == range(4):
        continue
    else:
        while pos_nans:  # while you have nans
            for pos_nan in pos_nans:
                if (pos_nans[0] + 1) in non_nan:
                    damages.loc[naam, '{}_Getroffenen'.format(corr[pos_nan])] = back_damage(
                        damages.loc[naam,
                                    '{}_Getroffenen'.format(corr[pos_nan + 1])],
                        decim['{}.{}_Getroffenen_per_cm'.format(
                            pos_nan, pos_nan + 1)],
                        damages.loc[naam, 'Decimeringshoogte (cm)'])

                    pos_nans.remove(pos_nan)
                    non_nan.append(pos_nan)
                elif (pos_nans[0] - 1) in non_nan:
                    damages.loc[naam, '{}_Getroffenen'.format(corr[pos_nan])] = forw_damage(
                        damages.loc[naam,
                                    '{}_Getroffenen'.format(corr[pos_nan - 1])],
                        decim['{}.{}_Getroffenen_per_cm'.format(
                            pos_nan - 1, pos_nan)],
                        damages.loc[naam, 'Decimeringshoogte (cm)'])

                    pos_nans.remove(pos_nan)
                    non_nan.append(pos_nan)

    deaths = damages.loc[naam, ['{}_Slachtoffers'.format(TP) for TP in [
                         '-1', '0', '+1', '+2']]]
#        print(deaths)

    losses = damages.loc[naam, ['{}_Schade'.format(TP) for TP in tp]]
#        print(losses)

    affected = damages.loc[naam, ['{}_Getroffenen'.format(TP) for TP in tp]]
#        print(affected)

    wl = output['{}_maxwl'.format(node)]
    data = np.array([area, volume, h, deaths.values,
                     losses.values, affected.values, wl]).T

    flname = './data/damages/{}_lossestable.txt'.format(node)
    np.savetxt(flname, data)
