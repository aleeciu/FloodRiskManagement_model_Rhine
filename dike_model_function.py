# -*- coding: utf-8 -*-
"""
Created on Tue Oct 31 13:18:05 2017

@author: ciullo
"""

from copy import deepcopy

from functions.funs_generate_network import get_network
from functions.funs_dikes import Lookuplin, dikefailure, init_node
from functions.funs_economy import cost_dike, discount
from functions.funs_hydrostat import werklijn_cdf, werklijn_inv
import numpy as np
import pandas as pd
from scipy.stats import norm


def Muskingum(C1, C2, C3, Qn0_t1, Qn0_t0, Qn1_t0):
    Qn1_t1 = C1 * Qn0_t1 + C2 * Qn0_t0 + C3 * Qn1_t0

    return Qn1_t1


class DikeNetwork(object):
    def __init__(self):
        # load network
        network = get_network(Musk_params=True, fragcurves=True, Losses=True)

        G = network[0]
        branches = network[1]
        upstream_node = network[2]
        transboundary_nodes = network[3]
        national_nodes = network[4]
        tree = network[5]

        # # of events:
        n = 10
        R0 = pd.read_excel(
            './data/pre_policies/zero_pol_wSB_n200_r35_{}.xlsx'.format(n))
        
        # maximum damage per dikering:
        maximaals = pd.read_excel('./data/damages/max_damages_dikerings.xlsx')
        
        # Load hydrologic statistic
        A = pd.read_excel('./data/hydraulics/werklijn_params.xlsx')

#        low, high = werklijn_inv([1-1/125.0, 1-1/12500.0], A)
#        Qpeak = np.unique(np.random.uniform(low, high, n))[::-1]
#        np.savetxt('./data/events/{}uniform_sampled_Qpeaks_125_12500.txt'.format(n), Qpeak)

#        sample_probs = np.random.uniform(1-1/125.0, 1-1/12500.0, n)
#        Qpeak = np.unique(werklijn_inv(sample_probs, A))[::-1]
#        np.savetxt('./data/events/{}sampled_Qpeaks_125_12500.txt'.format(n), Qpeak)

#        Qpeak = np.unique(np.loadtxt(
#                './preprocessing/10uniform_sampled_Qpeaks_125_12500.txt'))[::-1]

        Qpeak = np.unique(np.loadtxt(
            './data/events/{}sampled_Qpeaks_125_12500.txt'.format(n)))[::-1]

        p_exc = 1 - werklijn_cdf(Qpeak, A)  # probabiltiy of exceedence

        self.R0 = R0
        self.Qpeak = Qpeak
        self.p_exc = p_exc
        self.A = A
        self.G = G
        self.maximaals = maximaals
        self.transboundary_nodes = transboundary_nodes
        self.national_nodes = national_nodes
        self.dikenodes = transboundary_nodes + national_nodes
        self.branches = branches
        self.upstream_node = upstream_node
        self.tree = tree

        self.sb = True
        self.n = 200
        self.rate = 3.5  # discount rate
        self.step = 10  # dike increase step [cm]

        # Time step correction: Q is a mean daily value expressed in m3/s
        self.timestepcorr = 24 * 60 * 60

    def _initialize_hydroloads(self, node, time, Q_0):
        node['cumVol'], node['wl'], node['Qpol'], node['hbas'] = (
            init_node(0.0, time) for _ in range(4))
        node['Qin'], node['Qout'] = (init_node(Q_0, time) for _ in range(2))
        node['status'] = init_node(False, time)
        node['tbreach'] = np.nan
        return node

    def _initialize_rfr__ooi(self, G, dikenodes, country):
        for n in dikenodes:
            node = G.node[n]
            # Create a copy of the rating curve that will be used in the sim:
            node['rnew'] = deepcopy(node['r'])

            # Initialize ooi:
            node['losses_{}'.format(country)] = []
            node['deaths_{}'.format(country)] = []
        G.node['rfr']['totcost'] = 0
        return G

    def __call__(self, timestep=1, div=0.2, q=4, **kwargs):

        G = self.G
        upst_node = self.upstream_node
        tree = self.tree
        dikenodes = self.dikenodes

        self._initialize_rfr__ooi(G, self.national_nodes, 'nl')
        self._initialize_rfr__ooi(G, self.transboundary_nodes, 'nl')
        self._initialize_rfr__ooi(G, self.transboundary_nodes, 'de')

        # load all kwargs into network
        for item in kwargs:
            string1, string2 = item.split('_')
            # string1: usually dikename
            # string2: usually name of uncertainty or lever
            # exception only for rfr strategies where: 'project number id_rfr'

            # If it's a RfR measure:
            if string2 == 'rfr':
                rfr_project_node = G.node['rfr'][string1]
                # add costs current project to the whole costs
                G.node['rfr']['totcost'] += kwargs[item] * \
                    rfr_project_node['Costs_1e6'] * 1e6

                # apply wl reduction to the locations involved by this measure
                for key in rfr_project_node.keys():
                    if key != 'Costs_1e6':
                        G.node[key]['rnew'][:,
                                            1] += kwargs[item] * rfr_project_node[key]

            # If it's DikeIncrease or Diversion
            else:
                # get something like G.node[DikeIncrease][] or
                # G.node[Diversion][]
                G.node[string1][string2] = kwargs[item]

                # if it's dike rasing:
                if string2 == 'DikeIncrease':
                    node = G.node[string1]

                    node['fnew'] = deepcopy(node['f'])

                    node[string2] = (kwargs[item] * self.step) / 100.0  # cm
                    node['fnew']['mu'] += node[string2]

                    if node[string2] == 0:
                        node['dikecosts'] = 0
                    else:
                        node['dikecosts'] = cost_dike(
                            node['c1'], node['b1'], node['lambda1'], node['ratio1'],
                            node['c2'], node['b2'], node['lambda2'], node['ratio2'],
                            node['c3'], node['b3'], node['lambda3'], node['ratio3'],
                            node[string2])

        # Dictionary storing outputs:
        ring_output = {'{}_{}'.format(dr, key): [] for key in [
                       'Damage_de', 'Damage_nl', 'Deaths_de', 'Deaths_nl', 'Dike Inv Cost'
                       ] for dr in np.unique(tree['dikering'])}

        area_output = {'{}_{}'.format(a, key): [] for key in [
                       'Damage', 'Deaths', 'EAD', 'Dike Inv Cost'] for a in range(6)}

        # Possible extra outputs of interest:
#        extra_output_list = ['Qpol', 'Qout', 'Qin', 'wl', 'status', 'critWL']
#        # Selected extra outputs of interest
#        eooi = [3]
#        for k in eooi:
#            extra_output = {'{}_{}'.format(extra_output_list[k],
#                                            dike): [] for dike in dikenodes}

        for Q in self.Qpeak:
            # THE upstream event:
            time = np.arange(0, G.node[upst_node]['Qevents_shape'].loc[q].shape[0],
                             timestep)
            G.node[upst_node]['Qin'] = Q * \
                G.node[upst_node]['Qevents_shape'].loc[q]

            # Select a branch:
            for j in self.branches.keys():
                branchnodes = self.branches[j]

                # Initialize upstream hydrpgraph
                G.node[branchnodes[0]]['Qout'] = G.node[G.node[branchnodes[0]][
                    'prec_node']]['Qin']
                Q_0 = int(G.node[branchnodes[0]]['Qout'][0])

                # Impose initial conditions at each node:
                for n in branchnodes:
                    node = G.node[n]

                    if node['type'] == 'dike':
                        self._initialize_hydroloads(node, time, Q_0)
                    # this could go before the loop, but we are not sure if pfail
                    # is iterated before or after dikeincrease:
                        node['critWL'] = norm.ppf(node['pfail'], node['fnew']['mu'],
                                                  node['fnew']['sd'])

                    elif node['type'] == 'bifurcation':

                        # Assign the diversion measure to the corresponding
                        # branches:
                        Qbif = '{}_Q'.format(n)
                        # subs_node1 relate to the uneven branches (301, 501),
                        # subs_node2 to the even ones (201, 401)
                        r1 = '{}_ratio'.format(node['subs_node1'])
                        r2 = '{}_ratio'.format(node['subs_node2'])

                        G.node[node['subs_node2']]['div_measure'] = pd.DataFrame(
                            {k: list(G.node['Diversion'][k].values())
                             for k in [Qbif, r2]},
                            columns=[r2, Qbif]).values

                        # Apply the diversion policy to the even branches: a diversion
                        # policy is the percentage increase in the fraction of water sent
                        # to a branch. If status quo is 29% of the incoming water sent
                        # to a branch, a policy of +80% would set it at 52%
                        # (=0.29 * (1 + 0.8))
                        G.node[node['subs_node2']]['div_measure'][:, 0] *= 1 + G.node[
                            'Diversion']['{}'.format(node['subs_node2'])] / 10.0

                        G.node[node['subs_node1']]['div_measure'] = pd.DataFrame(
                            {k: list(G.node['Diversion'][k].values())
                             for k in [Qbif, r1]},
                            columns=[r1, Qbif]).values
                        G.node[node['subs_node1']]['div_measure'][:, 0] = 1 - G.node[
                            node['subs_node2']]['div_measure'][:, 0]

                        # Qin for the three nodes of the bifurcation:
                        # main node
                        node['Qin'] = (init_node(Q_0, time))

                        G.node[node['subs_node1']]['Qin'] = Lookuplin(
                            G.node[node['subs_node1']]['div_measure'], 1, 0,
                            node['Qin']) * node['Qin']

                        G.node[node['subs_node2']]['Qin'] = Lookuplin(
                            G.node[node['subs_node2']]['div_measure'], 1, 0,
                            node['Qin']) * node['Qin']

                # Run over each node of the branch:
                for n in range(0, len(branchnodes)):
                    node = G.node[branchnodes[n]]

                    # Propagate the discharge wave:
                    for t in range(1, len(time)):

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
                                node['rnew'], 0, 1, [node['Qin'][t]])

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

                            node['cumVol'][t] = np.trapz(
                                node['Qpol']) * self.timestepcorr

                            # the 1/125.0 area is considered the minimum one:
                            Area = Lookuplin(
                                node['table_nl'], 3, 0, [
                                    node['wl'][t]])

                            if branchnodes[n] in self.transboundary_nodes:
                                Area += Lookuplin(node['table_de'],
                                                  3, 0, [node['wl'][t]])

                            node['hbas'][t] = node['cumVol'][t] / float(Area)

                        elif node['type'] == 'bifurcation':

                            # Distribute Qin into the two branches. Note:
                            # Qins and Qouts identical, no breach can happen.

                            G.node[node['subs_node1']]['Qin'][t] = Lookuplin(
                                G.node[node['subs_node1']]['div_measure'],
                                1, 0, [node['Qin'][t]]) * node['Qin'][t]

                            G.node[node['subs_node2']]['Qin'][t] = Lookuplin(
                                G.node[node['subs_node2']]['div_measure'],
                                1, 0, [node['Qin'][t]]) * node['Qin'][t]

        # iterate over the dike_nodes and save outputs
            for n in dikenodes:
                node = G.node[n]

                # append losses for each Q event:
                # each location causes a damage in the NL
                node['losses_nl'].append(Lookuplin(node['table_nl'],
                             3, 1, [np.max(node['wl'])]) * node['status'][-1])

#                node['deaths_nl'].append(Lookuplin(node['table_nl'],
#                        3, 2, np.max(node['wl']))*node['status'][-1])

                # if you are on a transb. dikering, then you also have damage
                # in DE
                if n in self.transboundary_nodes:
                    node['losses_de'].append(Lookuplin(node['table_de'],
                         3, 1, [np.max(node['wl'])]) * node['status'][-1])

#                    node['deaths_de'].append(Lookuplin(node['table_de'],
#                         3, 2, np.max(node['wl']))*node['status'][-1])

                ## Extra outcomes of interest:                    
#                for k in eooi:
#                    o = extra_output_list[k]
#                    key = '{}_{}'.format(o,n)
#
#                    if not o in ['wl', 'status']:
#                        ## usually scalar, when series (wl) we want the max value
#                         extra_output[key].append(node[o])
#                    else:
#                         extra_output[key].append(np.max(node[o]))

        # TODO: the following two loops look pretty similar: make a function
        # iterate over the dike rings and save outputs
        for ring in np.unique(tree['dikering']):
            locations = tree['location'][tree['dikering'] == ring]
            for l in locations:

                ring_output['{}_Dike Inv Cost'.format(ring)].append(
                    G.node[l]['dikecosts'])

                # location on dike rings 42 and 48 will also cause damage in DE
                if l in self.transboundary_nodes:
                    # either you are on eg dk ring 42_nl or 42_de, save the german
                    # damage as you were in 42_de...

                    _ring = '{}_de'.format(ring.split('_')[0])
                    ring_output['{}_Damage_de'.format(_ring)].append(
                        G.node[l]['losses_de'])

                    # ..and the other way around with the dutch damage
                    _ring = '{}_nl'.format(ring.split('_')[0])
                    ring_output['{}_Damage_nl'.format(_ring)].append(
                        G.node[l]['losses_nl'])

                else:
                    ring_output['{}_Damage_nl'.format(ring)].append(
                        G.node[l]['losses_nl'])

            ring_output['{}_Dike Inv Cost'.format(ring)] = np.nansum(
                ring_output['{}_Dike Inv Cost'.format(ring)])

        # iterate over the aggregated areas and save outputs
        for area in range(6):
            rings = np.unique(tree['dikering'][tree['area'] == area])

            for ring in rings:
                area_output['{}_Dike Inv Cost'.format(area)].append(
                    ring_output['{}_Dike Inv Cost'.format(ring)])

                # german areas
                if area in [4, 5]:

                    damage_in_ring = np.minimum(
                                     self.maximaals.loc[str(ring)]['max_damage'],
                                     np.nansum(ring_output['{}_Damage_de'.format(ring)], 
                                                           axis=0))
                    
                    area_output['{}_Damage'.format(area)].append(damage_in_ring)

            # dutch areas
                else:
                    damage_in_ring = np.minimum(
                                     self.maximaals.loc[str(ring)]['max_damage'],
                                     np.nansum(ring_output['{}_Damage_nl'.format(ring)], 
                                                              axis=0))
                    
                    area_output['{}_Damage'.format(area)].append(damage_in_ring)

            area_output['{}_Dike Inv Cost'.format(area)] = np.nansum(
                                area_output['{}_Dike Inv Cost'.format(area)])
            
            area_output['{}_Damage'.format(area)] = np.nansum(area_output[
                                            '{}_Damage'.format(area)], axis=0)
                                
            # Expected Annual Damage:
            _damage = np.reshape(
                area_output['{}_Damage'.format(area)], self.p_exc.shape)
            
            ## Doing this doesnt make much sense, as you may have 
            ## ead (area under curve) > max risk
            
#            _risks = _damage*self.p_exc
#            
#            area_output['{}_minR'.format(area)] = np.min(_risks)
#            
#            area_output['{}_maxR'.format(area)] = np.max(_risks)  
            
            EAD = np.trapz(_damage, self.p_exc)
            # Discounted annual risk per dike ring:
            area_output['{}_EAD'.format(area)] = np.sum(discount(EAD,
                                                    rate=self.rate, n=self.n))
            
#        area_output.update(extra_output)

        area_output.update({'RfR Total Costs': G.node['rfr']['totcost']})

        return area_output







