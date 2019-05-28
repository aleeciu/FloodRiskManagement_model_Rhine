# -*- coding: utf-8 -*-
"""
Created on Thu Jul 26 14:58:51 2018

@author: ciullo
"""

from __future__ import division
from functions.funs_generate_network import get_network
from functions.funs_hydrostat import werklijn_cdf
from functions.funs_economy import discount
import numpy as np
import pandas as pd


def risk_shift_optimization(init_risk, *args):

    EADa1 = args[0]
    EADb1 = np.sum(args[1:])

    EADa0, EADb0 = init_risk

    r1 = ((EADa0 - EADa1) / EADa0)
    r2 = ((EADb0 - EADb1) / EADb0)

    # distance from the r1 = r2 line
    distance = abs(r1 - r2) / np.sqrt(2)

    return distance


''' build file for maximum risk tranfers per area to input in the
 hypervolume computation '''
 
folder = r'C:\\Users\\ciullo\\OneDrive - Stichting Deltares\\Documents\\GitHub\\FloodRiskManagement_model_Rhine'

maximaals = pd.read_excel('{}/data/damages/max_damages_dikerings.xlsx'.format(folder),
                          index_col=0)
tree = get_network(Musk_params=True, fragcurves=True, Losses=True)[5]
#R0 = pd.read_excel('{}/data/pre_policies/zero_pol_wSB_n200_r35.xlsx'.format(folder))
A = pd.read_excel('{}/data/hydraulics/werklijn_params.xlsx'.format(folder))
Qpeak = np.unique(np.loadtxt(
            '{}/data/events/{}sampled_Qpeaks_125_12500.txt'.format(folder,10)))[::-1]
p_exc = 1 - werklijn_cdf(Qpeak, A)

## risk transfers upper bound (the lower bound is 0):
#areas = [0, 1, 2, 3, 4, 5]
#max_dist = {key: 0 for key in areas}
#
## note: these limits do not assume any constraint to risk transfers, in order words
## it allows a location to get (much) worse thank in the status quo.
#
#for area in areas:
#    # We assume all areas but the one under consideration have reached zero risk.
#    # In that area the the 'maximaal' damage is always reached
#    rest = [a for a in areas if not a == area]
#    rings = np.unique(tree['dikering'][tree['area'] == area]).tolist()
#    damages = [maximaals.loc[rings, 'max_damage'].sum()]
#    ead_area = np.trapz(damages * len(p_exc), p_exc)
#    eaddisc_area = np.sum(discount(ead_area, rate=3.5, n=200))
#
#    # if you want to simulate a constraint, worst case is risk equal to init value
##    eaddisc_area = R0['{}_EAD'.format(area)].values
#
#    init_risk = R0['{}_EAD'.format(area)].values.tolist() + [np.sum(
#        [R0['{}_EAD'.format(a)][0] for a in rest])]
#    ead_rest = [0] * (len(areas) - 1)
#    max_dist[area] = risk_shift_optimization(
#        init_risk, *[eaddisc_area] + ead_rest)
#    print(area, max_dist[area], eaddisc_area / init_risk[0])

# note: maximum distances values are the same for PF2 and PF2, even though they
# are calculated differently in the optimization process.
# It represents the ''physical boundary'' (practically: In anycase the rest risk 
# level is set to zero)
#pd.DataFrame(max_dist, index = [0]).to_excel('./data/pre_policies/max_dists.xlsx')

## max investment costs (all RfR + dike raising = 1 m)
max_InvC = {u'Investment Costs_de': 323680579.1322484,
            u'Investment Costs_nl': 19029010167.40925,
            u'Total Costs_de': 685037758.2952368,
            u'Total Costs_nl': 18845199466.77497}

## by putting a contraint on inv costs:
#max_InvC = {u'Investment Costs_de': 2*1e8,
#            u'Investment Costs_nl': 2*1e8,
#            u'Total Costs_de': 685037758.2952368,
#            u'Total Costs_nl': 18845199466.77497}

# max expected annual damage per country:
nl = [0, 1, 2, 3]
de = [4, 5]

nl_ead = []
de_ead = []

for area in nl + de:
    # Assuming that the maximaal is always reached

    if area in nl:
        rings = np.unique(tree['dikering'][tree['area'] == area]).tolist()
        damages = [maximaals.loc[rings, 'max_damage'].sum()]
        ead_area = np.trapz(damages * len(p_exc), p_exc)
        nl_ead.append(np.sum(discount(ead_area, rate=3.5, n=200)))

    elif area in de:
        rings = np.unique(tree['dikering'][tree['area'] == area]).tolist()
        damages = [maximaals.loc[rings, 'max_damage'].sum()]
        ead_area = np.trapz(damages * len(p_exc), p_exc)
        de_ead.append(np.sum(discount(ead_area, rate=3.5, n=200)))

max_nl = np.sum(nl_ead)
max_de = np.sum(de_ead)

max_costs = {'de': np.max([max_InvC['Investment Costs_de'], max_de]),
             'nl': np.max([max_InvC['Investment Costs_nl'], max_nl])}

# In Germany is the maximum expected annual damage leading, in the netherlands
# are investment costs

#pd.DataFrame(max_costs, index = [0]).to_excel('./data/pre_policies/max_costs.xlsx')
