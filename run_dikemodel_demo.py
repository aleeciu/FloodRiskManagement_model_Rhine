# -*- coding: utf-8 -*-
"""
Created on Mon Sep 23 09:32:46 2019

@author: ciullo
"""

import numpy as np
import pandas as pd
from dike_model_function import DikeNetwork 

network = DikeNetwork()

G = network.G
# print(G.nodes) # to see list of nodes (locations)
                 # import networkx if you want to interrogate the network

#print(G.node['201.880']) # example

tree = network.tree

unc_values = {'Bmax': 150, 'Brate': 1.5, 'pfail': 0.5}
scenario = {}

for unc in unc_values.keys():
    for loc in tree.location:
        scenario.update({'{}_{}'.format(loc, unc): unc_values[unc]})

measure_values = {'DikeIncrease': 0, 'rfr': 0, '201.0': 0, '401.0': 0}
policy = {}

for pol in measure_values.keys():
    if pol == 'DikeIncrease':
        for loc in tree.location:
            policy.update({'{}_{}'.format(loc, pol): measure_values[pol]})
    elif pol == 'rfr':
        policy.update({'{}_rfr'.format(i): measure_values[pol] for i in range(156)})
    else:
        policy.update({'Diversion_{}'.format(pol): measure_values[pol]})

#print(scenario)
#print(policy)

model_inputs = {**scenario, **policy}

results = network.__call__(**model_inputs)






