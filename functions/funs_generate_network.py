# -*- coding: utf-8 -*-
"""
Created on Tue Oct 31 13:18:05 2017

@author: ciullo
"""

import numpy as np
import networkx as nx
import pandas as pd


def to_dict_dropna(data):
    return dict((str(k), v.dropna().to_dict())
                for k, v in pd.compat.iteritems(data))

def get_network(Musk_params=True, fragcurves=True, Losses=True):
    # Upload dike info
    df = pd.read_excel('./data/Rhine_network.xlsx', 
                       converters={'NodeName': str,
                                   'prec_node': str,
                                   'subs_node1': str,
                                   'subs_node2': str,
                                   'Brescode': int}, index_col=0)
    df = df.set_index('NodeName')

    nodes = df.to_dict('index')

    # Create network out of dike info
    G = nx.MultiDiGraph()
    for key, attr in nodes.items():
        G.add_node(key, **attr)

    # Nodes' type vectors:
    dike_nodes = df[df['type'] == 'dike'].index.values.tolist()
    bif_nodes = df[df['type'] == 'bifurcation'].index.values.tolist()

    # Branches dict:
    branches = {int(k / 100): df[df['BRANCH'] == k].index.values
                for k in np.unique(df['BRANCH'])}

    # Upload data:
    frag_curves = pd.read_excel('./data/fragility_curves/FC.xlsx', index_col=0)
    muskingum_params = pd.read_excel('./data/hydraulics/musk_calib_factor.xlsx',
                                     converters={'NodeName': str})
    muskingum_params = muskingum_params.set_index('NodeName')

    if fragcurves:
        cal_factors = pd.read_excel(
            './data/fragility_curves/calibration/FC_calfactors_VNK_3000_ISTrue.xlsx',
            converters={'NodeName': str}, index_col=0)
        cal_factors = cal_factors.set_index('NodeName')

    # Upload rfr:
    rfr = pd.read_excel('./data/measures/RfR.xlsx', index_col=0)

#    rfr.columns = range(rfr.shape[1])

    G.add_node('rfr', **to_dict_dropna(rfr))
    G.node['rfr']['type'] = 'measure'

    # Upload diversion policies:
    G.add_node('Diversion', **pd.read_excel(
        './data/measures/diversion_policies.xlsx', index_col=0).to_dict())
    G.node['Diversion']['type'] = 'measure'

    areas_to_rings = dict(zip(range(6), (['38', '24', '41', '42_nl'], ['43', '16'],
                                         ['15', '45', '44'],
                                         ['48_nl', '49', '50', '51', '52', '53'],
                                         ['42_de'], ['48_de'])))
    rings_to_areas = {}
    for i in range(6):
        [rings_to_areas.update({j: i}) for j in areas_to_rings[i]]

    tree = pd.DataFrame(columns=['dikering', 'area'], index=dike_nodes)

    # Fill network with crucial info:
    for curr_node in dike_nodes + bif_nodes:

        node = G.node[curr_node]
        # Assign Muskingum parameters per location
        node['C1'] = muskingum_params.loc[node['prec_node'], 'C1']
        node['C2'] = muskingum_params.loc[node['prec_node'], 'C2']
        node['C3'] = muskingum_params.loc[node['prec_node'], 'C3']

        if curr_node in dike_nodes:
            dikering = node['dikering']
            tree.loc[curr_node]['dikering'] = dikering
            tree.loc[curr_node]['area'] = rings_to_areas[dikering]

            # Assign fragility curves to breaching locations:
            if not np.isnan(node['Brescode']):
                node['f'] = frag_curves.loc[node['Brescode'], ['mu', 'sd']]

                # Adjust fragility curves
                if fragcurves:
                    node['f']['mu'] += cal_factors.loc[curr_node].values[0]

                node['dikelevel'] = node['f']['mu']

            # Assign stage-discharge curves
            filename = './data/hydraulics/rating_curves_extrap/{}_ratingcurve.txt'.format(
                                                                                curr_node)
            node['r'] = np.loadtxt(filename)

            if Losses:
                # Assign tables with damage in holland to all locations:
                name = './data/damages/{}_nl.txt'.format(curr_node)
                node['table_nl'] = np.loadtxt(name)

                # if trans. also attach damage in germany tables
                if dikering in ['42_de', '48_de', '42_nl', '48_nl']:
                    name = './data/damages/{}_de.txt'.format(curr_node)
                    node['table_de'] = np.loadtxt(name)

    # Plauisble wave-shapes for the most upstream location - from GRADE data
    upstream_node = df[df['type'] == 'upstream'].index.values[0]
    G.node[upstream_node]['Qevents_shape'] = pd.read_excel(
        './data/hydraulics/wave_shapes.xls', index_col=0)  # .to_dict('index')

    # Branch aggregation criteria:
#    agg_branches = ['BR+Waal', 'BR+Waal', 'PK+Lek', 'PK+Lek', 'IJssel']

#    corr = dict(zip([k*100+1 for k in branches.keys()], agg_branches))
    tree['location'] = tree.index
    tree.index = range(len(tree))

    transboundary_nodes = tree['location'][(tree['dikering'] == '42_de'
                               ) | (tree['dikering'] == '48_de'
                               ) | (tree['dikering'] == '42_nl'
                               ) | (tree['dikering'] == '48_nl')].values.tolist()

    national_nodes = tree['location'][(tree['dikering'] != '42_de'
                          ) & (tree['dikering'] != '48_de'
                          ) & (tree['dikering'] != '42_nl'
                          ) & (tree['dikering'] != '48_nl')].values.tolist()

    return G, branches, upstream_node, transboundary_nodes, national_nodes, tree
