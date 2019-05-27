# -*- coding: utf-8 -*-
"""
Created on Thu Feb 08 10:19:52 2018

@author: ciullo
"""

''' ADJUST EXCEL NETWORK FILE '''

import pandas as pd
import numpy as np

df = pd.read_excel('./data/Rhine_info.xls')

# order according to branch
df = df.sort_values('BRANCH')

# iteratively extract nodes of a given branch, order by km and concat them back
branches = np.unique(df['BRANCH'])

new_df = pd.DataFrame(columns=df.columns)

# make df where locations come in blocks of branches number and orderd by km
for branch in branches:
    dummy_df = df[df['BRANCH'] == branch].sort_values('KILOMETER')
    new_df = pd.concat([new_df, dummy_df])

# order indexes
new_df.index = range(len(new_df.index))

# load VNK probability of failure
VNKpf = pd.read_excel('./data/initial_pfs/VNK_Kans.xlsx', index_col=[0])

# assign types
upstream = ['91.750', '201.0', '301.0', '401.0', '501.0']
conjunction = ['101.0']
bifurcation = {'101.1000': {1: '201.0', 2: '301.0'},
               '301.1000': {1: '401.0', 2: '501.0'}}
downstream = ['201.960', '401.977', '501.977']

for i in list(new_df.index):

    new_df.loc[i, 'NodeName'] = '{}.{}'.format(new_df.loc[i, 'BRANCH'],
                                               new_df.loc[i, 'KILOMETER'])
    node = new_df.loc[i, 'NodeName']

    if new_df.loc[i, 'traject'] is not np.nan:
        if len(new_df.loc[i, 'traject'].split('_de')) > 1:
            new_df.loc[i, 'dikering'] = '{}_{}'.format(new_df.loc[i, 'traject'
                                                                  ].split('_')[0], 'de')

        elif len(new_df.loc[i, 'traject'].split('_nl')) > 1:
            new_df.loc[i, 'dikering'] = '{}_{}'.format(new_df.loc[i, 'traject'
                                                                  ].split('_')[0], 'nl')
        else:
            new_df.loc[i, 'dikering'] = new_df.loc[i, 'traject'].split('_')[0]
    else:
        new_df.loc[i, 'dikering'] = np.nan

    if node in upstream:
        new_df.loc[i, 'type'] = 'upstream'
        # mind trick: for all upstream ones, prec_node is the node itself!
        new_df.loc[i, 'prec_node'] = node
        new_df.loc[i, 'Brescode'] = np.nan

    elif node in downstream:
        new_df.loc[i, 'type'] = 'downstream'
        new_df.loc[i, 'Brescode'] = np.nan

    elif node in conjunction:
        new_df.loc[i, 'type'] = 'conjunction'
        new_df.loc[i, 'prec_node'] = new_df.loc[i - 1, 'NodeName']
        new_df.loc[i, 'Brescode'] = np.nan

    elif node in bifurcation.keys():

        new_df.loc[i, 'type'] = 'bifurcation'
        new_df.loc[i, 'prec_node'] = new_df.loc[i - 1, 'NodeName']

        # subs node 1 gets the uneven branches
        new_df.loc[i, 'subs_node1'] = bifurcation[node][2]

        # subs node 2 gets the even branches
        new_df.loc[i, 'subs_node2'] = bifurcation[node][1]
        new_df.loc[i, 'Brescode'] = np.nan

    elif node is not upstream + bifurcation.keys() + downstream:

        new_df.loc[i, 'type'] = 'dike'
        new_df.loc[i, 'prec_node'] = new_df.loc[i - 1, 'NodeName']

        # Assign VNK pf:
        pf = VNKpf.loc[new_df.loc[i, 'traject']].values[0]

        if pf == 0:
            new_df.loc[i, 'VNK_pf'] = 1250.0
        else:
            new_df.loc[i, 'VNK_pf'] = pf
    else:
        print('Unknown node type')

# There is one dike node with no brescode: assign one available
new_df.loc[new_df['NodeName'] == '201.907', 'Brescode'] = 50879

new_df = new_df[['NodeName', 'VNKNAAM', 'Brescode', 'traject', 'dikering', 'VNK_pf',
                 'type', 'prec_node', u'subs_node1', u'subs_node2', 'BRANCH',
                 'hground']]

#new_df.to_excel('./data/Rhine_network_nodikecost2.xlsx')
