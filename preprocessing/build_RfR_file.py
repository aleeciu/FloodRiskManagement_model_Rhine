# -*- coding: utf-8 -*-
"""
Created on Tue Mar 20 17:22:56 2018

@author: ciullo
"""
import pandas as pd
import numpy as np
from funs_generate_network import get_network


def codify_branch(row, corr):
    if row in corr.keys():
        o = str(corr[row])
    else:
        o = np.nan
    return o


dike_nodes = get_network(Musk_params=True,
                         fragcurves=False,
                         Losses=True)[3]

log_string = ['index == "{}" |'.format(NodeName)
              for NodeName in dike_nodes[:-1]]
#log_string = []
logical = ' '.join(log_string + ['index == "{}"'.format(dike_nodes[-1])])


branches = ['BR', 'PK', 'WA', 'BOME', 'IJ', 'NR']
codes = [101, 301, 201, 201, 501, 401]
corr = dict(zip(branches, codes))

rfr_projects = pd.DataFrame(index=dike_nodes)

'''Load RfR project'''
files_name = [
    'Bovenrijn',
    'PannerdenschKanaal',
    'Nederrrijn_Lek',
    'Waal',
    'IJssel']
# count number of prpjects you the profile for:
count = 0.0

for fl in files_name:
    projects = pd.read_excel('./data/measures/RfR_projects_data/{}.xls'.format(fl),
                             sheetname=None)
    for sheet in projects.keys():
        count += 1

        dummy_df = projects[sheet][['locatie', 'verschil', 'riviertak']]

        # skip projects with non-float values: I could not convert the
        # their object type into floats
        if dummy_df.locatie.dtype != 'float64':
            continue

        # I want rounded integer from the lacatie (km) field:
        dummy_df.loc[:, 'locatie'] = dummy_df.locatie.round().astype(int)
        # Corresponding branch codes:
        dummy_df.loc[:, 'branch'] = dummy_df.apply(
            lambda row: codify_branch(row.riviertak, corr=corr),
            axis=1)
        # Drop rows not belonging to the study area:
        dummy_df = dummy_df.loc[dummy_df[['branch']].dropna().index]

        # Build NodeNames as km.branch:
        dummy_df.loc[:, 'NodeName'] = dummy_df.apply(
            lambda row: '{}.{}'.format(row.branch, row.locatie),
            axis=1)
        # Sum up wl reduction for the same node. If you want something
        # that works per stretch, then you need a list of nodes instead
        # of a single node
        dummy_df = dummy_df.groupby(['NodeName']).sum()[['verschil']]
        # Rename col with the project name:
        dummy_df = dummy_df.rename(columns={'verschil': sheet})

        # Select only the nodes of interest:
        dummy_df = dummy_df.query(logical, engine='python')

        rfr_projects = pd.concat([rfr_projects, dummy_df], axis=1)

rfr_projects.replace(0, np.nan)
# Load costs:
projects_costs = pd.read_excel('./data/measures/RfR_projects_data/Costs.xlsx',
                               sheetname=1)[['code maatregel', 'm2', 'Investerings kosten']]

projects_costs['code maatregel'] = projects_costs['code maatregel'].astype(str)

for col in rfr_projects.columns:
    cost = projects_costs.loc[projects_costs['code maatregel'] == col,
                              'Investerings kosten'].values
    rfr_projects.loc['Costs_1e6', col] = cost

# Drop measures with no costs:
rfr_projects = rfr_projects.loc[:, rfr_projects.loc['Costs_1e6'] != 0]

# number of dropped plausible projects becuase of line 42-43:
n = 1 - len(rfr_projects.columns) / float(count)  # about 15 %

rfr_projects.to_excel('./data/measures/RfR.xlsx')
