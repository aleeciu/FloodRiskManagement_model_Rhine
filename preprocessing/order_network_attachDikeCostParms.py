# -*- coding: utf-8 -*-
"""
Created on Thu Apr 26 10:36:32 2018

@author: ciullo
"""

''' Attach cost parameters to nodes '''

import pandas as pd


location = 'C:\\Users\\ciullo\\OneDrive - Stichting Deltares\\Desktop\\main_D\\RhineModel_py3'

rhine_file = pd.read_excel(
    '{}/data/Rhine_network_nodikecost.xlsx'.format(location),
    dtype=object)
costs_param = pd.read_excel('{}/data/measures/WV21_costsparams.xlsx'.format(location),
                            sheet_name='data_restructured')

costs_param['ratio1'] = costs_param['pL1'] / costs_param['L1']
costs_param['ratio2'] = costs_param['pL2'] / costs_param['L2']
costs_param['ratio3'] = costs_param['pL3'] / costs_param['L3']

costs_param = costs_param.drop(['pL1', 'L1', 'pL2', 'L2', 'pL3', 'L3'], axis=1)
#rhine_file = rhine_file.drop(['c', 'b', 'lambda', 'ratio'], axis = 1)

rhine_file = rhine_file.reindex(
    columns=rhine_file.columns.append(costs_param.columns[1:]))

columns = costs_param.columns[1:]

for i in costs_param.CODE:
    attach = costs_param.loc[costs_param.CODE == i, columns].fillna(0)
    rhine_file.loc[rhine_file['Brescode'] == i, columns] = attach.values  # ???

# TODO: For some reasons, there is a mismatch between the VW21 and the costs
# trajects. No info at rows 8, 25.
rhine_file.loc[8, columns] = rhine_file.loc[7, columns]
rhine_file.loc[25, columns] = rhine_file.loc[26, columns]

rhine_file.to_excel('{}/data/Rhine_network.xlsx'.format(location))
