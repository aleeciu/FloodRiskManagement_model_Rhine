# -*- coding: utf-8 -*-
"""
Created on Thu Mar 15 11:52:39 2018

@author: ciullo
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


' BUILD FILE WITH A WV21 FC AT EACH LOCATION '

all_frags = pd.read_excel('./data/fragility_curves/FC_fromprojects/WV21.xls')

fragility_curves = pd.DataFrame(columns=['brescode', 'bresnaam', 'mu', 'sd'],
                                index=range(len(all_frags['brescode'].dropna())))

for brescode in all_frags['brescode'].dropna():
    bcode_index = all_frags['brescode'] == brescode

    bname = all_frags.loc[bcode_index, 'bresnaam2'].values[0]

    bname_index = all_frags['Bresnaam'] == bname
    if np.sum(bname_index) != 0:
        mu, sd = (all_frags['mu1'][bname_index].values[0],
                  all_frags['sd1'][bname_index].values[0])

    fragility_curves.loc[bcode_index, ['brescode', 'bresnaam', 'mu', 'sd']
                         ] = brescode, bname, mu, sd

fragility_curves.to_excel('./data/fragility_curves/FC.xlsx')

# note: following WV21, then all FC have the same shape (same std) but diff
# mu, which, however we are going to calibrate anyhow. For the BOA FC, it looks
# like you may get different shapes even for overtopping.
