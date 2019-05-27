# -*- coding: utf-8 -*-
"""
Created on Tue Nov 21 23:11:30 2017

@author: ciullo
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

data_net = pd.read_excel('./data/Rhine_network2.xlsx', dtype=object)

files = list('0ABCDE')  # BovnRijn, pnKan, Waal, Lek, IJssel

branches = [91, 101, 201, 301, 401, 501]

link = dict(zip(branches, files))

#plt.figure(figsize = (15,15))
for b in branches:
    fl = link[b]
    name = '{}_branch'.format(fl)
    wls = pd.read_excel(
        './data/hydraulics/SOBEK_water_levels/{}'.format('{}_wl.xlsx'.format(name))).iloc[:, 1:]
    Qs = pd.read_excel(
        './data/hydraulics/SOBEK_discharges/{}'.format('{}_Q.xlsx'.format(name))).iloc[:, 1:]

    Qs.columns, wls.columns = 2 * \
        [data_net['NodeName'][data_net['BRANCH'] == b]]

    for l in Qs.columns:
        # choose ascending limb only
        idmax = Qs[l].idxmax()

        # get the mean value between the two limbs,
        # since they are of different lenghts, interpolation of the limmbs is required.
        # a: acending, d: discending
        min_a, max_a = min(Qs.loc[:idmax, l]), max(Qs.loc[:idmax, l])
        min_d, max_d = min(Qs.loc[idmax:, l]), max(Qs.loc[idmax:, l])

        new_a_x = np.linspace(min_a, max_a, 500)
        new_d_x = np.linspace(min_d, max_d, 500)
        # degree five polynom
        a_coeff = np.polyfit(Qs.loc[:idmax, l].values,
                             wls.loc[:idmax, l].values, 5)
        d_coeff = np.polyfit(Qs.loc[idmax:, l].values,
                             wls.loc[idmax:, l].values, 5)

        new_a_y = np.polyval(a_coeff, new_a_x)
        new_d_y = np.polyval(d_coeff, new_d_x)

        mid_y = [np.mean([new_a_y[i], new_d_y[i]]) for i in range(500)]
        mid_x = [np.mean([new_a_x[i], new_d_x[i]]) for i in range(500)]

        # now interpolate-extrapolate (linearly) the mid line, extrapolation
        # up to very high Q values is needed in case, due to diversion policies,
        # Q higher than those of SOBEK rating curves are experienced.
        # This is particularly useful for the IJssel, which seems getting
        # much more water than usual.
        # alternative: run SOBEK with higher upstream discharge, eg 40000:

        # extrapolate based only on the higher discharges:
        m_coeff = np.polyfit(mid_x[len(mid_x) / 2:], mid_y[len(mid_x) / 2:], 1)
        f_m = np.poly1d(m_coeff)
        extr_range = np.linspace(mid_x[-5], 4 * mid_x[-1], 200)

#            plt.figure(figsize = (5, 5))
#            plt.plot(Qs[l], wls[l],
#                     label = 'Sobek Q-stage relationship', color = 'black')
#            plt.plot(mid_x, mid_y, 'bo', label = 'mid rating curve')
#            plt.plot(extr_range, f_m(extr_range), 'ro',
#                     label = 'extrap rating curve')
#            plt.show()

        Q = np.concatenate([mid_x[:-5], extr_range])
        wl = np.concatenate([mid_y[:-5], f_m(extr_range)])

        new_rating_curve = np.column_stack([Q, wl])

#            plt.plot(Q, wl, label = 'new rating curve', color = 'g', linewidth = 1)
##            plt.xlim(mid_x[-10], mid_x[-1])
##            plt.ylim(mid_y[-10], mid_y[-1])
#
##            plt.plot(Qs.loc[:idmax, l], wls.loc[:idmax, l], label = 'Ascending limb')
#
#            title = 'Location {}'.format(l)
#            plt.title(title)
#            plt.legend()
#            plt.show()
#
        filename = './data/hydraulics/rating_curves_extrap/{}_ratingcurve.txt'.format(
            l)
        np.savetxt(filename, new_rating_curve)

#            figname = './data/rating_curves/{}_ratingcurve.png'.format(dike = dike)
#            plt.savefig(figname)
