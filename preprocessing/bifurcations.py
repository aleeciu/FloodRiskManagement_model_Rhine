# -*- coding: utf-8 -*-
"""
Created on Wed Mar 14 15:42:21 2018

@author: ciullo
"""
from scipy.optimize import curve_fit
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

'''Branching function'''
" Rijn-PK-Waal "
PK_Waal = pd.read_excel(
    './data/hydraulics/SOBEK_discharges/1st_bifrc.xlsx').iloc[1:]
PK_Waal = PK_Waal.groupby((PK_Waal.index / 24).astype(int)).mean()
PK_Waal.columns = ['BR', 'PK', 'WA']

PK_Waal.plot(title='first bifurcation point', figsize=(8, 6))
plt.xlabel('days')
plt.ylabel('discharge [m3/s]')
# plt.savefig('./data/hydraulics/SOBEK_discharges/discharge_1stbif.png')

PK_Waal['to_theWaal'] = PK_Waal['WA'] / PK_Waal['BR']

PK_Waal.plot(kind='scatter', y='to_theWaal', x='BR',
             title='first bifurcation point', figsize=(8, 6))
plt.xlabel('discharge BR [m3/s]')
plt.ylabel('discharge fraction Waal/BR')

z_1 = np.polyfit(PK_Waal['BR'],
                 PK_Waal['to_theWaal'], 4)

r_1 = np.poly1d(z_1)
Q_1 = np.linspace(np.min(PK_Waal['BR']), np.max(PK_Waal['BR']), 30)

plt.plot(PK_Waal['BR'], PK_Waal['to_theWaal'], 'o')
plt.plot(Q_1, r_1(Q_1))
# plt.savefig('./data/hydraulics/SOBEK_discharges/dischargefraction_1stbif.png')

"PK-LEK-IJssel"
NK_IJ = pd.read_excel(
    './data/hydraulics/SOBEK_discharges/2nd_bifrc.xlsx').iloc[1:]
NK_IJ = NK_IJ.groupby((NK_IJ.index / 24).astype(int)).mean()
NK_IJ.columns = ['PK', 'Lek', 'IJ']

NK_IJ.plot(title='second bifurcation point', figsize=(8, 6))
plt.xlabel('days')
plt.ylabel('discharge [m3/s]')
# plt.savefig('./data/hydraulics/SOBEK_discharges/discharge_2ndbif.png')

NK_IJ['to_theLek'] = NK_IJ['Lek'] / NK_IJ['PK']

NK_IJ.plot(kind='scatter', y='to_theLek', x='PK',
           title='first bifurcation point', figsize=(8, 6))
plt.xlabel('discharge PK [m3/s]')
plt.ylabel('discharge fraction Lek/PK')

z_2 = np.polyfit(NK_IJ['PK'], NK_IJ['to_theLek'], 3)
r_2 = np.poly1d(z_2)
Q_2 = np.linspace(np.min(NK_IJ['PK']), np.max(NK_IJ['PK']), 30)

plt.plot(NK_IJ['PK'], NK_IJ['to_theLek'], 'o')
plt.plot(Q_2, r_2(Q_2))
# plt.savefig('./data/hydraulics/SOBEK_discharges/dischargefraction_2ndbif.png')

# table = pd.DataFrame([Q_1, Q_2,
#                      r_1(Q_1), 1-r_1(Q_1),
#                      r_2(Q_2), 1-r_2(Q_2)],
#             index = [['{}_Q'.format(bifname) for bifname in ['101.1000', '301.1000'
#                      ]] + ['{}_ratio'.format(bifname) for bifname in ['201.0', '301.0'
#                      ]] + ['{}_ratio'.format(bifname) for bifname in ['401.0', '501.0']]]
#                      ).T.to_excel('./data/hydraulics/distr_factors.xlsx')

plt.figure(figsize=(8, 6))
for i in np.linspace(-2, 2, 5) / 10.0:
    r = r_1(Q_1)
    r += r * i
    plt.plot(Q_1, r, label='change of {}0 %'.format(i))
#plt.plot(PK_Waal['BR'], PK_Waal['to_theWaal'], 'o')
plt.legend()
plt.xlabel('discharge BR [m3/s]')
plt.ylabel('discharge fraction Waal/BR')
# plt.savefig('./data/hydraulics/SOBEK_discharges/fasciodischargefraction_1stbif.png')


plt.figure(figsize=(8, 6))
for i in np.linspace(-2, 2, 5) / 10.0:
    r = r_2(Q_2)
    r += r * i
    plt.plot(Q_2, r, label='change of {}0 %'.format(i))
#plt.plot(PK_Waal['BR'], PK_Waal['to_theWaal'], 'o')
plt.legend()
plt.xlabel('discharge PK [m3/s]')
plt.ylabel('discharge fraction Lek/PK')
# plt.savefig('./data/hydraulics/SOBEK_discharges/fasciodischargefraction_2ndbif.png')


plt.figure(figsize=(8, 6))
# to the PK
i = -4 / 10.0
r = r_1(Q_1)
r += r * i
plt.plot(Q_1, (1 - r) * Q_1, label='change of {}0 %'.format(i))

plt.show()
plt.figure(figsize=(8, 6))

# to the lek
i = -8 / 10.0
r = r_2((1 - r) * Q_1)
r += r * i
plt.plot((1 - r) * Q_1, r * (1 - r) * Q_1, label='change of {}0 %'.format(i))
plt.show()
