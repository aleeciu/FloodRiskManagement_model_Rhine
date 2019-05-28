# -*- coding: utf-8 -*-
"""
Created on Wed Apr 24 11:18:24 2019

@author: ciullo
"""

from ema_workbench import (MultiprocessingEvaluator,
                           Scenario, SequentialEvaluator, Policy)

from ema_workbench.em_framework.samplers import sample_uncertainties, sample_levers
from ema_workbench.util import ema_logging
import pickle
import time
import pandas as pd
from problem_formulation_2 import get_model_for_problem_formulation


def build_Pol(df, pf):
    policies = []
    for i, row in df.iterrows():
        name = '{}_{}'.format(pf, str(i))
        decision = {lever: row[lever] for lever in df.columns}
        policies.append(Policy(name=name, **decision))
    return policies

if __name__ == '__main__':
    ema_logging.log_to_stderr(ema_logging.INFO)

    pf = 0

    dike_model = get_model_for_problem_formulation(pf)[0]
        
    average_values = {'Bmax': 150, 'Brate': 1.5, 'pfail': 0.5}
    scen1 = {}

    for key in dike_model.uncertainties:
        dikename, unc = key.name.split('_')
        scen1.update({key.name: average_values[unc]})

    ref_scenario = Scenario('reference', **scen1)
    policies = sample_levers(dike_model, 500)

    # No dike increase, none of the rfr project is selected, no damage
    # reduction:
    dike_increase = dict(zip((91, 101, 201, 301, 401, 501),
                             (0, 0, 0, 0, 0, 0)))

    names = {'DikeIncrease': 0, 'rfr': 0,
             '201.0': 0, '401.0': 0}
    ## zero policy:
    pol0 = {}

    for key in dike_model.levers:
        s1, s2 = key.name.split('_')
        if s2 == 'DikeIncrease':
            branch = s1.split('.')[0]
            pol0.update({key.name: dike_increase[int(branch)]})

        else:
            pol0.update({key.name: names[s2]})

    policy0 = Policy('Policy 0', **pol0)

##  set of policies
#    policies = []
#    for i in range(11):
#         dike_increase = dict(zip((101, 201, 301, 401, 501),
#                                  (i, 0, 0, i, 0)))
#         pol0 = {}
#
#         for key in dike_model.levers:
#             s1, s2 = key.name.split('_')
#             if s2 == 'DikeIncrease':
#                branch = s1.split('.')[0]
#                pol0.update({key.name: dike_increase[int(branch)]})
#
#             else:
#                pol0.update({key.name: names[s2]})
#
#         policies.append(Policy('Policy {}'.format(i), **pol0))
#
#    n_scenarios = 5
#    scenarios = sample_uncertainties(dike_model, n_scenarios)
#    n_policies = 1000

# pf1
#    df = pd.read_excel('{}/FRM_model_Rhine_RESULTS/2nd paper/1bn_nl/PF1.xlsx'.format(folder)
#                       )[[l.name for l in dike_model.levers]]
#    policies = build_Pol(df, 'pf1')

## pf1c
#    df = pd.read_excel('{}/FRM_model_Rhine_RESULTS/2nd paper/1bn_nl/PF1b.xlsx'.format(folder)
#                       )[[l.name for l in dike_model.levers]]
#    policies.extend(build_Pol(df, 'pf1c'))
#### zero policy
#    policies.extend(build_Pol(df.iloc[[0]]*0, '0'))
##
### pf2
#    df = pd.read_excel('{}/FRM_model_Rhine_RESULTS/2nd paper/1bn_nl/PF2.xlsx'.format(folder)
#                       )[[l.name for l in dike_model.levers]]
#    policies.extend(build_Pol(df, 'pf2'))
### pf2a
#    df = pd.read_excel('{}/FRM_model_Rhine_RESULTS/2nd paper/1bn_nl/PF2a.xlsx'.format(folder)
#                       )[[l.name for l in dike_model.levers]]
#    policies.extend(build_Pol(df, 'pf2a'))

## SIMULATION:
# single run
    start = time.time()
    dike_model.run_model(ref_scenario, policy0)
    end = time.time() - start
    print(end)
    results = dike_model.outcomes_output

##     series run
#    with SequentialEvaluator(dike_model) as evaluator:
#        results = evaluator.perform_experiments(ref_scenario, policies)

# multiprocessing
#    with MultiprocessingEvaluator(dike_model, n_processes = 3) as evaluator:
#        results = evaluator.perform_experiments(ref_scenario, policies)
###
#### save results
#    with open("./results//FirstDecentResults/simulation_noncooppolsx100.p".format(pf), "wb") as f:
#            pickle.dump(results, f)
