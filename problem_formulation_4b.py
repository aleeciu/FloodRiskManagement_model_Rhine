# -*- coding: utf-8 -*-
"""
Created on Fri Jun 14 15:26:31 2019

@author: ciullo
"""


from dike_model_function import DikeNetwork  # @UnresolvedImport
from ema_workbench import (Model, TimeSeriesOutcome, ScalarOutcome,
                           RealParameter, IntegerParameter, CategoricalParameter)
from ema_workbench.em_framework.outcomes import Constraint
from itertools import combinations
from functools import partial
import numpy as np

def sum_over(*args):
    return sum(args)

def relative_risk(init_risk, *args):
    # only show relative risk change of each location

    EADa0 = init_risk
    EADa1 = args[0]
    ratio = (EADa0 - EADa1) / EADa0

    return ratio

def risk_diff(init_risk, *args):

    EADa0 = init_risk
    EADa1 = args[0]
    ratio = (EADa0 - EADa1)

    return ratio

def risk_shift_distance(init_risks, *args):
    r1 = (init_risks[0] - args[0]) / (init_risks[0]+0.00001)
    r2 = (init_risks[1] - args[1]) / (init_risks[1]+0.00001)

    distance = abs(r1 - r2) / np.sqrt(2)
    return distance


def max_risk_shift_distance(init_risks, *args):
    distance = []
    for comb in combinations(range(len(init_risks)),2):
        a1, a2 = comb
        # if both areas have initinial risk higher than zero 
        if not (init_risks[a1]*init_risks[a2] == 0): 
            r1 = (init_risks[a1] - args[a1]) / (init_risks[a1])
            r2 = (init_risks[a2] - args[a2]) / (init_risks[a2])

            # distance from the r1 = r2 line
            distance.append(abs(r1 - r2) / np.sqrt(2))
        else:
            distance.append(0)
    return np.max(distance)

#def max_risk_shift_distance(init_risks, *args):
#    
#    distance = []
#    count = 0
#    for arg in args:
#        for i in range(len(init_risks)):
#            if not i == count:
#                r1 = (init_risks[count] - arg) / init_risks[count]
#                r2 = (init_risks[i] - args[i]) / init_risks[i]
#
#                # distance from the r1 = r2 line
#                distance.append(abs(r1 - r2) / np.sqrt(2))
#        count += 1
#    return np.max(distance)

def get_model_for_problem_formulation(problem_formulation_id):
    function = DikeNetwork()
    dike_model = Model('dikesnet', function=function)

    # specify uncertainties, levers
    # uncertainties
    uncert = {'Bmax': [5, 350], 'pfail': [0, 1]}
    cat_uncert = {'Brate': (1., 1.5, 10.)}

    # dike levers
    dike_lev = {'DikeIncrease': [0, 10]}

    # Rfr project code, diversion policies,
    rfr_lev = ['{}_rfr'.format(project_id) for project_id in list(range(0, 156))]
    div_lev = {'Diversion_201.0': [-3, 3],
               'Diversion_401.0': [-3, 3]}

    uncertainties = []
    levers = []
    for dike in function.dikenodes:

        for lev_name in dike_lev.keys():
            name = "{}_{}".format(dike, lev_name)
            levers.append(IntegerParameter(name, dike_lev[lev_name][0],
                                           dike_lev[lev_name][1]))

        for uncert_name in uncert.keys():
            name = "{}_{}".format(dike, uncert_name)
            lower, upper = uncert[uncert_name]
            uncertainties.append(RealParameter(name, lower, upper))

        for uncert_name in cat_uncert.keys():
            name = "{}_{}".format(dike, uncert_name)
            categories = cat_uncert[uncert_name]
            uncertainties.append(CategoricalParameter(name, categories))

    for lev_name in rfr_lev:
            # choose or not a RfR project:
        levers.append(IntegerParameter(lev_name, 0, 1))

    for lev_name in div_lev.keys():
        levers.append(IntegerParameter(lev_name, div_lev[lev_name][0],
                                       div_lev[lev_name][1]))

    dike_model.uncertainties = uncertainties
    dike_model.levers = levers

    direction = ScalarOutcome.MINIMIZE

    # extra outcome of interest
    output_list = ['Qpol', 'Qout', 'Qin', 'wl', 'status', 'critWL']
    eooi = []

    nl_areas = list(range(4))
    de_areas = [4, 5]
    dikerings = np.unique(DikeNetwork().tree.dikering.T)
#    nl_max_costs = 2 * 1e8

    risk_keys = {0: 'minR', 1: 'EAD', 2: 'maxR'}
    ri = risk_keys[1] # risk indicator
    
    # general
    if problem_formulation_id == 0:
        constraints = []
        mins, maxs = 0, 0
        epsilons = []
        
        outcomes = []
        variable_name = []

        [outcomes.append(ScalarOutcome('{}_{}'.format(r, ri))) for r in dikerings]

        for i in eooi:
            [outcomes.append(TimeSeriesOutcome('{}_{}'.format(output_list[i], n))
                             ) for n in function.dikenodes]

        dike_model.outcomes = outcomes

    # non-cooperative
    elif problem_formulation_id == 1:
        constraints = []
        epsilons = []
        outcomes = []
        
        # adapt this constraint
#        variable_name = []
#        [variable_name.append('{}_Dike Inv Cost'.format(a)) for a in nl_areas]
#        variable_name.append('RfR Total Costs')
#        outcomes.append(ScalarOutcome('Investment Costs_nl',
#                                      variable_name=variable_name,
#                                      function=sum_over))
#        constraints.append(Constraint('Investment Costs_nl',
#                                      outcome_names='Investment Costs_nl',
#                                      function=lambda x: max(0, x - nl_max_costs)))

        variable_name = []
        for e in ['{}'.format(ri), 'Dike Inv Cost']:
            variable_name.extend('{}_{}'.format(a, e) for a in nl_areas)
        variable_name.append('RfR Total Costs')

        outcomes.append(ScalarOutcome('Total Costs_nl',
                                      variable_name=variable_name,
                                      function=sum_over, kind=direction))
        epsilons.append(1e7)

        variable_name = []
        for e in ['{}'.format(ri), 'Dike Inv Cost']:
            variable_name.extend('{}_{}'.format(a, e) for a in de_areas)

        outcomes.append(ScalarOutcome('Total Costs_de',
                                      variable_name=variable_name,
                                      function=sum_over, kind=direction))
        epsilons.append(1e7)

        mins = [0.0] * 2
        maxs = [1.9 * 1e10, 2.16 * 1e9]

        dike_model.outcomes = outcomes

    # intra-country cooperation:
    elif problem_formulation_id == 2:
        constraints = []
        epsilons = []
        outcomes = []

#        variable_name = []
#        [variable_name.append('{}_Dike Inv Cost'.format(a)) for a in nl_areas]
#        variable_name.append('RfR Total Costs')
#        outcomes.append(ScalarOutcome('Investment Costs_nl',
#                                      variable_name=variable_name,
#                                      function=sum_over))
#        constraints.append(Constraint('Investment Costs_nl',
#                                      outcome_names='Investment Costs_nl',
#                                      function=lambda x: max(0, x - nl_max_costs)))

        variable_name = []

        for e in ['{}'.format(ri), 'Dike Inv Cost']:
            variable_name.extend('{}_{}'.format(a, e) for a in nl_areas)
        variable_name.append('RfR Total Costs')

        outcomes.append(ScalarOutcome('Total Costs_nl',
                                      variable_name=variable_name,
                                      function=sum_over, kind=direction))
        epsilons.append(1e7)

        variable_name = []
        for e in ['{}'.format(ri), 'Dike Inv Cost']:
            variable_name.extend('{}_{}'.format(a, e) for a in de_areas)

        outcomes.append(ScalarOutcome('Total Costs_de',
                                      variable_name=variable_name,
                                      function=sum_over, kind=direction))
        epsilons.append(1e7)


        init_risks = []
        variable_name = []
        
        for a in nl_areas:
            
            init_risks.append(function.R0['{}_{}'.format(a, ri)].values[0])
            variable_name.append('{}_{}'.format(a, ri))

            outcomes.append(ScalarOutcome('{}_EAD'.format(a),
                                          variable_name=['{}_EAD'.format(a)]))

            # when x is positive, the function gives 0, the constraint is met
            constraints.append(Constraint('{}_EAD'.format(a),
                                          outcome_names='{}_EAD'.format(a),
                                          function=lambda x: max(0, -x)))

        outcomes.append(ScalarOutcome('max_distance_nl',
                        function=partial(max_risk_shift_distance, init_risks),
                        variable_name=variable_name, kind=direction))

#        outcomes.append(ScalarOutcome('{}_u'.format(a),
#                        function=partial(relative_risk, init_risk[0]),
#                                          variable_name=[variable_name[0]]))
#
#        # when x is positive, the function gives 0, the constraint is met
#        constraints.append(Constraint('{}_u'.format(a),
#                                          outcome_names='{}_u'.format(a),
#                                          function=lambda x: max(0, -x)))

        init_risks = []
        variable_name = []

        for a in de_areas:
            
            init_risks.append(function.R0['{}_{}'.format(a, ri)].values[0])
            variable_name.append('{}_{}'.format(a, ri))

            outcomes.append(ScalarOutcome('{}_EAD'.format(a),
                                          variable_name=['{}_EAD'.format(a)]))

            # when x is positive, the function gives 0, the constraint is met
            constraints.append(Constraint('{}_EAD'.format(a),
                                          outcome_names='{}_EAD'.format(a),
                                          function=lambda x: max(0, -x)))

        outcomes.append(ScalarOutcome('max_distance_de',
                        function=partial(max_risk_shift_distance, init_risks),
                        variable_name=variable_name, kind=direction))

#        outcomes.append(ScalarOutcome('{}_u'.format(a),
#                        function=partial(relative_risk, init_risk[0]),
#                                          variable_name=[variable_name[0]]))
#
#        # when x is positive, the function gives 0, the constraint is met
#        constraints.append(Constraint('{}_u'.format(a),
#                                          outcome_names='{}_u'.format(a),
#                                          function=lambda x: max(0, -x)))

        epsilons.extend([0.05] * 5)

        mins = [0.0] * 7
        maxs = [1.9 * 1e10, 2.16 * 1e9] + 2 * [0.708]  # 1/np.sqrt(2)

        dike_model.outcomes = outcomes

    # intra-country cooperation:
    elif problem_formulation_id == '2b':
        constraints = []
        epsilons = []
        outcomes = []

#        variable_name = []
#        [variable_name.append('{}_Dike Inv Cost'.format(a)) for a in nl_areas]
#        variable_name.append('RfR Total Costs')
#        outcomes.append(ScalarOutcome('Investment Costs_nl',
#                                      variable_name=variable_name,
#                                      function=sum_over))
#        constraints.append(Constraint('Investment Costs_nl',
#                                      outcome_names='Investment Costs_nl',
#                                      function=lambda x: max(0, x - nl_max_costs)))

        variable_name = []

        for e in ['{}'.format(ri), 'Dike Inv Cost']:
            variable_name.extend('{}_{}'.format(a, e) for a in nl_areas)
        variable_name.append('RfR Total Costs')

        outcomes.append(ScalarOutcome('Total Costs_nl',
                                      variable_name=variable_name,
                                      function=sum_over, kind=direction))
        epsilons.append(1e7)

        variable_name = []
        for e in ['{}'.format(ri), 'Dike Inv Cost']:
            variable_name.extend('{}_{}'.format(a, e) for a in de_areas)

        outcomes.append(ScalarOutcome('Total Costs_de',
                                      variable_name=variable_name,
                                      function=sum_over, kind=direction))
        epsilons.append(1e7)


        init_risks = function.R0
        variable_name = []
        
        for r in dikerings:
            variable_name.append('{}_{}'.format(r, ri))

            outcomes.append(ScalarOutcome('{}_EADred'.format(r),
                            variable_name=['{}_EAD'.format(r)],
                            function=partial(risk_diff, 
                                             init_risks['{}_EAD'.format(r)].values[0])))

#            outcomes.append(ScalarOutcome('{}_EAD'.format(r),
#                            variable_name=['{}_EAD'.format(r)]))

            # when x is positive, the function gives 0, the constraint is met
            constraints.append(Constraint('{}_EADred'.format(r),
                                          outcome_names='{}_EADred'.format(r),
                                          function=lambda x: max(0, -x)))

        outcomes.append(ScalarOutcome('max_distance',
                        function=partial(max_risk_shift_distance, 
                                         init_risks.values[0]),
                        variable_name=variable_name, kind=direction))

#        outcomes.append(ScalarOutcome('{}_u'.format(a),
#                        function=partial(relative_risk, init_risk[0]),
#                                          variable_name=[variable_name[0]]))
#
#        # when x is positive, the function gives 0, the constraint is met
#        constraints.append(Constraint('{}_u'.format(a),
#                                          outcome_names='{}_u'.format(a),
#                                          function=lambda x: max(0, -x)))

        epsilons.extend([0.05] * 5)

        mins = [0.0] * 7
        maxs = [1.9 * 1e10, 2.16 * 1e9] + 5 * [0.708]  # 1/np.sqrt(2)

        dike_model.outcomes = outcomes

    # inter-country cooperation
    elif problem_formulation_id == 3:
        constraints = []
        epsilons = []
        outcomes = []

#        variable_name = []
#        [variable_name.append('{}_Dike Inv Cost'.format(a)) for a in nl_areas]
#        variable_name.append('RfR Total Costs')
#        outcomes.append(ScalarOutcome('Investment Costs_nl',
#                                      variable_name=variable_name,
#                                      function=sum_over))
#        constraints.append(Constraint('Investment Costs_nl',
#                                      outcome_names='Investment Costs_nl',
#                                      function=lambda x: max(0, x - nl_max_costs)))

        variable_name = []

        for e in ['EAD', 'Dike Inv Cost']:
            variable_name.extend('{}_{}'.format(a, e) for a in nl_areas+de_areas)
        variable_name.append('RfR Total Costs')

        outcomes.append(ScalarOutcome('Total Costs',
                                      variable_name=variable_name,
                                      function=sum_over, kind=direction))
        epsilons.append(1e7)


        init_risks = function.R0
        variable_name = []
        
        for r in dikerings:
            variable_name.append('{}_{}'.format(r, ri))

            outcomes.append(ScalarOutcome('{}_EADred'.format(r),
                            variable_name=['{}_EAD'.format(r)],
                            function=partial(risk_diff, 
                                             init_risks['{}_EAD'.format(r)].values[0])))

#            outcomes.append(ScalarOutcome('{}_EAD'.format(r),
#                            variable_name=['{}_EAD'.format(r)]))

            # when x is positive, the function gives 0, the constraint is met
            constraints.append(Constraint('{}_EADred'.format(r),
                                          outcome_names='{}_EADred'.format(r),
                                          function=lambda x: max(0, -x)))

        outcomes.append(ScalarOutcome('max_distance',
                        function=partial(max_risk_shift_distance, 
                                         init_risks.values[0]),
                        variable_name=variable_name, kind=direction))

        epsilons.extend([0.05] * 5)

        mins = [0.0] * 7
        maxs = [1.9 * 1e10, 2.16 * 1e9] + 5 * [0.708]  # 1/np.sqrt(2)

        dike_model.outcomes = outcomes

    # most-disaggregated version
    elif problem_formulation_id == 'disaggregated':
        constraints = []
        epsilons = []
        outcomes = []

        variable_name = []

        for e in [ri, 'Dike Inv Cost']:
            if e not in ['Dike Inv Cost']:
                for r in dikerings:
                    outcomes.append(ScalarOutcome('{}_{}'.format(r, e)))
                    
            variable_name.extend('{}_{}'.format(a, e) for a in nl_areas)
        variable_name.append('RfR Total Costs')

        outcomes.append(ScalarOutcome('Total Costs_nl',
                                      variable_name=variable_name,
                                      function=sum_over, kind=direction))

        variable_name = []
        for e in ['{}'.format(ri), 'Dike Inv Cost']:    
            variable_name.extend('{}_{}'.format(a, e) for a in de_areas)
        outcomes.append(ScalarOutcome('Total Costs_de',
                                      variable_name=variable_name,
                                      function=sum_over, kind=direction))

#        for e in ['EAD', 'Dike Inv Cost']:
#            for r in dikerings:
#                outcomes.append(ScalarOutcome('{}_{}'.format(r, e)))
#                variable_name.extend(['{}_{}'.format(r, e)])
#        variable_name.append('RfR Total Costs')
#
#        outcomes.append(ScalarOutcome('Total Costs',
#                                      variable_name=variable_name,
#                                      function=sum_over, kind=direction))

#        outcomes.append(ScalarOutcome('RfR Total Costs'))

#        for r in dikerings:
#            outcomes.append(ScalarOutcome('{}_Dike Inv Cost'.format(r)))
#            outcomes.append(ScalarOutcome('{}_EAD'.format(r)))

        init_risks = function.R0.values[0]            
        for comb in combinations(range(len(init_risks)),2):
            a1, a2 = comb
            n1, n2 = dikerings[a1], dikerings[a2]
            outcomes.append(ScalarOutcome('{}_{}_distance'.format(n1,n2),
                        function=partial(risk_shift_distance, 
                                         init_risks[[a1,a2]]),
                        variable_name=['{}_EAD'.format(n1), '{}_EAD'.format(n2)],
                        kind=direction))
                        
#        variable_name = []
#        [variable_name.append('{}_Dike Inv Cost'.format(a)) for a in nl_areas+de-areas]
#        variable_name.append('RfR Total Costs')
#        outcomes.append(ScalarOutcome('Investment Costs_nl',
#                                      variable_name=variable_name,
#                                      function=sum_over))
#        constraints.append(Constraint('Investment Costs_nl',
#                                      outcome_names='Investment Costs_nl',
#                                      function=lambda x: max(0, x - nl_max_costs)))

#        variable_name = []
#
#        for e in ['EAD', 'Dike Inv Cost']:
#            variable_name.extend('{}_{}'.format(a, e) for a in nl_areas+de_areas)
#        variable_name.append('RfR Total Costs')
#
#        outcomes.append(ScalarOutcome('Total Costs',
#                                      variable_name=variable_name,
#                                      function=sum_over, kind=direction))
#
#
#        init_risks = function.R0.values[0]
#        variable_name = []
#        
#        for a in nl_areas+de_areas:
#            variable_name.append('{}_{}'.format(a, ri))
#
#            outcomes.append(ScalarOutcome('{}_EAD'.format(a),
#                                          variable_name=['{}_EAD'.format(a)]))
#
#        outcomes.append(ScalarOutcome('max_distance',
#                        function=partial(max_risk_shift_distance, init_risks),
#                        variable_name=variable_name, kind=direction))

#        outcomes.append(ScalarOutcome('{}_u'.format(a),
#                        function=partial(relative_risk, init_risk[0]),
#                                          variable_name=[variable_name[0]]))
#
#        # when x is positive, the function gives 0, the constraint is met
#        constraints.append(Constraint('{}_u'.format(a),
#                                          outcome_names='{}_u'.format(a),
#                                          function=lambda x: max(0, -x)))

        dike_model.outcomes = outcomes
        mins, maxs = [0,0]

    else:
        raise TypeError('unknonw problem formulation')
    return dike_model, epsilons, (mins, maxs), constraints