# -*- coding: utf-8 -*-
"""
Created on Thu Mar 22 09:13:13 2018

@author: ciullo
"""

from dike_model_function import DikeNetwork  # @UnresolvedImport
from ema_workbench import (Model, TimeSeriesOutcome, ScalarOutcome,
                           RealParameter, IntegerParameter, CategoricalParameter)
from ema_workbench.em_framework.outcomes import Constraint
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


def absolute_risk(init_risk, *args):
    # only show relative risk change of each location

    EADa0 = init_risk
    EADa1 = args[0]
    ratio = (EADa0 - EADa1)

    return ratio


def risk_shift_distance(init_risk, *args):

    EADa1 = args[0]
    EADb1 = np.sum(args[1:])

    EADa0, EADb0 = init_risk

    r1 = (EADa0 - EADa1) / EADa0
    r2 = (EADb0 - EADb1) / EADb0

    # distance from the r1 = r2 line
    distance = abs(r1 - r2) / np.sqrt(2)

    return distance


def absolute_risk_shift_distance(init_risk, *args):

    EADa1 = args[0]
    EADb1 = np.sum(args[1:])

    EADa0, EADb0 = init_risk

    r1 = (EADa0 - EADa1) / (EADa0 + EADb0)
    r2 = (EADb0 - EADb1) / (EADa0 + EADb0)

    # distance from the r1 = r2 line
    distance = abs(r1 - r2) / np.sqrt(2)

    return distance
    
# Note on upperlim (maxs) of the risk tranfers criteria:
# no risk increase is allowed => 0<r<1 (i.e. when no risk reduction, r=0; when
# total risk reduction, r=1), thus 0<distance<1/sqrt(2) (=0.708)

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

    nl_max_costs = 2 * 1e8

    # general
    if problem_formulation_id == 0:
        mins, maxs = 0, 0
        constraints = []
        epsilons = []
        areas = nl_areas + de_areas

        outcomes = []
        [outcomes.append(ScalarOutcome('{}_EAD'.format(a))) for a in areas]

        variable_name = []
        [variable_name.append('{}_Dike Inv Cost'.format(a)) for a in nl_areas]
        variable_name.append('RfR Total Costs')

        outcomes.append(ScalarOutcome('Investment Costs_nl',
                                      variable_name=variable_name,
                                      function=sum_over))

        variable_name = []
        [variable_name.append('{}_Dike Inv Cost'.format(a)) for a in de_areas]

        outcomes.append(ScalarOutcome('Investment Costs_de',
                                      variable_name=variable_name,
                                      function=sum_over))

        for i in eooi:
            [outcomes.append(TimeSeriesOutcome('{}_{}'.format(output_list[i], n))
                             ) for n in function.dikenodes]

        dike_model.outcomes = outcomes

    # utilitarian
    elif problem_formulation_id == 1:
        constraints = []
        epsilons = []
        outcomes = []

        variable_name = []
        [variable_name.append('{}_Dike Inv Cost'.format(a)) for a in nl_areas]
        variable_name.append('RfR Total Costs')

        outcomes.append(ScalarOutcome('Investment Costs_nl',
                                      variable_name=variable_name,
                                      function=sum_over))
        constraints.append(Constraint('Investment Costs_nl',
                                      outcome_names='Investment Costs_nl',
                                      function=lambda x: max(0, x - nl_max_costs)))

        variable_name = []
        for e in ['EAD', 'Dike Inv Cost']:
            variable_name.extend('{}_{}'.format(a, e) for a in nl_areas)
        variable_name.append('RfR Total Costs')

        outcomes.append(ScalarOutcome('Total Costs_nl',
                                      variable_name=variable_name,
                                      function=sum_over, kind=direction))
        epsilons.append(1e7)

        variable_name = []
        for e in ['EAD', 'Dike Inv Cost']:
            variable_name.extend('{}_{}'.format(a, e) for a in de_areas)

        outcomes.append(ScalarOutcome('Total Costs_de',
                                      variable_name=variable_name,
                                      function=sum_over, kind=direction))
        epsilons.append(1e7)

        mins = [0.0] * 2
        maxs = [1.88 * 1e10, 1.71 * 1e9]

        dike_model.outcomes = outcomes

    # constrained utilitarian
    elif problem_formulation_id == '1c':
        constraints = []
        epsilons = []
        outcomes = []

        variable_name = []
        [variable_name.append('{}_Dike Inv Cost'.format(a)) for a in nl_areas]
        variable_name.append('RfR Total Costs')

        outcomes.append(ScalarOutcome('Investment Costs_nl',
                                      variable_name=variable_name,
                                      function=sum_over))
        constraints.append(Constraint('Investment Costs_nl',
                                      outcome_names='Investment Costs_nl',
                                      function=lambda x: max(0, x - nl_max_costs)))

        variable_name = []
        for e in ['EAD', 'Dike Inv Cost']:
            variable_name.extend('{}_{}'.format(a, e) for a in nl_areas)
        variable_name.append('RfR Total Costs')

        outcomes.append(ScalarOutcome('Total Costs_nl',
                                      variable_name=variable_name,
                                      function=sum_over, kind=direction))
        epsilons.append(1e7)

        variable_name = []
        for e in ['EAD', 'Dike Inv Cost']:
            variable_name.extend('{}_{}'.format(a, e) for a in de_areas)

        outcomes.append(ScalarOutcome('Total Costs_de',
                                      variable_name=variable_name,
                                      function=sum_over, kind=direction))
        epsilons.append(1e7)

        for a in nl_areas:
            init_risk = [function.R0['{}_EAD'.format(a)].values[0],
                         np.sum([function.R0['{}_EAD'.format(i)].values[0]
                                 for i in nl_areas if i != a])]

            variable_name = ['{}_EAD'.format(a)] + ['{}_EAD'.format(i
                                                                    ) for i in nl_areas if i != a]

            outcomes.append(ScalarOutcome('{}_u'.format(a),
                                          function=partial(
                                              relative_risk, init_risk[0]),
                                          variable_name=[variable_name[0]]))

            # when x is positive, the function gives 0, the constraint is met
            constraints.append(Constraint('{}_u'.format(a),
                                          outcome_names='{}_u'.format(a),
                                          function=lambda x: max(0, -x)))

        # only take one criteria in germany, cause there are only two locations.
        # A to B or B to A is the same
        for a in de_areas:

            init_risk = [function.R0['{}_EAD'.format(a)].values[0],
                         np.sum([function.R0['{}_EAD'.format(i)].values[0]
                                 for i in de_areas if i != a])]

            variable_name = ['{}_EAD'.format(a)] + ['{}_EAD'.format(i
                                                                    ) for i in de_areas if i != a]

            outcomes.append(ScalarOutcome('{}_u'.format(a),
                                          function=partial(
                                              relative_risk, init_risk[0]),
                                          variable_name=[variable_name[0]]))

            # when x is positive, the function gives 0, the constraint is met
            constraints.append(Constraint('{}_u'.format(a),
                                          outcome_names='{}_u'.format(a),
                                          function=lambda x: max(0, -x)))

        mins = [0.0] * 2
        maxs = [1.88 * 1e10, 1.71 * 1e9]

        dike_model.outcomes = outcomes

    # utilitarian + relative risk trnafers
    elif problem_formulation_id == 2:
        constraints = []
        epsilons = []
        outcomes = []

        variable_name = []
        [variable_name.append('{}_Dike Inv Cost'.format(a)) for a in nl_areas]
        variable_name.append('RfR Total Costs')

        outcomes.append(ScalarOutcome('Investment Costs_nl',
                                      variable_name=variable_name,
                                      function=sum_over))
        constraints.append(Constraint('Investment Costs_nl',
                                      outcome_names='Investment Costs_nl',
                                      function=lambda x: max(0, x - nl_max_costs)))

        variable_name = []

        for e in ['EAD', 'Dike Inv Cost']:
            variable_name.extend('{}_{}'.format(a, e) for a in nl_areas)
        variable_name.append('RfR Total Costs')

        outcomes.append(ScalarOutcome('Total Costs_nl',
                                      variable_name=variable_name,
                                      function=sum_over, kind=direction))
        epsilons.append(1e8)

        variable_name = []
        for e in ['EAD', 'Dike Inv Cost']:
            variable_name.extend('{}_{}'.format(a, e) for a in de_areas)

        outcomes.append(ScalarOutcome('Total Costs_de',
                                      variable_name=variable_name,
                                      function=sum_over, kind=direction))
        epsilons.append(1e8)

        for a in nl_areas:
            init_risk = [function.R0['{}_EAD'.format(a)].values[0],
                         np.sum([function.R0['{}_EAD'.format(i)].values[0]
                                 for i in nl_areas if i != a])]

            variable_name = ['{}_EAD'.format(a)] + ['{}_EAD'.format(i
                                                                    ) for i in nl_areas if i != a]

            outcomes.append(ScalarOutcome('{}_d'.format(a),
                                          function=partial(risk_shift_distance,
                                                           init_risk),
                                          variable_name=variable_name, kind=direction))

            outcomes.append(ScalarOutcome('{}_u'.format(a),
                                          function=partial(
                                              relative_risk, init_risk[0]),
                                          variable_name=[variable_name[0]]))

            # when x is positive, the function gives 0, the constraint is met
            constraints.append(Constraint('{}_u'.format(a),
                                          outcome_names='{}_u'.format(a),
                                          function=lambda x: max(0, -x)))

        # only take one criteria in germany, cause there are only two locations.
        # A to B or B to A is the same
        for a in de_areas:
            init_risk = [function.R0['{}_EAD'.format(a)].values[0],
                         np.sum([function.R0['{}_EAD'.format(i)].values[0]
                                 for i in de_areas if i != a])]

            variable_name = ['{}_EAD'.format(a)] + ['{}_EAD'.format(i
                                                                    ) for i in de_areas if i != a]

            if a == 5:
                outcomes.append(ScalarOutcome('{}_d'.format(a),
                                              function=partial(
                                                  risk_shift_distance, init_risk),
                                              variable_name=variable_name, kind=direction))

            outcomes.append(ScalarOutcome('{}_u'.format(a),
                                          function=partial(
                                              relative_risk, init_risk[0]),
                                          variable_name=[variable_name[0]]))

            # when x is positive, the function gives 0, the constraint is met
            constraints.append(Constraint('{}_u'.format(a),
                                          outcome_names='{}_u'.format(a),
                                          function=lambda x: max(0, -x)))

        epsilons.extend([0.05] * 5)

        mins = [0.0] * 7
        maxs = [1.88 * 1e10, 1.71 * 1e9] + 5 * [0.708]  # 1/np.sqrt(2)

        dike_model.outcomes = outcomes

    # utilitarian + absolute risk trnafers
    elif problem_formulation_id == '2a':
        constraints = []
        epsilons = []
        outcomes = []

        variable_name = []
        [variable_name.append('{}_Dike Inv Cost'.format(a)) for a in nl_areas]
        variable_name.append('RfR Total Costs')

        outcomes.append(ScalarOutcome('Investment Costs_nl',
                                      variable_name=variable_name,
                                      function=sum_over))
        constraints.append(Constraint('Investment Costs_nl',
                                      outcome_names='Investment Costs_nl',
                                      function=lambda x: max(0, x - nl_max_costs)))

        variable_name = []
        for e in ['EAD', 'Dike Inv Cost']:
            variable_name.extend('{}_{}'.format(a, e) for a in nl_areas)
        variable_name.append('RfR Total Costs')

        outcomes.append(ScalarOutcome('Total Costs_nl',
                                      variable_name=variable_name,
                                      function=sum_over, kind=direction))
        epsilons.append(1e8)

        variable_name = []
        for e in ['EAD', 'Dike Inv Cost']:
            variable_name.extend('{}_{}'.format(a, e) for a in de_areas)

        outcomes.append(ScalarOutcome('Total Costs_de',
                                      variable_name=variable_name,
                                      function=sum_over, kind=direction))
        epsilons.append(1e8)

        for a in nl_areas:
            init_risk = [function.R0['{}_EAD'.format(a)].values[0],
                         np.sum([function.R0['{}_EAD'.format(i)].values[0]
                                 for i in nl_areas if i != a])]

            variable_name = ['{}_EAD'.format(a)] + ['{}_EAD'.format(i
                                                                    ) for i in nl_areas if i != a]

            outcomes.append(ScalarOutcome('{}_d'.format(a),
                                          function=partial(absolute_risk_shift_distance,
                                                           init_risk),
                                          variable_name=variable_name, kind=direction))

            outcomes.append(ScalarOutcome('{}_u'.format(a),
                                          function=partial(
                                              absolute_risk, init_risk[0]),
                                          variable_name=[variable_name[0]]))

            # when x is positive, the function gives 0, the constraint is met
            constraints.append(Constraint('{}_u'.format(a),
                                          outcome_names='{}_u'.format(a),
                                          function=lambda x: max(0, -x)))

        # only take one criteria in germany, cause there are only two locations.
        # A to B or B to A is the same
        for a in de_areas:
            init_risk = [function.R0['{}_EAD'.format(a)].values[0],
                         np.sum([function.R0['{}_EAD'.format(i)].values[0]
                                 for i in de_areas if i != a])]

            variable_name = ['{}_EAD'.format(a)] + ['{}_EAD'.format(i
                                                                    ) for i in de_areas if i != a]

            if a == 5:
                outcomes.append(ScalarOutcome('{}_d'.format(a),
                                              function=partial(
                    absolute_risk_shift_distance, init_risk),
                    variable_name=variable_name, kind=direction))

            outcomes.append(ScalarOutcome('{}_u'.format(a),
                                          function=partial(
                                              absolute_risk, init_risk[0]),
                                          variable_name=[variable_name[0]]))

            # when x is positive, the function gives 0, the constraint is met
            constraints.append(Constraint('{}_u'.format(a),
                                          outcome_names='{}_u'.format(a),
                                          function=lambda x: max(0, -x)))

        epsilons.extend([0.05] * 5)

        mins = [0.0] * 7
        maxs = [1.88 * 1e10, 1.71 * 1e9] + 5 * [0.708]

        dike_model.outcomes = outcomes
        
    else:
        raise TypeError('unknonw identifier')
    return dike_model, epsilons, (mins, maxs), constraints
