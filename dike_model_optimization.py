# -*- coding: utf-8 -*-
"""
Created on Wed Apr 24 11:18:24 2019

@author: ciullo
"""

from ema_workbench import (MultiprocessingEvaluator, ema_logging,
                           Scenario, SequentialEvaluator)

from ema_workbench.em_framework.optimization import EpsilonProgress, HyperVolume
import matplotlib.pyplot as plt
from problem_formulation import get_model_for_problem_formulation

import time


if __name__ == '__main__':
    ema_logging.log_to_stderr(ema_logging.INFO)

    pf = '2a'

    problem = get_model_for_problem_formulation(pf)

    dike_model = problem[0]
    epsilons = problem[1]
    mins, maxs = problem[2]
    constraints = problem[3]

    average_values = {'Bmax': 150, 'Brate': 1.5, 'pfail': 0.5}
    scen1 = {}

    for key in dike_model.uncertainties:
        dikename, unc = key.name.split('_')
        scen1.update({key.name: average_values[unc]})

    ref_scenario = Scenario('reference', **scen1)

    nfe = 200000

    convergence_metrics = [HyperVolume(minimum=mins, maximum=maxs),
                           EpsilonProgress()]
    # OPTIMIZATION:
#    with SequentialEvaluator(dike_model) as evaluator:
#        results, convergence = evaluator.optimize(nfe = nfe,
#                                                 searchover = 'levers',
#                                                 epsilons = epsilons,
#                                                 constraints = constraints,
#                                                 convergence = convergence_metrics,
#                                                 reference = ref_scenario
#                                                 )

    start = time.time()
    with MultiprocessingEvaluator(dike_model, n_processes = 23) as evaluator:
         results, convergence = evaluator.optimize(
             nfe=nfe,
             searchover='levers',
             epsilons=epsilons,
             constraints=constraints,
             convergence=convergence_metrics,
             reference=ref_scenario
         )
    end = time.time() - start
    print(end)

    results.to_excel("./results/results_pf{}_{}nfe.xlsx".format(pf, nfe))
    convergence.to_excel("./results//convergence_pf{}_{}nfe.xlsx".format(pf, nfe))
    fig, axes = plt.subplots(1, 2)
    axes[0].plot(convergence.nfe, convergence.epsilon_progress)
    axes[1].plot(convergence.nfe, convergence.hypervolume)

    plt.savefig("./results//convergence_pf{}_{}nfe.png".format(pf, nfe))
