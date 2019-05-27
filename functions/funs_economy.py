# -*- coding: utf-8 -*-
"""
Created on Sun Oct 29 15:52:12 2017

@author: ciullo
"""
import numpy as np

# cost of raising a dike

def cost_dike_fun(c, b, lambd, ratio, dikeincrease, dikeinit=0):
    dikeincrease = dikeincrease * 100  # cm
    dikeinit = dikeinit * 100
    
    cost = ((c + b * dikeincrease) * np.exp(lambd * (dikeinit + dikeincrease)))
    return cost * 1e6 * ratio


def cost_dike(c1, b1, lambd1, ratio1,
              c2, b2, lambd2, ratio2,
              c3, b3, lambd3, ratio3,
              dikeincrease):
    costs = sum([cost_dike_fun(c1, b1, lambd1, ratio1, dikeincrease),
                 cost_dike_fun(c2, b2, lambd2, ratio2, dikeincrease),
                 cost_dike_fun(c3, b3, lambd3, ratio3, dikeincrease)
                 ])
    return costs

# discount overall the planning period n (years)


def discount(amount, rate, n):
    factor = 1 + rate / 100
    disc_amount = amount * 1 / (np.repeat(factor, n)**(range(1, n + 1)))
    return disc_amount


def cost_evacuation(N_evacuated, days_to_threat):
    # if days to threat is zero, then no evacuation happens, costs are zero
    cost = N_evacuated * 22 * (days_to_threat + 3) * (int(days_to_threat > 0))
    return cost
