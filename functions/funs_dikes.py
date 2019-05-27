# -*- coding: utf-8 -*-
"""
Created on Thu Jul 06 14:51:04 2017

@author: ciullo
"""
import numpy as np

# function splitting the the flow after bifurcation
def bifurcation(inflow,k):
    outflow1 = k*inflow
    outflow2 = inflow - outflow1
    
    return outflow1,outflow2

# function evaluating (1) if dike fails, (2) water volumes involved
#def dikefailure(sb, inflow, hriver, hground, hbas, hbreach, dikelevel ,
#                status_t1, Bmax, dtBmax, simtime, tbreach, critWL, Vinit, Area, 
#                timestep): 
def dikefailure(sb, inflow, hriver, hbas, hground, status_t1, Bmax, Brate, 
                simtime, tbreach, critWL):
        
    # inflow = flow coming into the node
    # status = if False the dike has not failed yet
    # critWL = water level above which we have failure
    
    tbr = tbreach
#    h1 = hriver - hbreach
#    h2 = (hbas + hground) - hbreach

# h river is a water level, hbas a water depth
    h1 = hriver - (hground + hbas)

    
    if status_t1 == True:
        # breach growth
        B = Bmax * (1 - np.exp(-Brate*(simtime - tbreach)))

#        if h2>0 and h1>0 and h1>h2: #h2>0: submerged weir
#            breachflow = ((1 - (h2/h1)**1.5)**0.385)*0.66*np.sqrt(9.81)*B*(h1)*np.sqrt(np.abs(h1))
#            
#            FloodingVolume = breachflow*timestep
#            V = Vinit + FloodingVolume
#
#        elif h2>0 and h1>0 and h2>h1: 
#            breachflow = -((1 - (h1/h2)**1.5)**0.385)*0.66*np.sqrt(9.81)*B*(
#                            h1)*np.sqrt(np.abs(h1))
#            
#            FloodingVolume = breachflow*timestep
#            V = max(0, Vinit + FloodingVolume, (hbreach-hground)*Area)
#
#        elif h2<0 and h1>0:
#            breachflow = max(0, 0.66*np.sqrt(9.81)*B*(h1)*np.sqrt(np.abs(h1)))
#            FloodingVolume = breachflow*timestep
#            V = Vinit + FloodingVolume

        if h1 > 0:
            # see: http://evidence.environment-agency.gov.uk/FCERM/en/FluvialDesignGuide/Chapter7.aspx?pagenum=4#    
            breachflow = 1.7*B*(h1)*1.5
            
        else: # h1 <0; no flow
            breachflow = 0

        outflow = np.max([0, inflow - breachflow])
        
        status_t2 = status_t1
        
    else: # status_t1 == False
        failure = hriver > critWL
        outflow = inflow
        breachflow = 0
        if failure:
            status_t2 = True
            tbr = simtime
        else:
            status_t2 = False
    
    if sb == False:
        outflow = inflow
		      
    return outflow, breachflow, status_t2, tbr

# look up linear function
def Lookuplin(MyFile, inputcol, searchcol, inputvalue):
#    inputvalue = np.asarray([inputvalue])
    minTableValue = np.min(MyFile[:,inputcol])
    maxTableValue = np.max(MyFile[:,inputcol])
    outputvalue = []
    for value in inputvalue:
        
        if value >= maxTableValue: 
            #print 'Overflow in table in procedure inlaatcap. Peil=',inputvalue
            value = maxTableValue-0.01 #decrease the seeking value
        
        if value < minTableValue:
            value = minTableValue+0.01 #increase the seeking value
        
        A = np.max(MyFile[MyFile[:,inputcol]<=value,inputcol])
        B = np.min(MyFile[MyFile[:,inputcol]>value,inputcol]) 
        C = np.max(MyFile[MyFile[:,inputcol]==A,searchcol]) 
        D = np.min(MyFile[MyFile[:,inputcol]==B,searchcol])
    
        outputvalue.extend([C-((D-C)*((value-A)/(A-B)))*1.0])
    return np.asarray(outputvalue)


def init_node(value, time):
        init = np.repeat(value, len(time)) # if not a list, divisions between e.g. cumVol/Area are wrong??!
        return init