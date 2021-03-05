import config as cf

import numpy as np
import math

def driftVelocity(E=cf.E_drift):
    """ for E > 0.5 kV/cm : 
    Walkowiak NIM A 449 (2000) 288
    for E < 0.5 kV/cm :
    ICARUS NIM A 516 (2004) 68
    """


    T = cf.LAr_Temperature
    #E = cf.E_drift


    T_walk = 90.371 #K
    walk = [-0.01481, -0.0075, 0.141, 12.4, 1.627, 0.317]
    icarus = [-0.03229, 6.231, -10.62, 12.74, -9.112, 2.83]

    vd = 0.
    dt = T-T_walk

    if(E > 0.5):
        vd = (walk[0]*dt+1.) * (walk[2]*E*math.log(1.+walk[3]/E)+walk[4]*pow(E,walk[5]))+walk[1]*dt

    else:
        for p in range(len(icarus)):
            vd += icarus[p]*pow(E,p)
        Etmp = 0.5
        tmp1 = 0.
        for p in range(len(icarus)):            
            tmp1 += icarus[p]*pow(Etmp,p)
        tmp2 = (walk[0]*dt+1.) * (walk[2]*Etmp*math.log(1.+walk[3]/Etmp)+walk[4]*pow(Etmp,walk[5]))+walk[1]*dt
        
        vd = vd * tmp2/tmp1

    return vd # in mm/us


def recombination(E=cf.E_drift):
    """ following ICARUS NIM A 523 (2004) 275 """
    A = 0.8
    k = 0.0486
    rho = 1.39
    alpha = 0.83
    dEdx = 2.11 #for muon at mip

    return A/(1.+(k/E)*dEdx/rho)
