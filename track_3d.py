import config as cf
import data_containers as dc
import numpy as np
import math
from scipy.interpolate import CubicSpline


def complete_trajectory(track, other, view):
    
    #reversed because cubic spline wants a increasing x only
    x_o = [x[0] for x in reversed(other.path)]
    z_o = [x[1] for x in reversed(other.path)]


    """ order lists according to z increasing """ 
    z_o, x_o = (list(t) for t in zip(*sorted(zip(z_o, x_o))))



    """ ugly fix to remove multiple hits at the same z """
    if(len(z_o) != len(set(z_o))) :
        i=1
        while(i < len(z_o)-1):
            if(z_o[i-1] >= z_o[i]):
                del z_o[i]
                del x_o[i]
            else:
                i += 1


    """at least 3 points for the spline """
    if(len(z_o) < 4):
        return -1, [], []


    spline = CubicSpline(z_o, x_o)
    deriv = spline.derivative()

    a0, a1 = 0., 0.
    dx, dy, dz = 0., 0., 0.

    trajectory = []
    dQds       = []
    length = 0.
    
    for i in range(len(track.path)):
        x = track.path[i][0]
        z = track.path[i][1]

        if( i == 0 ):
            a0 = 0. if track.ini_slope == 0 else 1./track.ini_slope
        else:
            dx = track.path[i][0] - track.path[i-1][0]
            dy = 0.
            dz = track.path[i][1] - track.path[i-1][1]
            
            a0 = 0. if dx == 0 else dz/dx

        y = float(spline(z))
        a1 = float(deriv(z))

        a1 = 0. if a1 == 0 else 1./a1
        
        dr = cf.ChanPitch
        if(a1 == 0):
            dr *= math.sqrt(1. + pow(a0,2))
        else : 
            dr *= math.sqrt(1. + pow(a0, 2)*(1./(pow(a1,2) + 1.)))

        length += dr
        dQds.append(track.dQ[i]/dr)

        if(view == 0):
            trajectory.append( (x, y, z) )
        else:
            trajectory.append( (y, x, z) )


    return length, trajectory, dQds


def find_tracks(ztol, qfrac):
    t_v0 = [x for x in dc.tracks2D_list if x.view == 0]
    t_v1 = [x for x in dc.tracks2D_list if x.view == 1]
    
    for ti in t_v0:

        if(ti.matched >= 0): 
            continue

        tbest = ti
        mindist = 9999.
        match = False
        for tj in t_v1:

            if(tj.matched >= 0):
                continue

            if( (ti.ini_crp == tj.ini_crp) and (ti.end_crp == tj.end_crp) ):

                d_start = math.fabs( ti.path[0][1] - tj.path[0][1] )
                d_stop  = math.fabs( ti.path[-1][1] - tj.path[-1][1] )

                qv0 = ti.tot_charge
                qv1 = tj.tot_charge
                balance =  math.fabs(qv0 - qv1)/(qv0 + qv1)

                if( d_start < ztol and d_stop < ztol and balance < qfrac):
                    d = d_start + d_stop
                    if(d < mindist):
                        tbest = tj
                        mindist = d
                        match = True

        if(match == True):            
            track = dc.trk3D(ti, tbest)

            l, t, q = complete_trajectory(ti, tbest, 0)
            

            if(l < 0):
                continue
            track.set_view0(l, t, q)
            
            l, t, q = complete_trajectory(tbest, ti, 1)

            
            if(l < 0):
                continue
            track.set_view1(l, t, q)
            track.matched(ti, tbest)

            dc.tracks3D_list.append(track)
            dc.evt_list[-1].nTracks3D += 1
