import config as cf
import data_containers as dc
import lar_param as lar

import numpy as np
#import numba as nb

"""numba slows down this code!"""
#@nb.jit(forceobj=True,nopython=True)
def HitSearch(data,rms,crp,view,channel,start):
    """search hit-shape in a list of points"""
    """algorithm from qscan"""
    
    ll = []
    npts = len(data)
    hitFlag = False

    """reco parameters"""
    dt_min = 10    
    thr1 = 3.*rms
    thr2 = 4.*rms

    if(thr1 < 1): thr1 = 1.
    if(thr2 < 1): thr2 = 1.

    i=0
    hitFlag = False
    minimum = 10000
    minSamp = -1
    singleHit = True



    while(i<npts):
        while(i < npts and data[i] >= thr1):
            val = data[i]        
            it = i+start

            if(hitFlag == False):
                hitFlag = True
                singleHit = True
                
                h = dc.hits(crp,view,channel,it,0,0.,it,val)
                minSamp = -1
                
            if(it > h.max_t and val < h.max_adc - thr2 and (minSamp==-1 or minimum >= val)):
                minSamp = it
                minimum = val

                
            if(minSamp >= 0 and it > minSamp and val > minimum + thr2 and (it-h.start) >= dt_min):
                h.stop = minSamp-1
                ll.append(h)
                hitFlag = True
                singleHit = False
                h = dc.hits(crp,view,channel,minSamp,0,0,it,val)
                minSamp = -1

                
            if(h.stop == 0 and val > h.max_adc):
                h.max_t = it
                h.max_adc = val
                if(minSamp >= 0):
                    minSamp = -1
                    minimum = val
                    
            h.charge += val
            i+=1
        if(hitFlag == True):
            hitFlag = False
            h.stop = it-1

            if((singleHit and (h.stop-h.start >= dt_min)) or not singleHit):
                ll.append(h)

        i+=1
    return ll


def HitFinder(rms): 
    
    dt_min = 10
    pad_left = 5
    pad_right = 10

    """ get boolean roi based on mask and alive channels """
    ROI = np.array(~dc.mask & dc.alive_chan, dtype=bool)

    """ adds 0 (False) and the start and end of each waveform """
    falses = np.zeros((cf.n_CRPUsed, cf.n_View, cf.n_ChanPerCRP,1),dtype=int)
    ROIs = np.r_['-1',falses,np.asarray(ROI,dtype=int),falses]
    d = np.diff(ROIs)

    """ a change from false to true in difference is = 1 """
    start = np.where(d==1)
    """ a change from true to false in difference is = -1 """
    end   = np.where(d==-1)
    """ look at long enough sequences of trues """
    gpe = (end[3]-start[3])>=dt_min

    assert len(start[0])==len(end[0]), " Mismatch in groups of hits"
    assert len(gpe)==len(start[0]), "Mismatch in groups of hits"
    
    ntr = 0
    ndble = 0
    for g in range(len(gpe)):
        if(gpe[g]):
            ntr += 1

            """ make sure starts and ends of hit group are in the same crp,view,channel """
            for i in range(3):
                assert start[i][g] == end[i][g], "Hit Mismatch"

            crp = start[0][g]
            view = start[1][g]
            channel = start[2][g]

            tdc_start = start[3][g]
            tdc_stop = end[3][g]
            
            """ add l/r paddings """
            for il in range(pad_left, 0, -1):
                if(tdc_start-1>=0 and not ROI[crp,view,channel,tdc_start-1]):
                    tdc_start -= 1
                else:
                    break

            for ir in range(0, pad_right):
                if(tdc_stop+1 < cf.n_Sample and not ROI[crp,view,channel,tdc_stop+1]):
                    tdc_stop += 1
                else:
                    break
                      
            
            adc = dc.data[crp,view,channel,tdc_start:tdc_stop+1]                
            hh = HitSearch(adc, rms[crp,view,channel], crp, view, channel, tdc_start)
            
            dc.hits_list.extend(hh)

    print("nb of hits found ",len(dc.hits_list))


    v = lar.driftVelocity()
    print("Drift Velocity : v = %.3f mm/mus"%v)

    """ transforms hit channel and tdc to positions """
    [x.GetDistances(v, cf.ChanPitch) for x in dc.hits_list]

