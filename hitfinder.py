import config as cf
import data_containers as dc
import lar_param as lar

import numpy as np
#import numba as nb

"""numba slows down this code!"""
#@nb.jit(forceobj=True,nopython=True)
def hit_search(data,crp,view,channel,start, dt_min, thr1, thr2):#, pad_l, pad_r):
    """search hit-shape in a list of points"""
    """algorithm from qscan"""
    
    ll = []
    npts = len(data)
    hitFlag = False

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

            #if((singleHit and (h.stop-h.start >= dt_min)) or not singleHit):
            if(h.stop-h.start >= dt_min):
                ll.append(h)

        i+=1
    return ll



def recompute_hit_charge(hit):
    crp, view, ch, start, stop = hit.crp, hit.view, hit.channel, hit.start, hit.stop
    val = 0.
    mean = dc.ped_mean[crp, view, ch]
    for t in range(start, stop):
        val += dc.data[crp, view, ch, t] - mean

    hit.charge = val

def hit_finder(pad_left, pad_right, dt_min, n_sig_1, n_sig_2): 
    
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
                      
            
            adc = dc.data[crp, view, channel, tdc_start:tdc_stop+1]                

            thr1 = dc.ped_mean[crp, view, channel] + n_sig_1 * dc.ped_rms[crp, view, channel]
            thr2 = dc.ped_mean[crp, view, channel] + n_sig_2 * dc.ped_rms[crp, view, channel]

            if(thr1 < 0.5): thr1 = 0.5
            if(thr2 < 0.5): thr2 = 0.5

            hh = hit_search(adc, crp, view, channel, tdc_start, dt_min, thr1, thr2)

            ''' this is wrong '''
            '''
            """add padding to found hits"""
            for i in range(len(hh)): 
                """ to the left """
                if(i == 0): 
                    if(hh[i].start > pad_left):
                        hh[i].start -= pad_left
                    else:
                        hh[i].start = 0
                else:
                    if(hh[i].start - pad_left > hh[i-1].stop):
                        hh[i].start -= pad_left
                    else:
                        hh[i].start = hh[i-1].stop + 1
                

                """ to the right """
                if(i == len(hh)-1):
                    if(hh[i].stop < cf.n_Sample - pad_right):
                        hh[i].stop += pad_right
                    else:
                        hh[i].stop = cf.n_Sample
                else:
                    if(hh[i].stop + pad_right < hh[i+1].start):
                        hh[i].stop += pad_right
                    else:
                        hh[i].stop = hh[i+1].start - 1

            '''



            dc.evt_list[-1].nHits[crp,view] += len(hh)
            dc.hits_list.extend(hh)

    print("nb of hits found ",len(dc.hits_list))

    """ add paddings for hit charge integration properly """
    """ the hits in dc.hits_list are ordered by channel and time """
    
    """ special treatment for the first and last hit of the list """
    h = dc.hits_list[0]
    if(h.start > pad_left):
        h.start -= pad_left
    else:
        h.start = 0            

    h = dc.hits_list[-1]    
    if(h.stop < cf.n_Sample - pad_right):
        h.stop += pad_right
    else:
        h.stop = cf.n_Sample


    for i in range(1, len(dc.hits_list)-1):
        h = dc.hits_list[i]        
        hprev = dc.hits_list[i-1]
        hnext = dc.hits_list[i+1]
        
        """ to the left """
        if(h.view == hprev.view and h.crp == hprev.crp and h.channel == hprev.channel):
            if(h.start - pad_left > hprev.stop):
                h.start -= pad_left
            else:
                h.start = hprev.stop + 1
        else:
            if(h.start > pad_left):
                h.start -= pad_left
            else:
                h.start = 0            

        """ to the right """
        if(h.view == hnext.view and h.crp == hnext.crp and h.channel == hnext.channel):
            if(h.stop + pad_right < hnext.start):
                h.stop += pad_right
            else:
                h.stop = hnext.start - 1
        else:            
            if(h.stop < cf.n_Sample - pad_right):
                h.stop += pad_right
            else:
                h.stop = cf.n_Sample



        

    v = lar.driftVelocity()
    #print("Drift Velocity : v = %.3f mm/mus"%v)

    """ transforms hit channel and tdc to positions """
    [x.hit_positions(v) for x in dc.hits_list]

    """ set hit an index number """
    [dc.hits_list[i].set_index(i) for i in range(len(dc.hits_list))]

    """ compute hit charge in fC """
    [recompute_hit_charge(x) for x in dc.hits_list]
    [x.hit_charge() for x in dc.hits_list]


