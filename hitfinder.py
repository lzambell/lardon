import config as cf
import numpy as np
import numba as nb

@nb.jit(forceobj=True,nopython=True)
def HitSearch(data, ROI, rms, view, crp, chan):
    test_hits = []
    dt_min = 10 #nb of time bins


    ROI_ind = np.nonzero(ROI[1:] != ROI[:-1])[0] + 1
    ROIval = np.split(data, ROI_ind)
    ROIval = ROIval[0::2] if ROI[0] else ROIval[1::2]

    for gpe in ROIval:
        npts = len(gpe)
        if(npts < dt_min): continue
        i = 0
        hitFlag = False
        while(i < npts):
            val = gpe[i]
            if(hitFlag == False):
                h = cf.hits(view,crp,chan,i,i,val,val,i)
                hitFlag = True
                i+=1
            else:
                if(val > h.max_adc):
                    h.max_t = i
                    h.max_adc = val
                h.charge += val
                h.stop = i
                i+=1

        if(h.stop-h.start > dt_min):
            test_hits.append(h)#cf.hits(view, crp, ichan, start, stop, s

    """
    i = 0
    threshold = 4.*rms
    hitFalg = False
    while (i < cf.n_Sample):
        val = data[i]            
        if(val > 0):
                if(hitFlag is False):
                    hitFlag = True
                    h = cf.hits(view,crp,chan,i,i,val,val,i)
                    i+=1
                while(i < cf.n_Sample and val > 0):
                    if(val > h.max_adc):
                        h.max_t = i
                        h.max_adc = val
                    h.charge += val
                    h.stop = i
                    i+=1
                if(h.stop-h.start > dt_min):
                    test_hits.append(h)#cf.hits(view, crp, ichan, start, stop, sum_adc, max_t, max_v))
            else:
               i+=1
    """
    return test_hits



def HitFinder(data, ROI, rms):
    for icrp in range(2):
        for iview in range(cf.n_View):
            for ichan in range(cf.n_ChanPerCRP):
                if(np.sum(ROI) == 0): continue
                hh = HitSearch(data[icrp,iview,ichan], ROI[icrp,iview,ichan], rms[icrp,iview,ichan], iview, icrp, ichan)
                cf.hits_list.extend(hh)
            print("view ", iview, " crp ", icrp, " -> ", len(cf.hits_list))
