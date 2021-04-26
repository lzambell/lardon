import config as cf
import data_containers as dc

import uproot 
import numpy as np
import warnings


class readmc:
    def __init__(self, name):
        self.mc_file = uproot.open(name)
        self.t_pts   = self.mc_file['npts']
        self.n_event = self.t_pts.numentries
        self.n_pts_prev = np.zeros((4,2), dtype=np.uint32)

    def read_event(self, ev):
        if(ev > self.n_event):
            warnings.warn("Warning........... MC FILE HAS ONLY ", self.n_event, ", this event will be empty")
            return


        for icrp in range(2):
            for iv in range(2):
                npts = self.t_pts.array("crp"+str(icrp)+"_v"+str(iv), 
                                        entrystart=ev,
                                        entrystop=ev+1)[0]
                
                p_start = self.n_pts_prev[icrp,iv]
                p_stop  = self.n_pts_prev[icrp,iv]+npts
                
                tree = self.mc_file['CRP'+str(icrp)+'_V'+str(iv)]
      
                
                tbin = tree.array('tbin', entrystart=p_start, entrystop=p_stop)
                ch = tree.array('channel', entrystart=p_start, entrystop=p_stop)
                adc = tree.array('ADC', entrystart=p_start, entrystop=p_stop)

                self.n_pts_prev[icrp, iv] += npts
                
                dc.data[icrp,iv,ch,tbin]+=adc

                '''
                crp, view, tbin, ch, adc = self.data.arrays(['crp', 'view', 'tbin', 'channel', 'ADC'], outputtype=tuple)
        
                view = [1 if x else 0 for x in view]
                
                dc.data[crp,view,ch,tbin]+=adc
                '''
