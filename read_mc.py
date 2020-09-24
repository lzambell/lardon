import config as cf
import data_containers as dc

import uproot 
import numpy as np



class readmc:
    def __init__(self, name):
        self.mc_file = uproot.open(name)
        self.data    = self.mc_file['mcdata']
        
    def read_event(self, ev):
        crp, view, tbin, ch, adc = self.data.arrays(['crp', 'view', 'tbin', 'channel', 'ADC'], outputtype=tuple)
        
        view = [1 if x else 0 for x in view]

        dc.data[crp,view,ch,tbin]+=adc
