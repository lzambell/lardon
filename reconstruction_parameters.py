import yaml

import config as cf
import channelmapper as cm

class Reco:
    def __init__(self, fparam):    
        with open(cf.default_reco) as f: 
            self.param = yaml.load(f, Loader=yaml.FullLoader)

            if(fparam != cf.default_reco):
                with open(fparam) as f: 
                    custom = yaml.load(f, Loader=yaml.FullLoader)
                
                for key, val in custom.items():
                    self.param[key] = val 


class Run:
    def __init__(self, run):
        """ defaut """
        self.n_CRPUsed = cf.n_CRP
        self.daq_broken_channels = []

        with open(cf.experiment+"/runs.yaml") as f:
            val = yaml.load(f, Loader=yaml.FullLoader)
            if(run in val):
                self.n_CRPUsed = val[run]['nCRPUsed']
                if(len(val[run]['daq_broken_channels']) == 1 and val[run]['daq_broken_channels'][0] == -1):
                    print("NO DAQ BROKEN CHANNELS SET FOR RUN ", run)
                    print(" ----->>> NONE IS ASSUMED")
                    print(" !!!! PLEASE CHECK AND UPDATE ACCORDINGLY !!!!")
                    print(" (if no channels are broken change [-1] to [])")
                else:
                    self.daq_broken_channels = val[run]['daq_broken_channels']
            else:
                print("NO RUN ", run, " EXISTS IN THE ", cf.experiment, " CONFIGURATION FILE ")
                print(" ----->>> DEFAULT IS ASSUMED")
                print(" !!!! PLEASE CHECK AND UPDATE ACCORDINGLY !!!!")
                
        if(self.n_CRPUsed == 4):
            """ remove all CRP2 v0 """
            self.daq_broken_channels += [cm.CRPToDAQ(2, 0, x) for x in range(0, 960)]

            """ remove all CRP2 v1 """
            self.daq_broken_channels += [cm.CRPToDAQ(2, 1, x) for x in range(0, 960)]

            """ remove un-instrumented CRP3 v0 """
            self.daq_broken_channels += [cm.CRPToDAQ(3, 0, x) for x in range(320, 960)]            

            """ remove un-instrumented CRP3 v1 """
            self.daq_broken_channels += [cm.CRPToDAQ(3, 1, x) for x in range(0, 640)]            
