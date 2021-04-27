#from ROOT import TTree, TFile
import ROOT as root
import data_containers as dc
import config as cf
import numpy as np
from array import array

class store_root:
    def __init__(self, name):
        self.froot = root.TFile(name, 'recreate')
        self.tree  = root.TTree('hit_map', 'tree')
        
        self.event = array('i', [0])
        self.tree.Branch('event', self.event, 'event/I')
        self.trackNB = array('i', [0])

        self.tree.Branch('trackNB', self.trackNB, 'trackNB/I')

        self.channel = root.std.vector('int')()
        self.t_start = root.std.vector('int')()
        self.t_max   = root.std.vector('int')()
        self.t_stop  = root.std.vector('int')()
        self.crp     = root.std.vector('int')()
        self.view    = root.std.vector('int')()
        self.adc_max = root.std.vector('float')()
        self.charge  = root.std.vector('float')()
        self.pos     = root.std.vector('float')()
        self.drift   = root.std.vector('float')()
        
        self.tree.Branch('channel', self.channel)
        self.tree.Branch('t_max', self.t_max)
        self.tree.Branch('t_start', self.t_start) 
        self.tree.Branch('t_stop', self.t_stop) 
        self.tree.Branch('crp', self.crp) 
        self.tree.Branch('view', self.view) 
        self.tree.Branch('adc_max', self.adc_max) 
        self.tree.Branch('charge', self.charge) 
        self.tree.Branch('pos', self.pos) 
        self.tree.Branch('drift', self.drift) 
        
        self.to_larsoft_crp  = [3, 1, 0, 2]
        self.to_larsoft_view = [1, 0]

    def clear(self):
        self.channel.clear()
        self.t_max.clear()
        self.t_start.clear()
        self.t_stop.clear()
        self.crp.clear()
        self.view.clear()
        self.adc_max.clear()
        self.charge.clear()
        self.pos.clear()
        self.drift.clear()

    def store_and_close(self):
        self.tree.Write()
        self.froot.Close()


        

    def store_found_hits(self):
        self.event[0] = int(dc.evt_list[-1].evt_nb_glob)
        self.clear()
        for hv in dc.hits_list:
            self.channel.push_back(int(hv.channel))
            self.t_max.push_back(int(hv.max_t))
            self.t_start.push_back(int(hv.start))
            self.t_stop.push_back(int(hv.stop))
            crp, view = self.to_larsoft_crp[hv.crp], self.to_larsoft_view[hv.view]
            self.crp.push_back(int(crp))
            self.view.push_back(int(view))
            
            self.adc_max.push_back(float(hv.max_adc))
            self.charge.push_back(float(hv.charge))

            if(view == 0):
                self.pos.push_back(float(-1.*hv.X))
            else:
                self.pos.push_back(float(hv.X+300.))
            self.drift.push_back(float(hv.Z))
        self.tree.Fill()            

    def store_hits_tracks3D(self):
        self.event[0] = int(dc.evt_list[-1].evt_nb_glob)

        for it in range(len(dc.tracks3D_list)):
            t = dc.tracks3D_list[it]
            self.clear()
            self.trackNB[0] = int(it)

            print("track ", it)

            for t2d in dc.tracks2D_list:
                if(t2d.matched == it):
                    trkid = t2d.trackID
                    view = t2d.view

                    print(" found a matching 2D track ", trkid)

                    for hv in dc.hits_list:
                        if(np.fabs(hv.matched) == trkid):
                            self.channel.push_back(int(hv.channel))
                            self.t_max.push_back(int(hv.max_t))
                            self.t_start.push_back(int(hv.start))
                            self.t_stop.push_back(int(hv.stop))
                            self.crp.push_back(int(hv.crp))
                            self.view.push_back(int(hv.view))
                            self.adc_max.push_back(float(hv.max_adc))
                            self.charge.push_back(float(hv.charge))
                            self.pos.push_back(float(hv.X))
                            self.drift.push_back(float(hv.Z))
            self.tree.Fill()
