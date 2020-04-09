import config as cf
import numpy as np
import time

map_ref = []
evt_list = []
hits_list = []
tracks2D_list = []


data = np.zeros((cf.n_CRPUsed, cf.n_View, cf.n_ChanPerCRP, cf.n_Sample)) #crp, view, vchan

"""
the mask will be used to differentiate background (True for noise processing) from signal (False for noise processing)
at first everything is considered background (all at True)
"""
mask = np.ones((cf.n_CRPUsed, cf.n_View, cf.n_ChanPerCRP, cf.n_Sample), dtype=bool)

"""
alive_chan mask intends to not take into account broken channels
True : not broken
False : broken
"""

alive_chan = np.ones((cf.n_CRPUsed,cf.n_View, cf.n_ChanPerCRP, cf.n_Sample), dtype=bool)

ncluster = np.zeros((cf.n_CRPUsed,cf.n_View),dtype=int)


def reset_event():
    data[:,:,:,:] = 0.
    mask[:,:,:,:] = True
    ncluster[:,:] = 0
    hits_list.clear()
    tracks2D_list.clear()
    


class pdmap: 
    crp    = -1
    view   = -1
    vchan  = -1
    ped    = 0.0
    rms    = 0.0
    def __init__(self, crp, view, vchan):
        self.crp   = crp
        self.view  = view
        self.vchan = vchan
        self.ped   = 0.
        self.rms   = 0.

    def set_pedestal(self, ped, rms):
        self.ped = ped
        self.rms = rms


    def get_ana_chan(self):
        return self.crp, self.view, self.vchan

class event:
    run_nb      = -1
    evt_nb_loc  = -1 #in that subfile
    evt_nb_glob = -1 #in the run
    time_s      = 0
    time_ns     = 0
    evt_flag     = False

    def __init__(self, run_nb, evt_glob, t_s, t_ns, flag):
        self.run_nb      = run_nb
        self.evt_nb_glob = evt_glob
        self.evt_nb_loc  = -1
        self.time_s      = t_s
        self.time_ns     = t_ns
        self.evt_flag    = flag
        
    def __eq__(self, other):
        return (self.run_nb, self.evt_nb_glob, self.evt_nb_loc, self.time_s, self.time_ns, self.evt_flag) == (other.run_nb, other.evt_nb_glob, other.evt_nb_loc, other.time_s, other.time_ns, other.evt_flag)

    def dump(self):
        print("RUN ",self.run_nb, " EVENT ", self.evt_nb_loc, " / ", self.evt_nb_glob,)
        print("Taken at ", time.ctime(self.time_s), " + ", self.time_ns, " ns ")




class hits:
    view    = -1
    crp     = -1
    channel = -1
    start   = -1
    stop    = -1
    charge  = -1
    max_t   = -1 
    max_adc = -1
    cluster = -1
    X       = -1
    Z       = -1
    def __init__(self, crp, view, channel, start, stop, charge, max_t, max_adc):
        self.crp     = crp
        self.view    = view
        self.channel = channel
        self.start   = start
        self.stop    = stop
        self.charge  = charge
        self.max_t   = max_t
        self.max_adc = max_adc
        self.cluster = -1 #cluster
        self.X       = -1
        self.Z       = -1

    def __lt__(self,other):
        #"""sort hits by increasing channel and increasing Z"""
        #return (self.X < other.X) or (self.X== other.X and self.Z < other.Z)
        """ sort hits by decreasing Z and increasing channel """
        return (self.Z > other.Z) or (self.Z == self.Z and self.X < other.X)
    def GetDistances(self, v, pitch):
        self.X = self.channel*pitch
        
        if(self.crp == 1 and self.view == 0):
            self.X -= 300
        self.Z = 300. - self.max_t*0.4*v*0.1

class trk2D:
    ini_crp = -1
    end_crp = -1
    view    = -1
    ini_slope   = -1
    ini_slope_err = -1
    end_slope   = -1
    end_slope_err = -1
    nHits   = -1
    path    = []
    dQ      = []
    chi2    = -1
    def __init__(self, ini_crp, view, ini_slope, ini_slope_err, x0, y0, q0, chi2):
        self.ini_crp = ini_crp
        self.end_crp = ini_crp
        self.view    = view
        self.ini_slope   = ini_slope
        self.ini_slope_err   = ini_slope_err
        self.end_slope   = ini_slope
        self.end_slope_err   = ini_slope_err
        self.nHits   = 1
        self.path    = [(x0,y0)]
        self.dQ      = [q0]
        self.chi2    = chi2

    def add_hit(self, slope, slope_err, x, y, q, chi2):
        self.end_slope = slope
        self.end_slope_err = slope_err
        self.nHits += 1
        self.path.append([x,y])
        self.dQ.append(q)
        self.chi2 = chi2
        
