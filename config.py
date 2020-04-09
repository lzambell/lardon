class pdmap: 
    view = -1
    crp    = -1
    vchan  = -1
    ped    = 0.0
    rms    = 0.0


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

    def reset(self):
        self.ini_crp = -1
        self.end_crp = -1
        self.view    = -1
        self.ini_slope   = -1
        self.ini_slope_err = -1
        self.end_slope   = -1
        self.end_slope_err = -1
        self.nHits   = 0
        self.path    = []
        self.dQ      = []
        self.chi2    = -1
        
data_path = "/eos/experiment/neutplatform/protodune/rawdata/np02/rawdata/"
calib_path = "/afs/cern.ch/user/n/np02onlp/public/calib/pedestals/"

map_ref = []
evt_list = []
hits_list = []
tracks2D_list = []

n_CRP = 4
n_View = 2
n_Sample = 10000
n_ChanPerCRP = 960
n_ChanTot = 7680 # = n_CRP * n_View * n_ChanPerCRP
n_ChanPerView = 3840 #=n_CRP * n_CharPerCRP
n_Sampling = 0.4 #in mu-seconds
ChanPitch = 0.3125 #cm
LAr_Temperature = 87.5
E_drift = 0.166 #in kV/cm


# in DAQ Channel
# should be turned into a function depending on run nb
# shall add lem borders ? 
broken_channels = [3508, 3507, 3505, 3504]

run_inv_signal = 1256
