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
    def __init__(self, view, crp, channel, start, stop, charge, max_t, max_adc,cluster):
        self.view    = view
        self.crp     = crp
        self.channel = channel
        self.start   = start
        self.stop    = stop
        self.charge  = charge
        self.max_t   = max_t
        self.max_adc = max_adc
        self.cluster = cluster

data_path = "/eos/experiment/neutplatform/protodune/rawdata/np02/rawdata/"
calib_path = "/afs/cern.ch/user/n/np02onlp/public/calib/pedestals/"

map_ref = []
evt_list = []
hits_list = []

n_CRP = 4
n_View = 2
n_Sample = 10000
n_ChanPerCRP = 960
n_ChanTot = 7680 # = n_CRP * n_View * n_ChanPerCRP
n_ChanPerView = 3840 #=n_CRP * n_CharPerCRP
n_Sampling = 0.4 #in mu-seconds


# in DAQ Channel
# should be turned into a function depending on run nb
# shall add lem borders ? 
broken_channels = [3508, 3507, 3505, 3504]

run_inv_signal = 1256
