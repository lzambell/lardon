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
    ev_flag     = False

map_ref = []
evt_list = []
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
