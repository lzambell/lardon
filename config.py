data_path = "/eos/experiment/neutplatform/protodune/rawdata/np02/rawdata/"
calib_path = "/afs/cern.ch/user/n/np02onlp/public/calib/pedestals/"

n_CRP = 4
n_CRPUsed = 2
n_View = 2
n_Sample = 10000
n_ChanPerCRP = 960
n_ChanTot = 7680 # = n_CRP * n_View * n_ChanPerCRP
n_ChanPerView = 3840 #=n_CRP * n_CharPerCRP
n_Sampling = 0.4 #in mu-seconds
ChanPitch = 0.3125 #cm
LAr_Temperature = 87.5
E_drift = 0.166 #in kV/cm


""" BROKEN CHANNELS TO BE REMOVED FROM THE ANALYSIS"""

""" provided in DAQ Channel """
daq_broken_channels = [] #3508, 3507, 3505, 3504]
""" or provided in (crp, view, channel) tuple """
crp_broken_channels = [(0,0,x) for x in range(946, 953)]

run_inv_signal = 1256
