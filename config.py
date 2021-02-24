data_path = "/eos/experiment/neutplatform/protodune/rawdata/np02/rawdata/"
calib_path = "/afs/cern.ch/user/n/np02onlp/public/calib/pedestals/"

store_path = "."
plot_path  = "."

n_CRP = 4
n_CRPUsed = 4
n_View = 2
n_Sample = 10000
n_ChanPerCRP = 960
n_ChanTot = 7680 # = n_CRP * n_View * n_ChanPerCRP
n_ChanPerView = 3840 #=n_CRP * n_CharPerCRP
n_Sampling = 0.4 #in mu-seconds
ChanPitch = 0.3125 #cm
LAr_Temperature = 87.5
E_drift = 0.166 #kV/cm
Anode_Z = 300. #cm
len_det_x = 600. #cm
len_det_y = 600. #cm


ADCtofC = 35.64 #from Qscan, to be cross checked in units of (ADC x us) / fC

""" BROKEN CHANNELS TO BE REMOVED FROM THE ANALYSIS"""

""" provided in DAQ Channel """
daq_broken_channels = [] #3508, 3507, 3505, 3504]

crp_broken_channels = []

if(n_CRPUsed == 4):
    crp_broken_channels += [(2,0,x) for x in range(0,960)]
    crp_broken_channels += [(2,1,x) for x in range(0,960)]
    crp_broken_channels += [(3,0,x) for x in range(320, 960)]
    crp_broken_channels += [(3,1,x) for x in range(0, 640)]

""" or provided in (crp, view, channel) tuple """
""" run 1250 """
#crp_broken_channels += [(1,0,x) for x in range(255,320)]

""" run 1323 """
#crp_broken_channels += [(0,0,x) for x in range(944, 960)]

""" run 1415 """
crp_broken_channels += [(1,0,x) for x in range(32, 42)]
crp_broken_channels += [(1,0,x) for x in range(136, 142)]
crp_broken_channels += [(1,1,x) for x in range(422, 426)]
crp_broken_channels += [(3,0,x) for x in range(255, 261)]

""" run 1417 """
#crp_broken_channels += [(1,0,x) for x in range(32, 63)]
#crp_broken_channels += [(1,0,x) for x in range(128, 159)]
#crp_broken_channels += [(1,1,x) for x in range(416, 447)]

""" run 1407 """
#crp_broken_channels += [(1,1,425)]
#crp_broken_channels += [(3,0,259)]



run_inv_signal = 1256
