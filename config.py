data_path = "/eos/experiment/neutplatform/protodune/rawdata/np02/rawdata/"
calib_path = "/afs/cern.ch/user/n/np02onlp/public/calib/pedestals/"

""" User's specific """
store_path = ""
plot_path  = ""
default_reco = "default_reco.yaml"

experiment="np02"

""" GENERAL PARAMETERS FOR NP02 RECONSTRUCTION """
n_CRP = 4
n_CRPUsed = -1 #overwritten by runs.yaml
n_View = 2
n_Sample = 10000

n_ChanPerCRP = 960
n_ChanTot = 7680
n_ChanPerView = 3840

n_Sampling = 0.4 #in mu-seconds
ChanPitch = 0.3125 #cm

LAr_Temperature = 87. #K

E_drift = 0.166 #kV/cm

Anode_Z = 300. #cm
len_det_x = 600. #cm
len_det_y = 600. #cm

run_inv_signal = 1256

ADCtofC = 35.64 # in units of (ADC x us) / fC

""" CRP gaps in cm"""
crp_01_x = 2.5
crp_01_y = -1.
crp_03_x = 0.
crp_03_y = 2.5
