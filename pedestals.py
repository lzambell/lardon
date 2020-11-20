import config as cf
import data_containers as dc
import glob 

import numba as nb
import numpy as np

@nb.jit(nopython = True)
def compute_pedestal_RMS_nb(data, mask):
    """ do not remove the @jit above, the code would then take ~40 s to run """
    shape = data.shape

    """ to avoid cases where the rms goes to 0"""
    min_val = 1e-5

    res   = np.zeros(shape[:-1])
    for idx,v in np.ndenumerate(res):
        ch = data[idx]
        ma  = mask[idx]
        """ use the assumed mean method """
        K = ch[0]
        n, Ex, Ex2, = 0., 0., 0.
        for x,v in zip(ch,ma):
            if(v == True):
                n += 1
                Ex += x-K
                Ex2 += (x-K)*(x-K)

        """cut at min. 10 pts to compute a proper RMS"""
        if( n < 10 ):
            res[idx] = -1.
        else:
            val = np.sqrt((Ex2 - (Ex*Ex)/n)/(n-1))
            res[idx] = min_val if val < min_val else val
    return res



def compute_pedestal_RMS():
    """ the numpy way is slower and cannot handle well dead channels """
    #dc.ped_rms =  np.std(dc.data[dc.mask], axis=3)
    dc.ped_rms = compute_pedestal_RMS_nb(dc.data, dc.mask)

def compute_pedestal_mean():
    """ the numpy way may be faster but do not handle dead channels """
    #dc.ped_mean = np.mean(dc.data[dc.mask], axis=3)

    with np.errstate(divide='ignore', invalid='ignore'):
        dc.ped_mean = np.einsum('ijkl,ijkl->ijk', dc.data, dc.mask)/dc.mask.sum(axis=3)
        """require at least 3 points to take into account the mean"""
        dc.ped_mean[dc.mask.sum(axis=3) < 3] = -999.




def get_reference_ped_mean(i):
    return dc.map_ped[i].ref_ped

def get_reference_ped_RMS(i):
    return dc.map_ped[i].ref_rms

def get_closest_ped_run(run):
    calib_runs = glob.glob(cf.calib_path+"noise*.dat")
    crop = len(cf.calib_path)+6 #for 'noise_'
    calib_runs_nb = [elem[crop:] for elem in calib_runs]
    calib_runs_nb = [int(elem[:elem.find("_")]) for elem in calib_runs_nb]

    if(run >= calib_runs_nb[-1]):
        return calib_runs[-1]
    else:
        for irun in range(len(calib_runs_nb)-1):
            if(calib_runs_nb[irun+1] > run):
                return calib_runs[irun]
    return calib_runs[0] #in case of problems

def map_reference_pedestal(run):
    reference_file = get_closest_ped_run(run)
    print("Will use calibration run : ", reference_file)

    with open(reference_file, "r") as pedcalib:
        for iline in pedcalib:
            li = iline.split()
            daqch = int(li[0])
            ped = float(li[7])
            rms = float(li[9])
            dc.map_ped[daqch].set_ref_pedestal(ped, rms)

def store_raw_ped_rms():
    
    compute_pedestal_mean()
    compute_pedestal_RMS()

    """ store the raw pedestal and rms """
    for i in range(cf.n_ChanTot):
        crp, view, ch = dc.map_ped[i].get_ana_chan()
        if(crp >= cf.n_CRPUsed): continue
        dc.map_ped[i].set_raw_pedestal(dc.ped_mean[crp,view,ch], dc.ped_rms[crp,view,ch])

def store_final_ped_rms():

    compute_pedestal_mean()
    compute_pedestal_RMS()

    for i in range(cf.n_ChanTot):
        crp, view, ch = dc.map_ped[i].get_ana_chan()
        if(crp >= cf.n_CRPUsed): continue
        dc.map_ped[i].set_evt_pedestal(dc.ped_mean[crp,view,ch], dc.ped_rms[crp,view,ch])

