import config as cf
import data_containers as dc
import glob 

"""
def dump(i):
    print(" DAQCH ", i)
    print(" CRP ", cf.reference[i].crp, " VIEW ", cf.reference[i].view, " CH ", cf.reference[i].vchan)
    print(" ped = ", cf.reference[i].ped, " +/- ", cf.reference[i].rms)
"""




def GetRefPed(i):
    return dc.map_ped[i].ref_ped

def GetRefPedRMS(i):
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

def MapRefPedestal(run):
    reference_file = get_closest_ped_run(run)
    print("Will use calibration run : ", reference_file)

    with open(reference_file, "r") as pedcalib:
        for iline in pedcalib:
            li = iline.split()
            daqch = int(li[0])
            ped = float(li[7])
            rms = float(li[9])
            dc.map_ped[daqch].set_ref_pedestal(ped, rms)

