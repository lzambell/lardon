import config as cf

def dump(i):
    print " DAQCH ", i
    print " CRP ", cf.reference[i].crp, " VIEW ", cf.reference[i].view, " CH ", cf.reference[i].vchan
    print " ped = ", cf.reference[i].ped, " +/- ", cf.reference[i].rms



def GetPed(i):
    return cf.map_ref[i].ped
def GetPedRMS(i):
    return cf.map_ref[i].rms


def MapRefPedestal():
    with open("/afs/cern.ch/user/n/np02onlp/public/calib/pedestals/noise_1317_5_b.dat", "r") as pedcalib:
        for iline in pedcalib:
            li = iline.split()
            daqch = int(li[0])
            cf.map_ref[daqch].ped = float(li[7])
            cf.map_ref[daqch].rms = float(li[9])
            

