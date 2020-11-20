import sys
import os
import numpy as np
import time 
import tables as tables

import config as cf
import data_containers as dc

import pedestals as ped
import channelmapper as cmap
import read_event as read
#import plot_event as plot_ev
import noise_filter as noise
import hitfinder as hf
import clustering as clus
import track_2d as trk2d
import track_3d as trk3d
import store as store
import read_mc as mc

import plotting as plot


def need_help():
    print("Usage: python reader.py ")
    print(" -run <run number ex:1323> ")
    print(" -sub <sub file ex: 10_a> ")
    print(" -n   <number of event to process>  [default (or -1) is all]")
    print(" -out <output name optn>")
    print(" -type <evt type cosmics/ped/...> [default is cosmics]")
    print(" -mc <simulation root file to be added to noise>") 
    print(" -h print this message")

    sys.exit()
    

if len(sys.argv) == 1:
    need_help()
else:
    for index, arg in enumerate(sys.argv):
        if arg in ['-h'] :
            need_help()
            

""" Reconstruction parameters """
lowpasscut       = 0.06 #0.1 #MHz    
freqlines        = []#0.0234, 0.0625, 0.0700] #in MHz
signal_thresh    = 4.
signal_thresh_2  = 2.5
adc_thresh       = 6.
coherent_groups  = [320, 64]#, 8]#[64,32,8]

outname_option = ""
nevent = -1 
evt_type = "cosmics"          
mc_file  = ""
addMC = False

for index, arg in enumerate(sys.argv):
    if arg in ['-run'] and len(sys.argv) > index + 1:
        run_n = sys.argv[index + 1]
    elif arg in ['-sub'] and len(sys.argv) > index + 1:
        evt_file = sys.argv[index + 1]
    elif arg in ['-n'] and len(sys.argv) > index + 1:
        nevent = int(sys.argv[index + 1])
    elif arg in ['-type'] and len(sys.argv) > index + 1:
        evt_type = sys.argv[index + 1]
    elif arg in ['-out'] and len(sys.argv) > index + 1:
        outname_option = sys.argv[index + 1]
    elif arg in ['-mc'] and len(sys.argv) > index + 1:
        addMC = True
        mc_file = sys.argv[index + 1]
    

tstart = time.time()

name_in = cf.data_path + run_n + "/" + run_n + "_" + evt_file + "." + evt_type

if(outname_option):
    outname_option = "_"+outname_option
else:
    outname_option = ""
    
#name_out = cf.store_path +"/" + run_n + "/" + run_n + "_" + evt_file + outname_option + ".h5"
name_out = cf.store_path + "/" + run_n + "_" + evt_file + outname_option + ".h5"

if(os.path.exists(name_in) is False):
    print(" ERROR ! file ", name_in, " do not exists ! ")
    sys.exit()

data = open(name_in, "rb")
output = tables.open_file(name_out, mode="w", title="Reconstruction Output")



""" Reading Run Header """
run_nb, nb_evt = np.fromfile(data, dtype='<u4', count=2)

""" Build DAQ Channel <-> Analysis Channel correspondance """
cmap.ChannelMapper()

""" Get Reference Pedestals """
ped.map_reference_pedestal(run_nb)


if(nevent > nb_evt or nevent < 0):
    nevent = nb_evt

print(" --->> Will process ", nevent, " events [ out of ", nb_evt, "] of run ", run_nb)

if(addMC):
    print(" Will add MC from file :", mc_file)
    the_mc = mc.readmc(mc_file)
    
    
store.store_infos(output, run_n, evt_file, nevent, time.time())
 
sequence = []
for i in range(nb_evt):
    seq  = np.fromfile( data, dtype='<u4', count=4)
    """4 uint of [event number - event total size with header- event data size - 0]"""
    sequence.append(seq[1])

event_pos = []
event_pos.append( data.tell() )
for i in range(nb_evt-1):
    data.seek(sequence[i], 1)
    """ get the byte position of each event """
    event_pos.append( data.tell() ) 
""" End of run header reading part """


for ibrok in cf.daq_broken_channels:
    crp, view, vch = cmap.DAQToCRP(ibrok)
    dc.alive_chan[crp, view, vch, : ] = False


for crp, view, vch in cf.crp_broken_channels:
    if(crp >= cf.n_CRPUsed): continue
    dc.alive_chan[crp, view, vch, : ] = False

plot.set_style()
for ievent in range(nevent):
    print("-*-*-*-*-*-*-*-*-*-*-")
    print(" READING EVENT ", ievent)
    print("-*-*-*-*-*-*-*-*-*-*-")

    dc.reset_event()

    tevtread = time.time()
    idx = event_pos[ievent]

    if( read.read_event( data, idx, ievent) < 0):
        print(" there is a pbm with the data")
        continue

    dc.evt_list[-1].dump()

    tevtdata = time.time()

    print(" -> Reading time %.2f s"%( tevtdata - tevtread))


    if(run_nb <= cf.run_inv_signal):
        dc.data *= -1.    
    
    elif(cf.n_CRPUsed > 2):
        #should ask dario if this is fixed and if so, when
        dc.data[3,:,:,:] *= -1.

    if(addMC):
        the_mc.read_event(ievent)

    t_ped_raw = time.time()

    ped.store_raw_ped_rms()
    
    print("time to compute pedestals : %.3f s"%(time.time() - t_ped_raw))
    #plot.plot_ed_data("test", False)
    #plot.plot_ed_one_crp(0, "test", False)
    #plot.plot_wvf_single_current([(1, 1, 424), (1,1,425),(1,1,426),(3,0,259), (3,0,260),(3,0,261)], option="test", to_be_shown=True)
    tfft = time.time()

    #ps = 
    noise.FFTLowPass(lowpasscut, freqlines)


    #noise.FFT2D()

    print(" time to fft %.2f"%( time.time() - tfft))
    tadc = time.time()
    #plot.plot_wvf_multi_current([(3,0,i) for i in range(60,91,5)], to_be_shown=True, option="test")


    """ 1st ROI attempt based on ADC cut + broken channels """
    noise.define_ROI_ADC(adc_thresh)
    print("adc cut roi based %.2f"%(time.time() - tadc))


    """ Update ROI based on ped rms """
    troi = time.time()
    noise.define_ROI(signal_thresh, 2)    
    t3 = time.time()
    print(" 1st ROi : %.2f"% (t3-troi))


    """Apply coherent filter(s) """
    noise.coherent_filter(coherent_groups)

    print(" time to coh filt %.2f"%( time.time() - t3))
    #wvf_filt = dc.data[1,0,380,:].copy()    
    #plot.plot_wvf_evo([wvf_raw, wvf_fft, wvf_filt], title="CRP 1 - View 0, Ch 380", legends=['raw', 'fft', 'filt'], to_be_shown=True)
    

    t4 = time.time()    

    """ Update ROI regions """
    noise.define_ROI(signal_thresh, 2)
    print(" time to ROI %.2f"%(time.time() - t4))

    t5 = time.time()
    noise.median_filter(400)

    print("time to median filter ", time.time()-t5)
    t6 = time.time()

    """ Update ROI regions """
    noise.define_ROI(signal_thresh_2, 2)
    print(" time to ROI %.2f"%(time.time() - t6))


    """ final pedestal mean and rms """
    ped.store_final_ped_rms()

    
    t6 = time.time()
    """ parameters : pad left (n ticks) pad right (n ticks), min dt, thr1, thr2 """
    hf.hit_finder(5, 10, 20, 3., 6.)

    print(" time to Hit Search %.3f"%(time.time() - t6))
    
    t7 = time.time()

    
    """ 1st search for most of the tracks"""
    """parameters : eps (cm), min pts, y axis squeeze"""
    #clus.dbscan(15, 10, 1.)

    """2nd search for vertical tracks not yet clustered """
    #clus.dbscan(30, 5, 0.05)

    
    """ parameters : search radius (cm), min nb of hits in cluster """
    clus.mst(10, 5)

    tclus = time.time()
    print("time to cluster %.3f"%(tclus-t7)) 
    #plot.plot_2dcrp_clusters(option="mst", to_be_shown=True)


    t8 = time.time()

    """parameters : min nb hits, rcut, chi2cut, y error, slope error, pbeta"""
    trk2d.find_tracks(10, 6., 8., 0.3125, 1., 3., 1)

    """parameters : min nb hits, rcut, chi2cut, y error, slope error, pbeta"""
    trk2d.find_tracks(4, 8., 12., 0.3125, 1., 3., 2)


    t9 = time.time()
    print("time to find tracks %.3f"%(t9-t8))
    

    tst = time.time()
    """ parameters are : min distance in between 2 tracks end points, slope error tolerance, extrapolated distance tolerance + filter input parameters"""
    trk2d.stitch_tracks(50., 10., 6., 0.3125, 1., 3., 1)

    print("time to stitch tracks %.3f"%(time.time()-tst))
    #plot.plot_2dview_2dtracks(option='test', to_be_shown=False)
    
    t3d = time.time()
    """ parameters are : z start/end agreement cut (cm), v0/v1 charge balance """

    trk3d.find_tracks(8., 0.25)
    print("Time build 3D tracks %.3f"%(time.time() - t3d))

    #plot.plot_2dview_hits_and_3dtracks(option="filter", to_be_shown=False)
    #plot.plot_3d("filter", False)

    dc.evt_list[-1].dump_reco()

    tstore = time.time()
    gr = store.new_event(output, ievent)

    store.store_event(output, gr)
    store.store_pedestal(output, gr)
    store.store_hits(output, gr)
    store.store_tracks2D(output, gr)
    store.store_tracks3D(output, gr)
    print("time to store %.3f"%(time.time()-tstore))

data.close()
output.close()
tottime = time.time() - tstart
print(" TOTAL RUNNING TIME %.2f s == %.2f s/evt"% (tottime, tottime/nevent))
