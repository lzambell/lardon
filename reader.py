import sys
import os
import numpy as np
import time 
#from tables import *
import tables as tables

import config as cf
import data_containers as dc

import pedestals as ped
import channelmapper as cmap
import read_event as read
import plot_event as plot_ev
import noise_filter as noise
import hitfinder as hf
import clustering as clus
import track_2d as trk
import store as store


def need_help():
    print("Usage: python reader.py ")
    print(" -run <run number ex:1323> ")
    print(" -sub <sub file ex: 10_a> ")
    print(" -n   <number of event to process>  [default (or -1) is all]")
    print(" -type <evt type cosmics/ped/...> [default is cosmics]")
    print(" -h print this message")
    sys.exit()
    

if len(sys.argv) == 1:
    need_help()
else:
    for index, arg in enumerate(sys.argv):
        if arg in ['-h'] :
            need_help()
            

""" Reconstruction parameters """
lowpasscut     = 0.1 #MHz    
freqlines      = [0.00125, 0.0234] #in MHz
signal_thresh  = 4.
signal_thresh_2  = 2.5
adc_thresh     = 6.
coherent_groups = [64,320]


nevent = -1 
evt_type = "cosmics"          

for index, arg in enumerate(sys.argv):
    if arg in ['-run'] and len(sys.argv) > index + 1:
        run_n = sys.argv[index + 1]
    elif arg in ['-sub'] and len(sys.argv) > index + 1:
        evt_file = sys.argv[index + 1]
    elif arg in ['-n'] and len(sys.argv) > index + 1:
        nevent = int(sys.argv[index + 1])
    elif arg in ['-type'] and len(sys.argv) > index + 1:
        evt_type = sys.argv[index + 1]
                

tstart = time.time()

name_in = cf.data_path + run_n + "/" + run_n + "_" + evt_file + "." + evt_type
name_out = cf.store_path + "/" + run_n + "_" + evt_file + ".h5"

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
ped.MapRefPedestal(run_nb)


if(nevent > nb_evt or nevent < 0):
    nevent = nb_evt

print(" --->> Will process ", nevent, " events [ out of ", nb_evt, "] of run ", run_nb)

sequence = []
for i in range(nb_evt):
    seq  = np.fromfile( data, dtype='<u4', count=4)
    # 4 uint of [event number - event total size with header- event data size - 0]
    sequence.append(seq[1])

event_pos = []
event_pos.append( data.tell() )
for i in range(nb_evt-1):
    data.seek(sequence[i], 1)
    event_pos.append( data.tell() ) #get the byte position of each event
""" End of run header reading part """


for ibrok in cf.daq_broken_channels:
    crp, view, vch = cmap.DAQToCRP(ibrok)
    dc.alive_chan[crp, view, vch, : ] = False


for crp, view, vch in cf.crp_broken_channels:
    dc.alive_chan[crp, view, vch, : ] = False


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

    #plot_ev.plot_event_display("raw")

    if(run_nb <= cf.run_inv_signal):
        dc.data *= -1.    
    

    tfft = time.time()
    noise.FFTLowPass(lowpasscut, freqlines)    
    print(" time to fft %.2f"%( time.time() - tfft))

    tadc = time.time()
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


    t4 = time.time()

    """ Update ROI regions """
    noise.define_ROI(signal_thresh, 2)

    print(" time to ROI %.2f"%(time.time() - t4))


    """invert the mask, so now True is signal (and not broken channels)"""    
    #ROI = np.array(~mask & alive_chan, dtype=bool)

    
    t6 = time.time()
    hf.HitFinder(noise.get_RMS())


    print(" time to Hit Search %.3f"%(time.time() - t6))

    
    t7 = time.time()

    
    """ 1st search for most of the tracks"""
    """parameters : eps, min pts, y axis squeeze"""
    clus.dbscan(20, 15, 0.1)

    """2nd search for vertical tracks not yet clustered """
    clus.dbscan(30, 5, 0.1)

    print("time to cluster %.3f"%(time.time()-t7)) 

    #plot_ev.plot_hits_clustered()
    #plot_ev.plot_hits_var()
    #plot_ev.plot_event_display("filt")
    #plot_ev.plot_hits_view()

    t8 = time.time()
    
    """parameters : rcut, chi2cut, y error, slope error, pbeta"""
    trk.FindTracks(6., 12., 0.3125, 1., 3.)


    print("time to find tracks %.3f"%(time.time()-t8))
    #plot_ev.plot_tracks2D()
    #plot_ev.plot_track2D_var()

    store.store_tracks(output, ievent)


    """
    ped_fin = noise.get_RMS(npalldata*mask)
    plot_ev.plot_pedestal([ped_ref, ped_ini, ped_fin], ['Ref.', 'Raw', 'Final'],['r','k','b'])
    """

data.close()
output.close()
tottime = time.time() - tstart
print(" TOTAL RUNNING TIME %.2f s == %.2f evt/s"% (tottime, tottime/nevent))
