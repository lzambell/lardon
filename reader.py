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
import plot_event as plot_ev
import noise_filter as noise
import hitfinder as hf
import clustering as clus
import track_2d as trk2d
import track_3d as trk3d
import store as store


def need_help():
    print("Usage: python reader.py ")
    print(" -run <run number ex:1323> ")
    print(" -sub <sub file ex: 10_a> ")
    print(" -n   <number of event to process>  [default (or -1) is all]")
    print(" -out <output name optn>")
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
lowpasscut       = 0.06 #0.1 #MHz    
freqlines        = []#0.0234, 0.0625, 0.0700] #in MHz
signal_thresh    = 4.
signal_thresh_2  = 2.5
adc_thresh       = 6.
coherent_groups  = [64, 32, 8]
outname_option = ""

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
    elif arg in ['-out'] and len(sys.argv) > index + 1:
        outname_option = sys.argv[index + 1]
                

tstart = time.time()

name_in = cf.data_path + run_n + "/" + run_n + "_" + evt_file + "." + evt_type

if(outname_option):
    outname_option = "_"+outname_option
else:
    outname_option = ""
    
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
ped.MapRefPedestal(run_nb)


if(nevent > nb_evt or nevent < 0):
    nevent = nb_evt

print(" --->> Will process ", nevent, " events [ out of ", nb_evt, "] of run ", run_nb)
store.store_infos(output, run_n, evt_file, nevent, time.time())
 
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

    t_ped_raw = time.time()

    noise.compute_pedestal_mean()
    noise.compute_pedestal_RMS()

    for i in range(cf.n_ChanTot):
        crp, view, ch = dc.map_ped[i].get_ana_chan()
        if(crp >= cf.n_CRPUsed): continue
        dc.map_ped[i].set_raw_pedestal(dc.ped_mean[crp,view,ch], dc.ped_rms[crp,view,ch])
    

    print("time to compute pedestals : %.3f s"%(time.time() - t_ped_raw))

    """
    wvf_raw_crp0_v0 = [dc.data[0,0,100,:].copy(), dc.data[0,0,400,:].copy(), dc.data[0,0,850,:].copy()]
    wvf_raw_crp0_v1 = [dc.data[0,1,100,:].copy(), dc.data[0,1,400,:].copy(), dc.data[0,1,850,:].copy()]
    wvf_raw_crp1_v0 = [dc.data[1,0,100,:].copy(), dc.data[1,0,400,:].copy(), dc.data[1,0,850,:].copy()]
    wvf_raw_crp1_v1 = [ dc.data[1,1,100,:].copy(), dc.data[1,1,400,:].copy(), dc.data[1,1,850,:].copy()]
    """

    tfft = time.time()
    ps = noise.FFTLowPass(lowpasscut, freqlines)

    """
    if(ievent==0):
        ps_avg = ps/nevent #noise.FFTLowPass(lowpasscut, freqlines)/nevent
    else:
        ps_avg += ps/nevent #noise.FFTLowPass(lowpasscut, freqlines)/nevent
    """

    #plot_ev.plot_event_display("fft"+outname_option)

    """
    wvf_fft_crp0_v0 = [dc.data[0,0,100,:].copy(), dc.data[0,0,400,:].copy(), dc.data[0,0,850,:].copy()]
    wvf_fft_crp0_v1 = [dc.data[0,1,100,:].copy(), dc.data[0,1,400,:].copy(), dc.data[0,1,850,:].copy()]
    wvf_fft_crp1_v0 = [dc.data[1,0,100,:].copy(), dc.data[1,0,400,:].copy(), dc.data[1,0,850,:].copy()]
    wvf_fft_crp1_v1 = [ dc.data[1,1,100,:].copy(), dc.data[1,1,400,:].copy(), dc.data[1,1,850,:].copy()]
    """

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

    """
    wvf_coh_crp0_v0 = [dc.data[0,0,100,:].copy(), dc.data[0,0,400,:].copy(), dc.data[0,0,850,:].copy()]
    wvf_coh_crp0_v1 = [dc.data[0,1,100,:].copy(), dc.data[0,1,400,:].copy(), dc.data[0,1,850,:].copy()]
    wvf_coh_crp1_v0 = [dc.data[1,0,100,:].copy(), dc.data[1,0,400,:].copy(), dc.data[1,0,850,:].copy()]
    wvf_coh_crp1_v1 = [ dc.data[1,1,100,:].copy(), dc.data[1,1,400,:].copy(), dc.data[1,1,850,:].copy()]
    """
    
    """
    plot_ev.plot_waveform_evo([wvf_raw_crp0_v0, wvf_fft_crp0_v0, wvf_coh_crp0_v0], ['raw','fft','coh'], ['black','cyan','orange'],"crp0_v0"+outname_option)
    plot_ev.plot_waveform_evo([wvf_raw_crp0_v1, wvf_fft_crp0_v1, wvf_coh_crp0_v1], ['raw','fft','coh'], ['black','cyan','orange'],"crp0_v1"+outname_option)

    plot_ev.plot_waveform_evo([wvf_raw_crp1_v0, wvf_fft_crp1_v0, wvf_coh_crp1_v0], ['raw','fft','coh'], ['black','cyan','orange'],"crp1_v0"+outname_option)
    plot_ev.plot_waveform_evo([wvf_raw_crp1_v1, wvf_fft_crp1_v1, wvf_coh_crp1_v1], ['raw','fft','coh'], ['black','cyan','orange'],"crp1_v1"+outname_option)
    """

    print(" time to coh filt %.2f"%( time.time() - t3))

    t4 = time.time()

    """ Update ROI regions """
    noise.define_ROI(signal_thresh_2, 2)
    print(" time to ROI %.2f"%(time.time() - t4))

    """ final pedestal mean and rms """
    noise.compute_pedestal_mean()
    noise.compute_pedestal_RMS()
    for i in range(cf.n_ChanTot):
        crp, view, ch = dc.map_ped[i].get_ana_chan()
        if(crp >= cf.n_CRPUsed): continue
        dc.map_ped[i].set_evt_pedestal(dc.ped_mean[crp,view,ch], dc.ped_rms[crp,view,ch])




    """invert the mask, so now True is signal (and not broken channels)"""    
    #ROI = np.array(~mask & alive_chan, dtype=bool)

    
    t6 = time.time()

    """ parameters : pad left (n ticks) pad right (n ticks), min dt, thr1, thr2 """
    hf.hit_finder(5, 10, 10, 3., 4.)

    print(" time to Hit Search %.3f"%(time.time() - t6))
    #plot_ev.plot_waveform_hits(0,1,750)
    #plot_ev.plot_waveform_hits(1,0,100)
    
    t7 = time.time()

    
    """ 1st search for most of the tracks"""
    """parameters : eps, min pts, y axis squeeze"""
    clus.dbscan(20, 15, 0.05)

    """2nd search for vertical tracks not yet clustered """
    clus.dbscan(30, 5, 0.05)

    print("time to cluster %.3f"%(time.time()-t7)) 

    #plot_ev.plot_hits_clustered()
    #plot_ev.plot_hits_var()
    #plot_ev.plot_event_display("filt_fft0.06")
    #plot_ev.plot_hits_view()

    t8 = time.time()

    """parameters : min nb hits, rcut, chi2cut, y error, slope error, pbeta"""
    trk2d.find_tracks(10, 6., 8., 0.3125, 1., 3., 1)

    """parameters : min nb hits, rcut, chi2cut, y error, slope error, pbeta"""
    trk2d.find_tracks(4, 2., 12., 0.3125, 1., 3., 2)


    t9 = time.time()
    print("time to find tracks %.3f"%(t9-t8))

    #plot_ev.plot_tracks2D("raw")

    tst = time.time()
    """ parameters are : min distance in between 2 tracks end points, slope error tolerance, extrapolated distance tolerance"""
    trk2d.stitch_tracks(50., 3., 5.)
    print("time to stitch tracks %.3f"%(time.time()-tst))



    #plot_ev.plot_tracks2D("stitch")
    #plot_ev.plot_track2D_var()

    t3d = time.time()
    """ parameters are : z start/end agreement cut (cm), v0/v1 charge balance """
    trk3d.find_tracks(8., 0.25)
    print("Time build 3D tracks ", time.time() - t3d)
    #plot_ev.plot_tracks3D_proj()
    #plot_ev.plot_tracks3D()
    dc.evt_list[-1].dump_reco()    


    gr = store.new_event(output, ievent)

    store.store_event(output, gr)
    store.store_pedestal(output, gr)
    store.store_hits(output, gr)
    store.store_tracks2D(output, gr)
    store.store_tracks3D(output, gr)


    """
    ped_fin = noise.get_RMS(npalldata*mask)
    plot_ev.plot_pedestal([ped_ref, ped_ini, ped_fin], ['Ref.', 'Raw', 'Final'],['r','k','b'])
    """
#plot_ev.plot_event_fft(ps_avg, zmax=1.25, option="all")
#plot_ev.plot_event_fft(ps_avg, zmax=0.25, option="zoom")

"""
ps_avg = np.reshape(ps_avg, (2, 2, 15, 64, 5001))
ps_avg = np.einsum('ijklm->ijkm', ps_avg)/64.#ps_avg.sum(axis=3)


for icrp in range(2):
    for iview in range(2):
        np.savetxt("fft/"+run_n + "_" + evt_file +"_fft_f_crp"+str(icrp)+"_v"+str(iview)+".txt",ps_avg[icrp,iview,:,:], delimiter=',')

plot_ev.plot_event_fft(ps_avg, zmax=0.2, option="zoom_avg")
"""
data.close()
output.close()
tottime = time.time() - tstart
print(" TOTAL RUNNING TIME %.2f s == %.2f evt/s"% (tottime, tottime/nevent))
