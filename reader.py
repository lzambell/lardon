import sys
import os
import time 

import numpy as np
import tables as tables

import config as cf
import data_containers as dc

import pedestals as ped
import channelmapper as cmap
import read_event as read
import noise_filter as noise
import hitfinder as hf
import clustering as clus
import track_2d as trk2d
import track_3d as trk3d
import store as store
import read_mc as mc

import plotting as plot

import reconstruction_parameters as reco

def need_help():
    print("Usage: python reader.py ")
    print(" -run <run number ex:1323> ")
    print(" -sub <sub file ex: 10_a> ")
    print(" -n   <number of event to process>  [default (or -1) is all]")
    print(" -out <output name optn>")
    print(" -reco <reconstruction parameter's file> [default is default_reco.yaml]")
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
            

outname_option = ""
nevent = -1 
evt_type = "cosmics"          
mc_file  = ""
addMC = False
reco_param = cf.default_reco

for index, arg in enumerate(sys.argv):
    if arg in ['-run'] and len(sys.argv) > index + 1:
        run_n = sys.argv[index + 1]
    elif arg in ['-sub'] and len(sys.argv) > index + 1:
        evt_file = sys.argv[index + 1]
    elif arg in ['-n'] and len(sys.argv) > index + 1:
        nevent = int(sys.argv[index + 1])
    elif arg in ['-reco'] and len(sys.argv) > index + 1:
        reco_param = sys.argv[index + 1]
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
    print(" ERROR ! data file ", name_in, " do not exists ! ")
    sys.exit()



if(os.path.exists(reco_param) is False):
    print(" ERROR ! reco file ", reco_param, " do not exists ! ")
    sys.exit()
else:
    myreco = reco.Reco(reco_param)

data = open(name_in, "rb")
output = tables.open_file(name_out, mode="w", title="Reconstruction Output")


#output_test = tables.open_file(name_out.replace(".h5", "_lite.h5"), mode="w", title="Reconstruction Output Light")
#store.create_temp(output_test)



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

    n_bad_ch = ped.store_raw_ped_rms(myreco.param['bad_evt_thresh'])
    print("Nb of bad channels : ", n_bad_ch)

    if(n_bad_ch > myreco.param['bad_evt_nchan']):
        print("! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ")
        print(" EVENT LOOKS BAD .... SKIPPING")
        print("! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ")
        print("nb of bad channels ", n_bad_ch)
        #plot.plot_ed_data(option="microspark", to_be_shown=False)
        gr = store.new_event(output, ievent)
        store.store_event(output, gr)
        store.store_pedestal(output, gr)
        store.store_hits(output, gr)
        store.store_tracks2D(output, gr)
        store.store_tracks3D(output, gr)
        continue

    
    
    print("time to compute pedestals : %.3f s"%(time.time() - t_ped_raw))
    #plot.plot_ed_data(option="", to_be_shown=True)
    #plot.plot_ed_one_crp(0, "test", True)
    #plot.plot_ed_one_crp(1, "test", True)
    #plot.plot_ed_one_crp(3, "test", True)
    #plot.plot_wvf_single_current([(1, 1, 424), (1,1,425),(1,1,426),(3,0,259), (3,0,260),(3,0,261)], option="test", to_be_shown=True)
    #continue
    tfft = time.time()

    #ps = 
    noise.FFTLowPass(myreco.param['fft_lowpasscut'], myreco.param['fft_freqlines'])


    #noise.FFT2D()

    print(" time to fft %.2f"%( time.time() - tfft))
    tadc = time.time()
    #plot.plot_wvf_multi_current([(3,0,i) for i in range(60,91,5)], to_be_shown=True, option="test")


    """ 1st ROI attempt based on ADC cut + broken channels """
    noise.define_ROI_ADC(myreco.param['roi_adc_thresh'])
    print("adc cut roi based %.2f"%(time.time() - tadc))


    """ Update ROI based on ped rms """
    troi = time.time()
    noise.define_ROI(myreco.param['roi_signal_thresh'], myreco.param['roi_n_iter']) 
    t3 = time.time()
    print(" 1st ROi : %.2f"% (t3-troi))


    """Apply coherent filter(s) """
    noise.coherent_filter(myreco.param['coherent_groups'])

    print(" time to coh filt %.2f"%( time.time() - t3))
    #wvf_filt = dc.data[1,0,380,:].copy()    
    #plot.plot_wvf_evo([wvf_raw, wvf_fft, wvf_filt], title="CRP 1 - View 0, Ch 380", legends=['raw', 'fft', 'filt'], to_be_shown=True)
    

    t4 = time.time()    

    """ Update ROI regions """
    noise.define_ROI(myreco.param['roi_signal_thresh'], myreco.param['roi_n_iter'])
    print(" time to ROI %.2f"%(time.time() - t4))

    t5 = time.time()
    noise.median_filter(myreco.param['median_window'])

    print("time to median filter ", time.time()-t5)
    t6 = time.time()

    #plot.plot_ed_data(option="pres_filt", to_be_shown=False)

    """ Update ROI regions """
    noise.define_ROI(myreco.param['roi_signal_thresh_2'], myreco.param['roi_n_iter'])
    print(" time to ROI %.2f"%(time.time() - t6))


    """ final pedestal mean and rms """
    ped.store_final_ped_rms()

    t6 = time.time()
    """ parameters : pad left (n ticks) pad right (n ticks), min dt, thr1, thr2 """
    hf.hit_finder(myreco.param['hit_pad_left'], 
                  myreco.param['hit_pad_right'], 
                  myreco.param['hit_min_dt'], 
                  myreco.param['hit_thr1'], 
                  myreco.param['hit_thr2'])

    print(" time to Hit Search %.3f"%(time.time() - t6))
    
    t7 = time.time()

    
    """ 1st search for most of the tracks"""
    if(myreco.param['clus_use_dbscan']):
        """parameters : eps (cm), min pts, y axis squeeze"""
        for e, n, s in zip(myreco.param['clus_dbscan_eps'], 
                           myreco.param['clus_dbscan_npts'], 
                           myreco.param['clus_dbscan_ysqueez']):
            clus.dbscan(e, n, s)

    if(myreco.param['clus_use_mst']):
        """ parameters : search radius (cm), min nb of hits in cluster """
        clus.mst(myreco.param['clus_mst_radius'],
                 myreco.param['clus_mst_npts'])

    tclus = time.time()
    print("time to cluster %.3f"%(tclus-t7)) 
    #plot.plot_2dcrp_clusters(option="mst", to_be_shown=True)


    t8 = time.time()


    for n, r, c in zip(myreco.param['trk2d_npts'],
                       myreco.param['trk2d_rcut'],
                       myreco.param['trk2d_chi2']):

        """parameters : min nb hits, rcut, chi2cut, y error, slope error, pbeta"""
        trk2d.find_tracks(n, r, c, 
                          myreco.param['pfilt_posErr'], 
                          myreco.param['pfilt_slopeErr'], 
                          myreco.param['pfilt_pbeta'])


    t9 = time.time()
    print("time to find tracks %.3f"%(t9-t8))
    

    tst = time.time()
    """ parameters are : min distance in between 2 tracks end points, slope error tolerance, extrapolated distance tolerance + filter input parameters"""
    trk2d.stitch_tracks(myreco.param['trk2d_stitch_dmin'],
                        myreco.param['trk2d_stitch_slope'],
                        myreco.param['trk2d_stitch_dma'],
                        myreco.param['pfilt_posErr'], 
                        myreco.param['pfilt_slopeErr'], 
                        myreco.param['pfilt_pbeta'])

    print("time to stitch tracks %.3f"%(time.time()-tst))
    #plot.plot_2dview_2dtracks(option='filter', to_be_shown=False)
    
    t3d = time.time()
    """ parameters are : z start/end agreement cut (cm), v0/v1 charge balance, distance to detector boundaries for time correction (cm) """
    trk3d.find_tracks(myreco.param['trk3d_ztol'],
                      myreco.param['trk3d_qbal'],
                      myreco.param['trk3d_dbound'])
    
    print("Time build 3D tracks %.3f"%(time.time() - t3d))

    #plot.plot_2dview_hits_and_3dtracks(option="pres", to_be_shown=False)
    #plot.plot_3d("pres", False)

    dc.evt_list[-1].dump_reco()

    tstore = time.time()
    gr = store.new_event(output, ievent)

    store.store_event(output, gr)
    store.store_pedestal(output, gr)
    store.store_hits(output, gr)
    store.store_tracks2D(output, gr)
    store.store_tracks3D(output, gr)
    #store.store_tracks3D_test(output_test)
    #print("time to store %.3f"%(time.time()-tstore))

data.close()
output.close()
#output_test.close()
tottime = time.time() - tstart
print(" TOTAL RUNNING TIME %.2f s == %.2f s/evt"% (tottime, tottime/nevent))
