import sys
import numpy as np
import numpy.ma as ma
import struct
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import time


import config as cf
import pedestals as ped
import channelmapper as cmap
import read_event as read
import plot_event as plot_ev
import noise_filter as noise



tstart = time.time()
data = open("1323_10_a.cosmics", "rb")

nevent = 1

topen = time.time()


cmap.ChannelMapper()
ped.MapRefPedestal()

#Run Header
fmt_h='II'
size_h = struct.calcsize(fmt_h)
run_nb, nb_evt = struct.unpack(fmt_h, data.read(size_h))
print "Run Nb: ", run_nb, " Nb of Events: ", nb_evt

size_i32 = struct.calcsize('I')
evtsz = nb_evt * 4 * size_i32
print evtsz


fmt_seq = 'IIII'
size_seq = struct.calcsize(fmt_seq)

sequence = []
for i in range(nb_evt):
    seq = struct.unpack(fmt_seq, data.read(size_seq))
    sequence.append(seq[1])



event_pos = []
event_pos.append( data.tell() )
for i in range(nb_evt-1):
    data.seek(sequence[i], 1)
    event_pos.append( data.tell() )

#End of run header reading part


trunhead = time.time()

size_header = 44
npalldata = np.zeros((2, 2, 960, 10000)) #crp, view, vchan



for ievent in range(nevent):
    print "-*-*-*-*-*-*-*-*-*-*-"
    print " READING EVENT ", ievent
    print "-*-*-*-*-*-*-*-*-*-*-"


    tevtread = time.time()
    idx = event_pos[ievent]

    v0, lro, cro = read.read_event_header(data, idx)
    v0.evt_nb_loc = 0

    tevthead = time.time()

    idx = idx + size_header
    npdatav0 = read.read_evt_uint32( data.read(cro) )

    idx = idx + cro + 1 #add the "bruno" byte

    v1, lro, cro = read.read_event_header(data, idx)
    v1.evt_nb_loc = 0
    npdatav1 = read.read_evt_uint32( data.read(cro) )

    tevtdata = time.time()

    print " time to open ", topen - tstart
    print " time to read run header ", trunhead - topen
    print " time to evt data ", tevtdata - tevtread
    print " TOTAL time ", tevtdata - tstart


    if( len(npdatav0)/cf.n_Sample != cf.n_ChanPerView):
        print " PBM OF Nb of CHANNELS in V0 !!! ", len(npdatav0)/cf.n_Sample , " vs ", cf.n_ChanPerView
    if( len(npdatav1)/cf.n_Sample != cf.n_ChanPerView):
        print " PBM OF Nb of CHANNELS in V1 !!! ", len(npdatav1)/cf.n_Sample , " vs ", cf.n_ChanPerView
        

    npdatav0 = np.split(npdatav0,cf.n_ChanPerView)
    npdatav1 = np.split(npdatav1,cf.n_ChanPerView)

    
    #reshape the array and subtract reference pedestal
    for idq in range(cf.n_ChanTot):
        crp, view, vch = cmap.DAQToCRP(idq)
        if(crp > 1): continue
        pedval = ped.GetPed(idq)
        if(crp < 0 or view < 0 or vch < 0):
            print " ERROR ? ", idq
        if(view==0):
            npalldata[crp,view,vch] = npdatav0[idq] - pedval            
        elif(view==1):
            npalldata[crp,view,vch] = npdatav1[idq-3840] - pedval
            
            
    print " done getting event ", ievent, " ! "
    traw = time.time()

    #npmask = ma.array(npalldata, mask=False, copy=False)

    mask = np.ones((2,2,960,10000), dtype=bool)

    #do we need to separate it ??? Can set npalldata at 0 directly
    alive_chan = np.ones((2,2,960,10000), dtype=bool)

    for ibrok in cf.broken_channels:
        crp, view, vch = cmap.DAQToCRP(ibrok)
        #npmask[crp, view, vch,:] = ma.masked
        alive_chan[crp, view, vch, : ] = False


    #broken_mask = ma.getmask(npmask).copy()
        
    signal_thresh = 3.
    adc_thresh    = 6.

    tfft = time.time()
    lowpasscut = 0.1 #MHz    
    freqlines  = [0.0234]

    npalldata = noise.FFTLowPass(npalldata, lowpasscut, freqlines)
    
    #npmask = noise.FFTLowPass(npmask, lowpasscut, freqlines)    
    print " time to fft %.2f"%( time.time() - tfft)

    """
    ##Simple ROI finder, only ADC threshold-based
    # commented method as it was super long
    npmask = ma.masked_where( (npmask >= adc_thresh) | broken_mask, npmask)
    ##Improve ROI finder based on pedrms
    ped_rms = noise.get_RMS(npmask)
    for ik in range(2): #two iterations seems to be enough
        npmask = ma.masked_where( (npmask >= signal_thresh*ped_rms[:,:,:,None]) | broken_mask, npmask)
        ped_rms = noise.get_RMS(npmask)    
    npmask = noise.coherent_filter_prev(npmask, 64)
    """
    

    t3 = time.time()
    #Apply coherent filter - ROI is done in the function
    npalldata, mask = noise.coherent_filter(npalldata, alive_chan, 64, signal_thresh)
    

    
    print " time to coh filt %.2f"%( time.time() - t3)


    #recompute the mask based on threshold above the pedestal    
    """
    #npmask.mask = tmp_mask
    npmask = ma.masked_where( (npmask >= adc_thresh) | broken_mask, npmask)
    ped_rms = noise.get_RMS(npmask)    
    for ik in range(2): #two iterations seems to be enough
        npmask = ma.masked_where( (npmask >= signal_thresh*ped_rms[:,:,:,None]) | broken_mask, npmask)                
        ped_rms = noise.get_RMS(npmask)        
    #print " time to converge on pedestal/masking ", time.time() - tped
    npmask.mask = (~npmask.mask | broken_mask)
    #at this point we have only hit-like points
    """


    t4 = time.time()
    #mask = np.where( npalldata > adc_thresh, True, False)
    ped_rms = noise.get_RMS(npalldata*mask)
    for it in range(2):
        #mask = np.where( (npalldata > signal_thresh*ped_rms[:,:,:,None]) | (mask), False, True) ## do not work somehow
        mask = np.where( (npalldata > signal_thresh*ped_rms[:,:,:,None]) , False, True)
        mask = mask | alive_chan
        ped_rms = noise.get_RMS(npalldata*mask)

    print " time to ROI %.2f"%(time.time() - t4)
    
    ROI = np.array(~mask, dtype=int)
    ROI *= 50
    #plot_ev.plot_event_display(npalldata, v0.run_nb, v0.evt_nb_glob)
    plot_ev.plot_event_display(ROI, v0.run_nb, v0.evt_nb_glob)
    #plot_ev.plot_pedestal([ped_rms_raw, ped_rms_fft, ped_rms_coh], v0.run_nb, v0.evt_nb_glob)
    plot_ev.plot_pedestal([ped_rms], ['Final'],['r'], v0.run_nb, v0.evt_nb_glob)

tottime = time.time() - tstart

print " TOTAL RUNNING TIME %.2f s == %.2f evt/s"% (tottime, tottime/nevent)
