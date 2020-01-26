import sys
import numpy as np
#import numpy.ma as ma
#import struct
#import matplotlib.pyplot as plt
#from matplotlib.colors import LinearSegmentedColormap
import time
#from datetime import datetime

import config as cf
import pedestals as ped
import channelmapper as cmap
import read_event as read
import plot_event as plot_ev
import noise_filter as noise


tstart = time.time()
data = open("1323_10_a.cosmics", "rb")

nevent = 1
lowpasscut = 0.1 #MHz    
freqlines  = [0.0234]
signal_thresh = 3.5
adc_thresh    = 6.



topen = time.time()


cmap.ChannelMapper()
ped.MapRefPedestal()

#Run Header

run_nb, nb_evt = np.fromfile(data, dtype='<u4', count=2)
"""

data.seek(0)
fmt_h='II'
size_h = struct.calcsize(fmt_h)
run_nb, nb_evt = struct.unpack(fmt_h, data.read(size_h))
print "Run Nb: ", run_nb, " Nb of Events: ", nb_evt
"""
"""
size_i32 = struct.calcsize('I')
evtsz = nb_evt * 4 * size_i32
print evtsz
"""
"""
fmt_seq = 'IIII'
size_seq = struct.calcsize(fmt_seq)
"""

sequence = []
for i in range(nb_evt):
    #seq = struct.unpack(fmt_seq, data.read(size_seq))
    seq  = np.fromfile( data, dtype='<u4', count=4)
    # 4 uint of [event number - event total size with header- event data size - 0]
    sequence.append(seq[1])

event_pos = []
event_pos.append( data.tell() )
for i in range(nb_evt-1):
    data.seek(sequence[i], 1)
    event_pos.append( data.tell() ) #get the byte position of each event

#End of run header reading part


trunhead = time.time()

size_header = 44
npalldata = np.zeros((2, cf.n_View, cf.n_ChanPerCRP, cf.n_Sample)) #crp, view, vchan



for ievent in range(nevent):
    print "-*-*-*-*-*-*-*-*-*-*-"
    print " READING EVENT ", ievent
    print "-*-*-*-*-*-*-*-*-*-*-"


    tevtread = time.time()
    idx = event_pos[ievent]

    #v0, lro, cro = read.read_event_header_np(data, idx)
    #v0, lro, cro = read.read_event_header(data, idx)
    #v0.evt_nb_loc = 0



    #tevthead = time.time()

    #idx = idx + size_header
    #npdatav0 = read.read_evt_uint32( data.read(cro) )

    #idx = idx + cro + 1 #add the "bruno" byte

    #v1, lro, cro = read.read_event_header_np(data, idx)
    #v1.evt_nb_loc = 0
    #npdatav1 = read.read_evt_uint32( data.read(cro) )
    #data.seek(idx,0)
    evt, npdatav0, npdatav1 = read.read_event( data, idx)
    evt.evt_nb_loc = ievent
    
    print "RUN ",evt.run_nb, " EVENT ", evt.evt_nb_loc, " / ", evt.evt_nb_glob," is ", "good " if evt.evt_flag==True else "bad "
    print "Taken at ", time.ctime(evt.time_s), " + ", evt.time_ns, " ns "


    tevtdata = time.time()
    ### this part could  be simplified / merged

    #print " time to open ", topen - tstart
    #print " time to read run header ", trunhead - topen
    #print " time to evt data ", tevtdata - tevtread
    print " TOTAL time ", tevtdata - tevtread


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
    totlength = 2*2*960*10000

    #the mask will be used to differentiate background (True for noise processing) from signal (False for noise processing)
    #at first everything is considered background (all at True)
    mask = np.ones((2, cf.n_View, cf.n_ChanPerCRP, cf.n_Sample), dtype=bool)
    #print "initial mask has ", totlength-np.count_nonzero(mask), " hit-like points"

    #alive_chan mask intends to not take into account broken channels
    #True : not broken
    #False : broken
    #do we need to separate it ??? Can set npalldata at 0 directly -->NO would be included in noise removal then
    alive_chan = np.ones((2,2,960,10000), dtype=bool)

    for ibrok in cf.broken_channels:
        crp, view, vch = cmap.DAQToCRP(ibrok)
        #npmask[crp, view, vch,:] = ma.masked
        alive_chan[crp, view, vch, : ] = False
    
    #print " dead channel mask has ", totlength - np.count_nonzero(alive_chan), " not ok points, ", np.count_nonzero(alive_chan), " ok points"

    #broken_mask = ma.getmask(npmask).copy()
        

    tfft = time.time()

    npalldata = noise.FFTLowPass(npalldata, lowpasscut, freqlines)
    
    print " time to fft %.2f"%( time.time() - tfft)
    

    t3 = time.time()
    #Apply coherent filter - ROI is done in the function
    npalldata, mask = noise.coherent_filter(npalldata, alive_chan, 64, adc_thresh, signal_thresh)
    

    
    print " time to coh filt %.2f"%( time.time() - t3)


    #recompute the mask based on threshold above the pedestal    


    t4 = time.time()
    ped_rms = noise.get_RMS(npalldata*mask)
    for it in range(2):
        mask = np.where( (npalldata > signal_thresh*ped_rms[:,:,:,None]) | (~alive_chan) , False, True)
        #print " iteration ", it, " nb of roi ", totlength-np.count_nonzero(mask)
        #mask = mask | alive_chan
        #print " adding broken channels nb of roi ", totlength-np.count_nonzero(mask)
        ped_rms = noise.get_RMS(npalldata*mask)

    print " time to ROI %.2f"%(time.time() - t4)
    
    ROI = np.array(~mask, dtype=int)
    ROI = ROI & alive_chan
    ROI *= 50
    plot_ev.plot_event_display(npalldata, evt.run_nb, evt.evt_nb_glob,"filt")
    plot_ev.plot_event_display(ROI, evt.run_nb, evt.evt_nb_glob, "ROI")
    plot_ev.plot_pedestal([ped_rms], ['Final'],['r'], evt.run_nb, evt.evt_nb_glob)



tottime = time.time() - tstart
print " TOTAL RUNNING TIME %.2f s == %.2f evt/s"% (tottime, tottime/nevent)
