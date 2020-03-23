import sys
import os
import numpy as np
import time 


import matplotlib.pyplot as plt
import config as cf
import pedestals as ped
import channelmapper as cmap
import read_event as read
import plot_event as plot_ev
import noise_filter as noise
import deconv as dec

"""this part needs sklearn"""
#import clustering as clus
#import hitfinder as hf


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
            

""" Analysis parameters """
lowpasscut     = 0.1 #MHz    
freqlines      = [0.00125, 0.0234] #in MHz
signal_thresh  = 4.
adc_thresh     = 6.
coherent_group = 64


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
name = cf.data_path + run_n + "/" + run_n + "_" + evt_file + "." + evt_type

if(os.path.exists(name) is False):
    print(" ERROR ! file ", name, " do not exists ! ")
    sys.exit()

data = open(name, "rb")

cmap.ChannelMapper()




""" Reading Run Header """
run_nb, nb_evt = np.fromfile(data, dtype='<u4', count=2)
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



npalldata = np.zeros((2, cf.n_View, cf.n_ChanPerCRP, cf.n_Sample)) #crp, view, vchan


"""
the mask will be used to differentiate background (True for noise processing) from signal (False for noise processing)
at first everything is considered background (all at True)
"""
mask = np.ones((2, cf.n_View, cf.n_ChanPerCRP, cf.n_Sample), dtype=bool)

"""
alive_chan mask intends to not take into account broken channels
True : not broken
False : broken
TO DO : make a run-dependent definition of broken channels
"""

alive_chan = np.ones((2,cf.n_View, cf.n_ChanPerCRP, cf.n_Sample), dtype=bool)
for ibrok in cf.broken_channels:
    crp, view, vch = cmap.DAQToCRP(ibrok)
    alive_chan[crp, view, vch, : ] = False

"""
#could be needed for comparisons
ped_ref = np.zeros((2, cf.n_View, cf.n_ChanPerCRP)) #crp, view, vchan
for idaq in range(cf.n_ChanTot):
crp, view, vch = cmap.DAQToCRP(idaq)    
if(crp > 1): continue
ped_ref[crp, view, vch] = ped.GetPedRMS(idaq)
"""


for ievent in range(nevent):
    print("-*-*-*-*-*-*-*-*-*-*-")
    print(" READING EVENT ", ievent)
    print("-*-*-*-*-*-*-*-*-*-*-")


    tevtread = time.time()
    idx = event_pos[ievent]
    evt, npdatav0, npdatav1 = read.read_event( data, idx)
    evt.evt_nb_loc = ievent
    
    print("RUN ",evt.run_nb, " EVENT ", evt.evt_nb_loc, " / ", evt.evt_nb_glob,)
    print("Taken at ", time.ctime(evt.time_s), " + ", evt.time_ns, " ns ")


    tevtdata = time.time()

    print(" -> Reading time %.2f s"%( tevtdata - tevtread))

    if( len(npdatav0)/cf.n_Sample != cf.n_ChanPerView):
            print(" PBM OF Nb of CHANNELS in V0 !!! ", len(npdatav0)/cf.n_Sample , " vs ", cf.n_ChanPerView)
    if( len(npdatav1)/cf.n_Sample != cf.n_ChanPerView):
            print(" PBM OF Nb of CHANNELS in V1 !!! ", len(npdatav1)/cf.n_Sample , " vs ", cf.n_ChanPerView)
    

    npdatav0 = np.split(npdatav0, cf.n_ChanPerView)
    npdatav1 = np.split(npdatav1, cf.n_ChanPerView)

        
    """reset"""
    npalldata[:,:,:,:] = 0.
    mask[:,:,:,:] = True


    """reshape the array and subtract reference pedestal"""
    """for noise processing, the reshaping is useless 
    but makes the code readable - change in the future ? """
    
    for idq in range(cf.n_ChanTot):
        crp, view, vch = cmap.DAQToCRP(idq)
        if(crp > 1): continue #Do not care about CRP 2 & 3 ATM
        pedval = ped.GetPed(idq)
        if(crp < 0 or view < 0 or vch < 0):
            print(" ERROR ? ", idq)
        if(view==0):
            npalldata[crp,view,vch] = npdatav0[idq] - pedval
        elif(view==1):
            npalldata[crp,view,vch] = npdatav1[idq-3840] - pedval
    

        if(run_nb <= cf.run_inv_signal):
            npalldata *= -1.    
        
    print(" done getting event ", ievent, " ! %.2f"%(time.time() - tevtdata))


    tfft = time.time()
    npalldata = noise.FFTLowPass(npalldata, lowpasscut, freqlines)    
    print(" time to fft %.2f"%( time.time() - tfft))
    
    """ 1st ROI attempt based on ADC cut + broken channels """
    mask = np.where( (npalldata > adc_thresh) | ~alive_chan, False, True)
    ped_ini = noise.get_RMS(npalldata*mask)

    """ Update ROI based on ped rms """
    troi = time.time()
    mask = noise.define_ROI(npalldata, mask, signal_thresh, 2)
    
    t3 = time.time()
    print(" 1st ROi : %.2f"% (t3-troi))

    """Apply coherent filter(s) """
    npalldata = noise.coherent_filter(npalldata, mask, coherent_group)
    npalldata = noise.coherent_filter(npalldata, mask, 320)
    
    print(" time to coh filt %.2f"%( time.time() - t3))


    t4 = time.time()

    """ Update ROI regions """
    mask = noise.define_ROI(npalldata, mask, signal_thresh, 2)

    print(" time to ROI %.2f"%(time.time() - t4))


    """invert the mask, so now True is signal (and not broken channels)"""    
    ROI = np.array(~mask & alive_chan, dtype=bool)
    print("Nb of points ", np.sum(ROI))
    #only for plotting purpose
    #ROI *= npalldata     
    #ROI = clus.rebin(ROI, 2, 5)
    """
    t5 = time.time()
    hf.HitSearch(ROI[0,0,300], 2.)
    print(" time to Hit Search %.2f"%(time.time() - t5))
    """
    

    #t6 = time.time()
    #hf.HitFinder(npalldata, ROI, noise.get_RMS(npalldata*mask))
    #print(" time to Hit Search %.3f"%(time.time() - t6))
    #print(" cross check ", len(cf.hits_list))
    #dec.deconvolute(npalldata[0,0])
    plot_ev.plot_event_display(npalldata, evt.run_nb, evt.evt_nb_glob,"filt")
    #plot_ev.plot_event_display(ROI, evt.run_nb, evt.evt_nb_glob, "ROI_rb")
    """
    ped_fin = noise.get_RMS(npalldata*mask)
    plot_ev.plot_pedestal([ped_ref, ped_ini, ped_fin], ['Ref.', 'Raw', 'Final'],['r','k','b'], evt.run_nb, evt.evt_nb_glob)
    """

data.close()
tottime = time.time() - tstart
print(" TOTAL RUNNING TIME %.2f s == %.2f evt/s"% (tottime, tottime/nevent))
