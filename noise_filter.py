import config as cf
import data_containers as dc

import numpy as np
import numexpr as ne 

def define_ROI_ADC(thresh):
    #dc.mask = np.where( (dc.data > thresh) | ~dc.alive_chan, False, True)
    dc.mask = ne.evaluate( "where((data > thresh) | ~alive_chan, 0, 1)", global_dict={'data':dc.data, 'alive_chan':dc.alive_chan}).astype(bool)

def define_ROI(sig_thresh, iteration):
    #ne.set_num_threads(4)
    """ Update the mask based on pedestal RMS """    
    for it in range(iteration):

        rms = get_RMS()
        rms = rms[:,:,:,None]
        dc.mask = ne.evaluate( "where((data > sig_thresh*rms) | (~mask), 0, 1)", global_dict={'data':dc.data, 'mask':dc.mask}).astype(bool)



def get_RMS():
    return np.std(dc.data * dc.mask, axis=3)


def coherent_filter(groupings):
    """
    1. Computes the mean along group of channels for non ROI points
    2. Subtract mean to all points
    """

    for group in groupings:
        if( (cf.n_ChanPerCRP % group) > 0):
            print(" Coherent Noise Filter in groups of ", group, " is not a possible ! ")
            return

        nslices = int(cf.n_ChanPerCRP / group)
        
        dc.data = np.reshape(dc.data, (cf.n_CRPUsed, cf.n_View, group, nslices, cf.n_Sample))
        dc.mask = np.reshape(dc.mask, (cf.n_CRPUsed, cf.n_View, group, nslices, cf.n_Sample))


    
        """sum data if mask is true"""
        with np.errstate(divide='ignore', invalid='ignore'):
            mean = np.einsum('ijklm,ijklm->ijkm', dc.data, dc.mask)/dc.mask.sum(axis=3)

            """require at least 3 points to take into account the mean"""
            mean[dc.mask.sum(axis=3) < 3] = 0.
        

        """Apply the correction to all data points"""
        dc.data -= mean[:,:,:,None,:]


        """ restore original data shape """
        dc.data = np.reshape(dc.data, (cf.n_CRPUsed, cf.n_View, cf.n_ChanPerCRP, cf.n_Sample))
        dc.mask = np.reshape(dc.mask, (cf.n_CRPUsed, cf.n_View, cf.n_ChanPerCRP, cf.n_Sample))




def gaussian(x, mu, sig):
    return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))


def FFTLowPass(lowpass, freqlines) :

    """it's 5001 (n/2+1) points from 0Hz to 1./(2*sampling) = 1.25MHz (nyquist freq)"""
    
    n    = int(cf.n_Sample/2) + 1
    rate = 1./cf.n_Sampling #in MHz
    freq = np.linspace(0, rate/2., n)

    """define gaussian low pass filter"""
    gauss_cut = np.where(freq < lowpass, 1., gaussian(freq, lowpass, 0.02))
    
    """remove specific ferquencies"""
    for f in freqlines:
        fbin = int(f * cf.n_Sample * cf.n_Sampling)
        #print "frequency ", f, " at bin ", fbin
        gauss_cut[max(fbin-2,0):min(fbin+3,n)] = 0.1
        gauss_cut[max(fbin-1,0):min(fbin+2,n)] = 0.

    """go to frequency domain"""
    fdata = np.fft.rfft(dc.data)
    
    """Apply filter"""
    fdata *= gauss_cut[None, None, None, :]

    """go back to time"""
    dc.data = np.fft.irfft(fdata)

    """get power spectrum"""
    #ps = 10.*np.log10(np.abs(fdata)) 
    #return ps
