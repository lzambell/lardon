import config as cf
import numpy as np
import numexpr as ne 


def define_ROI(data, mask, sig_thresh, iteration):
    #ne.set_num_threads(4)
    """ Update the mask based on pedestal RMS """    
    for it in range(iteration):
        rms = get_RMS(data*mask)
        #mask = np.where( (data > sig_thresh*rms[:,:,:,None]) | (~mask), False, True)    
        rms = rms[:,:,:,None]
        mask = ne.evaluate( "where((data > sig_thresh*rms) | (~mask), 0, 1)" ).astype(bool)
    return mask


def get_RMS(data):
    # Assume the data are organized as np array of (crp, view, anachan, tdc)
    # might change later
    return np.std(data, axis=3)



def coherent_filter(data, mask, group):
    """
    Assume the data are organized as np array of (crp, view, anachan, tdc)
    (might change later)
    1. Computes the mean along group of channels for non ROI points
    2. Subtract mean to all points
    """

    if( (cf.n_ChanPerCRP % group) > 0):
        print(" Coherent Noise Filter in groups of ", group, " is not a possible ! ")
        return

    nslices = int(cf.n_ChanPerCRP / group)
        
    data = np.reshape(data, (2, cf.n_View, group, nslices, cf.n_Sample))
    mask = np.reshape(mask, (2, cf.n_View, group, nslices, cf.n_Sample))


    
    #sum data if mask is true
    with np.errstate(divide='ignore', invalid='ignore'):
        mean = np.einsum('ijklm,ijklm->ijkm',data,mask)/mask.sum(axis=3)

        """require at least 3 points to taken into account the mean"""
        mean[mask.sum(axis=3) < 3] = 0.
        

    """Apply the correction to all data points"""
    data -= mean[:,:,:,None,:]

    data = np.reshape(data, (2, cf.n_View, cf.n_ChanPerCRP, cf.n_Sample))

    return data



def gaussian(x, mu, sig):
    return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))

def FFTLowPass(data, lowpass, freqlines) :

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
    fdata = np.fft.rfft(data)
    
    """Apply filter"""
    fdata *= gauss_cut[None, None, None, :]

    """go back to time"""
    data = np.fft.irfft(fdata)

    """get power spectrum"""
    #ps = 10.*np.log10(np.abs(fdata)) 
    return data #,ps
