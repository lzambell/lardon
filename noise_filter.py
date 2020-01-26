import config as cf
import numpy as np
import numpy.ma as ma



def get_RMS(data):
    # Assume the data are organized as np array of (crp, view, anachan, tdc)
    # might change later
    return np.std(data, axis=3)



def coherent_filter_prev(data, group):
    # Assume the data are organized as np array of (crp, view, anachan, tdc)
    # might change later

    if( (cf.n_ChanPerCRP % group) > 0):
        print " Coherent Noise Filter in groups of ", group, " is not a possible ! "
        return

    nslices = cf.n_ChanPerCRP / group

    data = np.reshape(data, (2, cf.n_View, group, nslices, cf.n_Sample))
    mean = np.mean(data, axis=3)

    data.mask = ma.nomask #Apply the correction to all data points
    data -= mean[:,:,:,None,:]
    data = np.reshape(data, (2, cf.n_View, cf.n_ChanPerCRP, cf.n_Sample))

    return data



def coherent_filter(data, alive, group, adc_thresh, sig_thresh):
    # Assume the data are organized as np array of (crp, view, anachan, tdc)
    # might change later

    if( (cf.n_ChanPerCRP % group) > 0):
        print " Coherent Noise Filter in groups of ", group, " is not a possible ! "
        return

    nslices = cf.n_ChanPerCRP / group
    mask = np.ones((2,cf.n_View,cf.n_ChanPerCRP,cf.n_Sample), dtype=bool)


    totlength = 2*2*960*10000
    #print " initial mask has ", totlength - np.count_nonzero(mask), " signal-like "

    mask = np.where( (data > adc_thresh) | ~alive, False, True)
    #print " after adc thrsh OR alive, mask has ", totlength - np.count_nonzero(mask), " signal-like "

    for it in range(2):
        rms = get_RMS(data*mask)
        mask = np.where( (data > sig_thresh*rms[:,:,:,None]) | ~mask, False, True)

        #print " after ", it, " iteration, mask has ", totlength - np.count_nonzero(mask), " signal-like "  

    data = np.reshape(data, (2, cf.n_View, group, nslices, cf.n_Sample))
    mask = np.reshape(mask, (2, cf.n_View, group, nslices, cf.n_Sample))


    
    #sum data if mask is true
    with np.errstate(divide='ignore', invalid='ignore'):
        mean = np.einsum('ijklm,ijklm->ijkm',data,mask)/mask.sum(axis=3)
        #mean[mask.sum(axis=3)==0] = 0.
        mean[mask.sum(axis=3) < 3] = 0.
        #require at least 3 points to taken into account the mean


    
    #Apply the correction to all data points
    data -= mean[:,:,:,None,:]

    data = np.reshape(data, (2, cf.n_View, cf.n_ChanPerCRP, cf.n_Sample))
    mask = np.reshape(data, (2, cf.n_View, cf.n_ChanPerCRP, cf.n_Sample))

    return data, mask



def gaussian(x, mu, sig):
    return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))

def FFTLowPass(data, lowpass, freqlines) :

    #it's 5001 (n/2+1) points from 0Hz to 1./(2*sampling) = 1.25MHz (nyquist freq)
    
    n    = cf.n_Sample/2 + 1
    rate = 1./cf.n_Sampling #in MHz
    freq = np.linspace(0, rate/2., n)

    #define gaussian low pass filter
    gauss_cut = np.where(freq < lowpass, 1., gaussian(freq, lowpass, 0.02))
    
    for f in freqlines:
        fbin = f * cf.n_Sample * cf.n_Sampling
        gauss_cut[fbin-2:fbin+3] = 0.1
        gauss_cut[fbin-1:fbin+2] = 0.

    #go to frequency domain
    fdata = np.fft.rfft(data)
    
    #Apply filter
    fdata *= gauss_cut[None, None, None, :]

    #go back to time
    data = np.fft.irfft(fdata)

    #get power spectrum
    #ps = 10.*np.log10(np.abs(fdata)) 
    return data #,ps
