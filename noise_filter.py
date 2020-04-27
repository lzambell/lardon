import config as cf
import data_containers as dc

import numpy as np
import numexpr as ne 
import numba as nb
import time

def define_ROI_ADC(thresh):
    #dc.mask = np.where( (dc.data > thresh) | ~dc.alive_chan, False, True)
    dc.mask = ne.evaluate( "where((data > thresh) | ~alive_chan, 0, 1)", global_dict={'data':dc.data, 'alive_chan':dc.alive_chan}).astype(bool)


def define_ROI(sig_thresh, iteration):
    #ne.set_num_threads(4) #does not speed up things

    """ Update the mask based on pedestal RMS """    
    for it in range(iteration):
        compute_pedestal_RMS()

        dc.ped_rms = dc.ped_rms[:,:,:,None]
        dc.mask = ne.evaluate( "where((data > sig_thresh*rms) | (~mask), 0, 1)", global_dict={'data':dc.data, 'mask':dc.mask, 'rms':dc.ped_rms}).astype(bool)
        dc.ped_rms = np.squeeze(dc.ped_rms, axis=3)


def compute_pedestal_RMS():
    dc.ped_rms =  np.std(dc.data * dc.mask, axis=3)

def compute_pedestal_mean():
    dc.ped_mean = np.mean(dc.data * dc.mask, axis=3)

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
    
def FFT2D() :
    
	for icrp in range(2):
		for iview in range(2):
		
			"""go to the 2D frequency domain"""
			fft2D=np.fft.fft2(dc.data[icrp,iview,:,:])
		
    
			'''x = np.linspace(0, 10000, 10000 )
			y = np.linspace(0, 960, 960)
			X, Y = np.meshgrid(x, y)'''
			gmask = np.ones((960,10000))
			
			"""Low Pass filter"""
			for i in range(960):
				gmask[i][400:9600]=0.
			t=time.time()
			"""Removing specific frequencies (different for each view and crp)"""
			if icrp==0 and iview==0:			
				for i in range(5):
					gmask[i][20:9980]=0
					gmask[959-i][20:9980]=0
				for i in range(13):
            				gmask[i][86:100]=0
            				gmask[i][9900:9913]=0
            				gmask[959-i][86:100]=0
            				gmask[959-i][9900:9913]=0
				for i in range(40):
            				gmask[i][91:96]=0
            				gmask[i][9905:9909]=0
            				gmask[959-i][91:96]=0
            				gmask[959-i][9905:9909]=0
            				gmask[i][249:251]=0
            				gmask[i][9749:9751]=0
            				gmask[959-i][249:251]=0
            				gmask[959-i][9749:9751]=0
				for i in range(100):
            				gmask[i][280:282]=0
            				gmask[i][9718:9721]=0
            				gmask[959-i][280:282]=0
            				gmask[959-i][9718:9721]=0
				for i in range(40,920):
            				gmask[i][92:95]=0
            				gmask[i][9906:9908]=0
				for i in range(3):
            				gmask[i][:]=0
            				gmask[959-i][:]=0
			
			if icrp==0 and iview==1:
				for i in range(5):
            				gmask[i][20:9980]=0
            				gmask[959-i][20:9980]=0
				for i in range(30):
					gmask[i][92:95]=0
					gmask[i][9906:9908]=0
					gmask[959-i][92:95]=0
					gmask[959-i][9906:9908]=0
				for i in range(2):
					gmask[i][:]=0
					gmask[959-i][:]=0
			if icrp==1 and iview==0:
				for i in range(5):
					gmask[i][20:9980]=0
					gmask[959-i][20:9980]=0
				for i in range(30):
					gmask[i][92:95]=0
					gmask[i][9906:9908]=0
					gmask[959-i][92:95]=0
					gmask[959-i][9906:9908]=0
				for i in range(2):
					gmask[i][:]=0
					gmask[959-i][:]=0
			if icrp==1 and iview==1:
				for i in range(5):
					gmask[i][20:9980]=0
					gmask[959-i][20:9980]=0
				for i in range(30):
					gmask[i][92:95]=0
					gmask[i][9906:9908]=0
					gmask[959-i][92:95]=0
					gmask[959-i][9906:9908]=0
				for i in range(2):
					gmask[i][:]=0
					gmask[959-i][:]=0
			
		
			#gmask = np.where(gmask<0,0,gmask)
			"""Apply the cuts"""
			fft2D = fft2D * gmask
			"""Go back in real space"""
			filt = np.fft.ifft2(fft2D).real
			dc.data[icrp,iview,:,:] = filt
			print("Time to FFT2D: "+str(time.time()-t))
		
	#return data #fft2D to see 2D frequency domain
