import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap


cdict1 = {
    'red': ((0.,    65./255.,  65./255.),
            (0.15, 123./255., 123./255.),
            (0.25, 160./255., 160./255.),
            (0.375, 222./255.,222./255.),
            (0.5, 214./255., 214./255.),
            (0.625, 199./255., 199./255.),
            (0.75, 183./255., 183./255.),
            (0.875, 153./255., 153./255.),
            (1., 78./255., 78./255.)),

    'green':  ((0.,  90./255.,  90./255.),
              (0.15, 171./255., 171./255.),
              (0.25,  211./255., 211./255.),
              (0.375,  220./255.,  220./255.),
              (0.5, 190./255., 190./255.),
              (0.625,  132./255., 132./255.),
              (0.75,  65./255.,  65./255.),
              (0.875, 0./255., 0./255.),
               (1.,  0./255., 0./255.)),
        
    
    
    'blue':   ((0.,  148./255., 148./255.),
               (0.15, 228./255., 228./255.),
               (0.25, 222./255., 222./255.),
               (0.375,  160./255.,  160./255.),
               (0.5, 105./255., 105./255.),
               (0.625, 60./255., 60./255.),
               (0.75, 34./255., 34./255.),
               (0.875, 0./255., 0./255.),
               (1.,  0./255., 0./255.))
    
}
newcmp = LinearSegmentedColormap('MyOwn', cdict1)

adcmin = -10
adcmax = 70


def plot_waveform(data, legtitle, colors, run_nb, evt_nb):
    nplot = len(data)
    if(nplot > 9):
        print " ooops, I will only plot 9 waveforms"
    fig = plt.figure(figsize=(12,9))
    ax = []
    for ip in range(nplot):
        ax.append(fig.add_subplot(nplot,1,ip+1))
        d = data[ip]                
        plt.plot(d, colors[ip], label=legtitle[ip])
        ax[-1].set_xlabel('Time [tdc]')
        ax[-1].set_ylabel('ADC')
        ax[-1].set_ylim(-8., 35.)
        plt.legend()
        plt.subplots_adjust(bottom=0.05, top=.98, hspace=0.27)
    plt.savefig('ED/waveform_run_'+str(run_nb)+'_evt_'+str(evt_nb)+'.png')
    



def plot_event_display(data, run_nb, evt_nb, option=None):    
    fig = plt.figure(figsize=(12,9))
    ax = []
    im = []
    iplot = 0
    for icrp in range(2):
        for iview in range(2):
            iplot = iplot + 1
            ax.append(fig.add_subplot(2,2,iplot))
            im.append(plt.imshow(data[icrp,iview,:,:].transpose(), origin='lower',aspect='auto',cmap=newcmp, vmin=adcmin, vmax=adcmax))
            plt.colorbar(im[-1])
            ax[-1].set_xlabel('View Channel')
            ax[-1].set_ylabel('Time')
            ax[-1].set_title('CRP '+str(icrp)+' - View '+str(iview))    
    plt.subplots_adjust(bottom=0.08, top=0.95)
    if(option):
        option = "_"+option
    plt.savefig('ED/ed'+option+'_run_'+str(run_nb)+'_evt_'+str(evt_nb)+'.png')
    #plt.show()




def plot_pedestal(datas, legtitle, colors, run_nb, evt_nb):
    nplot = len(datas)
    #colors = ['k','b','r']
    #legtitle = ['Raw', 'FFT', 'Coherent']
    fig = plt.figure(figsize=(12,9))
    ax = []
    iplot = 0
    for icrp in range(2):
        for iview in range(2):
            iplot = iplot + 1
            ax.append(fig.add_subplot(2,2,iplot))
            for ip in range(nplot):
                d = datas[ip]                
                plt.plot(d[icrp, iview,:], colors[ip], label=legtitle[ip])
            ax[-1].set_xlabel('View Channel')
            ax[-1].set_ylabel('Pedestal RMS (ADC)')
            ax[-1].set_title('CRP '+str(icrp)+' - View '+str(iview))  
            ax[-1].set_ylim(0., 6.)
            plt.legend()
    plt.subplots_adjust(bottom=0.08, top=0.95)
    plt.savefig('ED/pedrms_run_'+str(run_nb)+'_evt_'+str(evt_nb)+'.png')
    #plt.show()




def plot_event_fft(data, run_nb, evt_nb):    

    fig = plt.figure(figsize=(12,9))
    ax = []
    im = []
    iplot = 0


    for icrp in range(2):
        for iview in range(2):
            iplot += 1
            ax.append(fig.add_subplot(2,2,iplot))
            im.append(plt.imshow(data[icrp,iview,:,:].transpose(), origin='lower',aspect='auto',cmap=newcmp, extent=[0,959,0.,1.25], vmin=0.))
            plt.colorbar(im[-1])
            plt.ylim(0.,0.25)
            #plt.ylim(0.06, 0.08)
            ax[-1].set_xlabel('View Channel')
            ax[-1].set_ylabel('Signal Frequencies [MHz]')
            ax[-1].set_title('CRP '+str(icrp)+' - View '+str(iview))    

    plt.subplots_adjust(bottom=0.08, top=0.95)
    plt.savefig('ED/fft_run_'+str(run_nb)+'_evt_'+str(evt_nb)+'.png')
    #plt.show()



    
