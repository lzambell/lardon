import config as cf
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import itertools as itr

    
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
        print(" ooops, I will only plot 9 waveforms")
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
    plt.close()    



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
    plt.close()



def plot_pedestal(datas, legtitle, colors, run_nb, evt_nb):
    nplot = len(datas)

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
    plt.close()



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
    plt.close()


    
def plot_hits_ed(ncluster, run_nb, evt_nb, option=None):
    fig = plt.figure(figsize=(12,9))
    ax  = []
    iplot=0

    for iview in range(2):
        for icrp in range(2):
            iplot += 1
            ax.append(fig.add_subplot(2,2,iplot))
            
            hit_pos = [x.channel for x in cf.hits_list if x.view == iview and x.crp == icrp]
            hit_tdc = [x.max_t for x in cf.hits_list if x.view == iview and x.crp == icrp]
            hit_cls = [x.cluster for x in cf.hits_list if x.view == iview and x.crp == icrp]

            

            colors = np.array(list(itr.islice(itr.cycle(['#377eb8', '#ff7f00', '#4daf4a',
                                                 '#f781bf', '#a65628', '#984ea3',
                                                 '#79f2ff', '#e41a1c', '#dede00']),
                                          int(ncluster[icrp,iview] + 1))))
            # add grey color for outliers (if any)
            colors = np.append(colors, ["#c7c7c7"])

            
            ax[-1].scatter(hit_pos, hit_tdc, c=colors[hit_cls],s=2)
            
            ax[-1].set_xlim(0,959)
            ax[-1].set_ylim(0,9999)
            ax[-1].set_xlabel('View Channel')
            ax[-1].set_ylabel('Time')
            ax[-1].set_title('CRP '+str(icrp)+' - View '+str(iview))    
            
            print("CRP ", icrp, " View ", iview, " -> Nhits : ", len(hit_tdc))

    plt.subplots_adjust(bottom=0.08, top=0.95, hspace=0.3)
    if(option):
        option = "_"+option
    else:
        option = ""
    plt.savefig('ED/hit'+option+'_run_'+str(run_nb)+'_evt_'+str(evt_nb)+'.png')
    #plt.show()
    plt.close()


def plot_hits_var(run_nb,evt_nb,option=None):
    fig = plt.figure(figsize=(12,9))
    ax  = []
    
    
    ax.append(fig.add_subplot(3,2,1))
    charge = [x.charge for x in cf.hits_list]
    ax[-1].hist(charge, 200, color='orange', histtype="stepfilled", edgecolor='r',log=True,range=(0., 5000.))
    ax[-1].set_xlabel('Hit charge [ADC]')

    ax.append(fig.add_subplot(3,2,2))
    tmax = [x.max_t for x in cf.hits_list]
    ax[-1].hist(tmax, 250, color='c',histtype='stepfilled',edgecolor='b',range=(0, 10000))
    ax[-1].set_xlabel('Hit Time [tdc]')

    ax.append(fig.add_subplot(3,2,3))
    dt = [x.stop-x.start for x in cf.hits_list]
    ax[-1].hist(dt, 200, color='b',histtype='stepfilled',edgecolor='k',log=True, range=(0., 200.))
    ax[-1].set_xlabel('Hit Length [tdc]')

    ax.append(fig.add_subplot(3,2,4))
    amp = [x.max_adc for x in cf.hits_list]
    ax[-1].hist(amp, 200, color='r',histtype='stepfilled',edgecolor='k',log=True, range=(0., 60.))
    ax[-1].set_xlabel('Hit Amplitude [ADC]')

    ax.append(fig.add_subplot(3,2,5))
    dtstart = [x.max_t-x.start for x in cf.hits_list]
    ax[-1].hist(dtstart, 50, color='y',histtype='stepfilled',edgecolor='k',log=True,range=(0., 50))
    ax[-1].set_xlabel('Hit start-max [tdc]')

    ax.append(fig.add_subplot(3,2,6))
    dtstop = [x.stop-x.max_t for x in cf.hits_list]
    ax[-1].hist(dtstop, 50, color='g',histtype='stepfilled',edgecolor='k',log=True,range=(0.,50.))
    ax[-1].set_xlabel('Hit max-stop [tdc]')


    plt.subplots_adjust(bottom=0.08, top=0.95, hspace=0.3)
    if(option):
        option = "_"+option
    else:
        option = ""
    plt.savefig('ED/hit_properties_'+option+'_run_'+str(run_nb)+'_evt_'+str(evt_nb)+'.png')
    #plt.show()
    plt.close()
    
