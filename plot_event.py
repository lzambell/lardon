import config as cf
import data_containers as dc


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.gridspec as gridspec
from matplotlib import collections  as mc
import itertools as itr
import math

    
light_blue_red_dict = {
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
lbr_cmp = LinearSegmentedColormap('lightBR', light_blue_red_dict)

adcmin = -10
adcmax = 35


def plot_waveform(data, legtitle, colors, option=None):

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



    if(option):
        option = "_"+option
    else:
        option = ""


    run_nb = str(dc.evt_list[-1].run_nb)
    evt_nb = str(dc.evt_list[-1].evt_nb_glob)


    plt.savefig('ED/waveform'+option+'_run_'+run_nb+'_evt_'+evt_nb+'.png')
    plt.close()    



def plot_event_display(option=None):    

    fig = plt.figure(figsize=(12,9))
    ax = []
    im = []
    iplot = 0

    for icrp in range(2):
        for iview in range(2):
            iplot = iplot + 1
            ax.append(fig.add_subplot(2,2,iplot))
            im.append(plt.imshow(dc.data[icrp,iview,:,:].transpose(), origin='lower',aspect='auto',cmap=lbr_cmp, vmin=adcmin, vmax=adcmax))
            plt.colorbar(im[-1])

            ax[-1].set_xlabel('View Channel')
            ax[-1].set_ylabel('Time')
            ax[-1].set_title('CRP '+str(icrp)+' - View '+str(iview))    

    plt.subplots_adjust(bottom=0.08, top=0.95)

    if(option):
        option = "_"+option
    else:
        option = ""

    run_nb = str(dc.evt_list[-1].run_nb)
    evt_nb = str(dc.evt_list[-1].evt_nb_glob)

    plt.savefig('ED/ed'+option+'_run_'+run_nb+'_evt_'+evt_nb+'.png')

    #plt.show()
    plt.close()



def plot_pedestal(datas, legtitle, colors, option=None):

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


    if(option):
        option = "_"+option
    else:
        option = ""

    run_nb = str(dc.evt_list[-1].run_nb)
    evt_nb = str(dc.evt_list[-1].evt_nb_glob)


    plt.savefig('ED/pedrms'+option+'_run_'+run_nb+'_evt_'+evt_nb+'.png')
    #plt.show()
    plt.close()



def plot_event_fft(data):    

    fig = plt.figure(figsize=(12,9))
    ax = []
    im = []
    iplot = 0


    for icrp in range(2):
        for iview in range(2):
            iplot += 1

            ax.append(fig.add_subplot(2,2,iplot))

            im.append(plt.imshow(data[icrp,iview,:,:].transpose(), origin='lower',aspect='auto',cmap=lbr_cmp, extent=[0, 959, 0., 1.25], vmin=0.))

            plt.colorbar(im[-1])

            plt.ylim(0.,0.25)
            #plt.ylim(0.06, 0.08)

            ax[-1].set_xlabel('View Channel')
            ax[-1].set_ylabel('Signal Frequencies [MHz]')
            ax[-1].set_title('CRP '+str(icrp)+' - View '+str(iview))    

    plt.subplots_adjust(bottom=0.08, top=0.95)


    if(option):
        option = "_"+option
    else:
        option = ""

    run_nb = str(dc.evt_list[-1].run_nb)
    evt_nb = str(dc.evt_list[-1].evt_nb_glob)



    plt.savefig('ED/fft_run_'+run_nb+'_evt_'+evt_nb+'.png')
    #plt.show()
    plt.close()


    
def plot_hits_clustered(option=None):

    fig = plt.figure(figsize=(12,9))
    ax  = []
    iplot=0

    for icrp in range(2):
        for iview in range(2):
            iplot += 1
            ax.append(fig.add_subplot(2,2,iplot))
            
            hit_pos = [x.channel for x in dc.hits_list if x.view == iview and x.crp == icrp]
            hit_tdc = [x.max_t for x in dc.hits_list if x.view == iview and x.crp == icrp]
            hit_cls = [x.cluster for x in dc.hits_list if x.view == iview and x.crp == icrp]

            

            colors = np.array(list(itr.islice(itr.cycle(['#377eb8', '#ff7f00', '#4daf4a',
                                                         '#f781bf', '#a65628', '#984ea3',
                                                         '#79f2ff', '#e41a1c', '#dede00']),
                                              int(dc.ncluster[icrp,iview] + 1))))

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


    run_nb = str(dc.evt_list[-1].run_nb)
    evt_nb = str(dc.evt_list[-1].evt_nb_glob)


    plt.savefig('ED/hit'+option+'_run_'+run_nb+'_evt_'+evt_nb+'.png')
    #plt.show()
    plt.close()


def plot_hits_var(option=None):

    fig = plt.figure(figsize=(12,9))
    ax  = []
    
    
    ax.append(fig.add_subplot(3,2,1))
    charge = [x.charge for x in dc.hits_list]
    ax[-1].hist(charge, 200, color='orange', histtype="stepfilled", edgecolor='r',log=True,range=(0., 5000.))
    ax[-1].set_xlabel('Hit charge [ADC]')

    ax.append(fig.add_subplot(3,2,2))
    tmax = [x.max_t for x in dc.hits_list]
    ax[-1].hist(tmax, 250, color='c',histtype='stepfilled',edgecolor='b',range=(0, 10000))
    ax[-1].set_xlabel('Hit Time [tdc]')

    ax.append(fig.add_subplot(3,2,3))
    dt = [x.stop-x.start for x in dc.hits_list]
    ax[-1].hist(dt, 200, color='b',histtype='stepfilled',edgecolor='k',log=True, range=(0., 200.))
    ax[-1].set_xlabel('Hit Length [tdc]')

    ax.append(fig.add_subplot(3,2,4))
    amp = [x.max_adc for x in dc.hits_list]
    ax[-1].hist(amp, 200, color='r',histtype='stepfilled',edgecolor='k',log=True, range=(0., 60.))
    ax[-1].set_xlabel('Hit Amplitude [ADC]')

    ax.append(fig.add_subplot(3,2,5))
    dtstart = [x.max_t-x.start for x in dc.hits_list]
    ax[-1].hist(dtstart, 50, color='y',histtype='stepfilled',edgecolor='k',log=True,range=(0., 50))
    ax[-1].set_xlabel('Hit start-max [tdc]')

    ax.append(fig.add_subplot(3,2,6))
    dtstop = [x.stop-x.max_t for x in dc.hits_list]
    ax[-1].hist(dtstop, 50, color='g',histtype='stepfilled',edgecolor='k',log=True,range=(0.,50.))
    ax[-1].set_xlabel('Hit max-stop [tdc]')


    plt.subplots_adjust(bottom=0.08, top=0.95, hspace=0.3)
    if(option):
        option = "_"+option
    else:
        option = ""


    run_nb = str(dc.evt_list[-1].run_nb)
    evt_nb = str(dc.evt_list[-1].evt_nb_glob)

    plt.savefig('ED/hit_properties_'+option+'_run_'+run_nb+'_evt_'+evt_nb+'.png')
    #plt.show()
    plt.close()
    



def plot_hits_view(option=None):

    fig = plt.figure(figsize=(12,6))
    gs = gridspec.GridSpec(nrows=2, ncols=3, height_ratios=[1,10])

    ax_col = fig.add_subplot(gs[0,:])
    ax_v0 = fig.add_subplot(gs[1, 0:2])
    ax_v1 = fig.add_subplot(gs[1, 2], sharey=ax_v0)

    max_adc_range=50
    hit_x_v0 = [x.X for x in dc.hits_list if x.view == 0]
    hit_z_v0 = [x.Z for x in dc.hits_list if x.view == 0]
    hit_q_v0 = [x.max_adc for x in dc.hits_list if x.view == 0]


    hit_x_v1 = [x.X for x in dc.hits_list if x.view == 1]
    hit_z_v1 = [x.Z for x in dc.hits_list if x.view == 1]
    hit_q_v1 = [x.max_adc for x in dc.hits_list if x.view == 1]

    sc0 = ax_v0.scatter(hit_x_v0, hit_z_v0, c=hit_q_v0, cmap=lbr_cmp, s=2, vmin=0, vmax=max_adc_range)

    ax_v0.set_title('View 0')
    ax_v0.set_ylabel('Z [cm]')
    ax_v0.set_xlabel('Y [cm]')
    ax_v0.set_xlim([-300., 300])
    ax_v0.set_ylim([-30., 300])

    sc1 = ax_v1.scatter(hit_x_v1, hit_z_v1, c=hit_q_v1, cmap=lbr_cmp, s=2,vmin=0,vmax=max_adc_range)
    ax_v1.set_title('View 1')
    #ax_v1.set_ylabel('Z [cm]')
    ax_v1.set_xlabel('X [cm]')
    ax_v1.set_xlim([0,300])
    ax_v1.set_ylim([-30., 300])
    ax_col.set_title('Hit Max ADC')
    cb = fig.colorbar(sc1, cax=ax_col, orientation='horizontal')
    cb.ax.xaxis.set_ticks_position('top')
    cb.ax.xaxis.set_label_position('top')
    
    if(option):
        option = "_"+option
    else:
        option = ""


    run_nb = str(dc.evt_list[-1].run_nb)
    evt_nb = str(dc.evt_list[-1].evt_nb_glob)


    plt.savefig('ED/hit_view'+option+'_run_'+run_nb+'_evt_'+evt_nb+'.png')
    #plt.show()
    plt.close()

def plot_tracks2D(option=None):

    fig = plt.figure(figsize=(12,6))
    gs = gridspec.GridSpec(nrows=1, ncols=3)#, height_ratios=[1,10])

    #ax_col = fig.add_subplot(gs[0,:])
    ax_v0 = fig.add_subplot(gs[0, 0:2])
    ax_v1 = fig.add_subplot(gs[0, 2], sharey=ax_v0)

    #max_adc_range=50
    hit_x_v0 = [x.X for x in dc.hits_list if x.view == 0]
    hit_z_v0 = [x.Z for x in dc.hits_list if x.view == 0]

    hit_x_v0_noise = [x.X for x in dc.hits_list if x.view == 0 and x.cluster==-1]
    hit_z_v0_noise = [x.Z for x in dc.hits_list if x.view == 0 and x.cluster==-1]

    tracks_hits_x_v0 = [[p[0] for p in t.path] for t in dc.tracks2D_list if t.view==0]
    tracks_hits_z_v0 = [[p[1] for p in t.path] for t in dc.tracks2D_list if t.view==0]

    hit_x_v1 = [x.X for x in dc.hits_list if x.view == 1]
    hit_z_v1 = [x.Z for x in dc.hits_list if x.view == 1]

    hit_x_v1_noise = [x.X for x in dc.hits_list if x.view == 1 and x.cluster==-1]
    hit_z_v1_noise = [x.Z for x in dc.hits_list if x.view == 1 and x.cluster==-1]

    tracks_hits_x_v1 = [[p[0] for p in t.path] for t in dc.tracks2D_list if t.view==1]
    tracks_hits_z_v1 = [[p[1] for p in t.path] for t in dc.tracks2D_list if t.view==1]



    ax_v0.scatter(hit_x_v0, hit_z_v0, c="#ffb77d", s=2)#ffb77d
    ax_v0.scatter(hit_x_v0_noise, hit_z_v0_noise, c="#d8cfd6", s=2)
    

    for tx,tz in zip(tracks_hits_x_v0, tracks_hits_z_v0):
        ax_v0.scatter(tx, tz, c="#28568f",s=2) ##6d40cf
        ax_v0.plot(tx,tz, c="#de425b",linewidth=1)#f65789

    ax_v0.set_title('View 0')
    ax_v0.set_ylabel('Z [cm]')
    ax_v0.set_xlabel('Y [cm]')
    ax_v0.set_xlim([-300., 300])
    ax_v0.set_ylim([-30., 300])

    ax_v1.scatter(hit_x_v1, hit_z_v1, c="#ffb77d", s=2)#ffb77d
    ax_v1.scatter(hit_x_v1_noise, hit_z_v1_noise, c="#d8cfd6", s=2)

    for tx,tz in zip(tracks_hits_x_v1, tracks_hits_z_v1):
        ax_v1.scatter(tx, tz, c="#28568f",s=2)#6d40cf
        ax_v1.plot(tx,tz, c="#de425b",linewidth=1)#f65789


    ax_v1.set_title('View 1')
    ax_v1.set_xlabel('X [cm]')
    ax_v1.set_xlim([0,300])
    ax_v1.set_ylim([-30., 300])
    
    if(option):
        option = "_"+option
    else:
        option = ""


    run_nb = str(dc.evt_list[-1].run_nb)
    evt_nb = str(dc.evt_list[-1].evt_nb_glob)


    plt.savefig('ED/track2D'+option+'_run_'+run_nb+'_evt_'+evt_nb+'.png')
    plt.show()
    plt.close()


def plot_track2D_var(option=None):
    
    fig = plt.figure(figsize=(12,6))
    gs = gridspec.GridSpec(nrows=3, ncols=2)

    ax_nhits = fig.add_subplot(gs[0, 0])
    ax_length = fig.add_subplot(gs[0, 1])
    ax_angles = fig.add_subplot(gs[1:, :])

    nb_hits_per_track = [x.nHits for x in dc.tracks2D_list]
    ax_nhits.hist(nb_hits_per_track, bins=100, range=(0,100), histtype="stepfilled", color="#9bad51", edgecolor='#68723d')
    ax_nhits.set_xlabel('Nb of Hits attached to tracks')

    track_length = [math.sqrt(pow(x.path[0][0]-x.path[-1][0], 2) + pow(x.path[0][1]-x.path[-1][1], 2)) for x in dc.tracks2D_list]
    ax_length.hist(track_length, bins = 100, range=(0, 60.), histtype="stepfilled", color="#e16543", edgecolor="#4a0d0d") 
    ax_length.set_xlabel('Track 2D length [cm]')

    angle_start = [x.ini_slope for x in dc.tracks2D_list]
    angle_stop  = [x.end_slope for x in dc.tracks2D_list]
    ax_angles.scatter(angle_start, angle_stop, c="#8aa7cf", s=1)
    ax_angles.set_xlabel('initial slope')
    ax_angles.set_ylabel('final slope')

    plt.subplots_adjust(bottom=0.08, top=0.95, hspace=0.3)

    if(option):
        option = "_"+option
    else:
        option = ""


    run_nb = str(dc.evt_list[-1].run_nb)
    evt_nb = str(dc.evt_list[-1].evt_nb_glob)


    plt.savefig('ED/track2D_properties'+option+'_run_'+run_nb+'_evt_'+evt_nb+'.png')
    #plt.show()
    plt.close()



