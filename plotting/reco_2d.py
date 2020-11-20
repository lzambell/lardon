import config as cf
import data_containers as dc

from .select_hits import *
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import collections  as mc
from matplotlib.legend_handler import HandlerTuple
import itertools as itr
import math
import colorcet as cc


adc_min = -5
adc_max = 20
cmap_ed = cc.cm.kbc_r
marker_size = 5
"""define default color cycle from colorcet"""
mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=cc.glasbey_warm)

color_noise     = "#c7c7c7"
color_clustered = "#ffb77d"
color_matched1  = "#28568f"
color_matched2  = "#abdf7f"
color_track2d   = "#de425b"
color_track3d   = "#00a9b2"

def draw_hits(pos, time, z=[], ax=None, **kwargs):

    ax = plt.gca() if ax is None else ax

    if(len(pos) == 0):
        return ax
    
    if(len(z) > 0):
        ax.scatter(pos, time, c=z, **kwargs)
    else:
        ax.scatter(pos, time, **kwargs)
    return ax

def draw_all_hits(ax_v0, ax_v1, sel='True', adc=False, charge=False, **kwargs):

    axs = [ax_v0, ax_v1]

    for iview in range(2):
        z = []
        if(adc==True):
            z = get_hits_adc(iview, sel)
        elif(charge==True):
            z = get_hits_charge(iview, sel)

        axs[iview] = draw_hits(pos=get_hits_pos(iview, sel), 
                               time=get_hits_z(iview, sel), 
                               z=z,
                               ax=axs[iview], **kwargs)
    return axs[0], axs[1]



def draw_tracks(pos, time, ax=None, legend="", **kwargs):

    ax = plt.gca() if ax is None else ax
    
    if(len(pos) == 0):
        return ax

    if(len(legend)>0):
        ax.plot(pos[0], time[0], label=legend, **kwargs)
        
    for tx,tz in zip(pos, time):
        ax.plot(tx,tz, **kwargs)
    
    return ax


def draw_all_tracks(ax_v0, ax_v1, sel='True', legend="", **kwargs):

    axs = [ax_v0, ax_v1]

    for iview in range(2):
        axs[iview] = draw_tracks(pos=get_2dtracks_pos(iview,sel), 
                                 time=get_2dtracks_z(iview,sel), 
                                 ax=axs[iview], 
                                 legend=legend,
                                 **kwargs)
    return axs[0], axs[1]



def template_data_view():

    fig = plt.figure(figsize=(12,6))
    gs = gridspec.GridSpec(nrows=2, ncols=2, 
                           height_ratios=[1,20], 
                           width_ratios=[6,4])
    
    ax_col = fig.add_subplot(gs[0,:])
    ax_v0 = fig.add_subplot(gs[1, 0])
    ax_v1 = fig.add_subplot(gs[1, 1], sharey=ax_v0)

    ax_v0.set_title('View 0')
    ax_v0.set_ylabel('Z [cm]')
    ax_v0.set_xlabel('X [cm]')
    ax_v0.set_xlim([-300., 300])
    ax_v0.set_ylim([-30., 300])

    ax_v1.set_title('View 1')
    ax_v1.set_ylabel('Z [cm]')
    ax_v1.yaxis.tick_right()
    ax_v1.yaxis.set_label_position("right")
    ax_v1.set_xlabel('Y [cm]')
    ax_v1.set_xlim([-100., 300])
    ax_v1.set_ylim([-30., 300])

    plt.subplots_adjust(top=0.9,
                        bottom=0.11,
                        left=0.08,
                        right=0.92,
                        hspace=0.3,
                        wspace=0.1)
    
    return fig, ax_col, ax_v0, ax_v1


def template_data_crp():
    fig = plt.figure(figsize=(9,9))    
    gs = gridspec.GridSpec(nrows=4, 
                           ncols=2, 
                           height_ratios=[1, 10, 10, 10])
    ax_col = fig.add_subplot(gs[0,:])
    ax_v0 = [fig.add_subplot(gs[i+1, 0]) for i in range(3) ]
    ax_v1 = [fig.add_subplot(gs[i+1, 1], sharey=ax_v0[i]) for i in range(3) ]


    """ CRP/view min/max positions (as in [crp,view]) """
    #pos_min = np.array([[0, 0],[-300,0],[-300,-300],[0,-300]])
    #pos_max = np.array([[300, 300],[0,300],[0,0],[300,0]])


    for icrp in range(3):
        the_crp = icrp+1 if icrp==2 else icrp
        ax_v0[icrp].set_title('CRP '+str(the_crp)+' - View 0')
        ax_v0[icrp].set_ylabel('Time')
        ax_v0[icrp].set_xlim([0, 960])#[pos_min[the_crp,0], pos_max[the_crp,0]])
        ax_v0[icrp].set_ylim([0, 10000])#-30., 300])
        
        ax_v1[icrp].set_title('CRP '+str(the_crp)+' - View 1')
        ax_v1[icrp].set_ylabel('Time')
        ax_v1[icrp].yaxis.tick_right()
        ax_v1[icrp].yaxis.set_label_position("right")
        #ax_v1[icrp].set_xlabel('Y [cm]')
        #ax_v1[icrp].set_xlim([pos_min[the_crp,1], pos_max[the_crp,1]])
        #ax_v1[icrp].set_ylim([-30., 300])
        ax_v1[icrp].set_xlim([0, 960])
        ax_v1[icrp].set_ylim([0, 10000])

    """ special care for x-labels """        
    for a,b in zip(ax_v0[:-1], ax_v1[:-1]):
        a.tick_params(labelbottom=False)
        b.tick_params(labelbottom=False)
        
    ax_v0[-1].set_xlabel('View Channel')
    ax_v1[-1].set_xlabel('View Channel')

    """
    for ax in ax_v0:
        ax.yaxis.set_major_locator(plt.MaxNLocator(4))
    """

    """
    plt.subplots_adjust(top=0.903,
                        bottom=0.098,
                        left=0.097,
                        right=0.903,
                        hspace=1.0,
                        wspace=0.186)
    """
    
    return fig, ax_col, ax_v0, ax_v1



def plot_2dcrp_hits(option=None, to_be_shown=False):

    fig, ax_col, ax_v0, ax_v1 = template_data_crp()

    max_adc=50
    i=0
    for icrp in range(4):
        if(icrp==2): continue
        sel = 'x.crp=='+str(icrp)
        ax_v0[i] = draw_hits(pos=get_hits_ch(0, sel), time=get_hits_tdc(0, sel), z= get_hits_adc(0, sel), ax=ax_v0[i], cmap=cmap_ed, s=marker_size, vmin=0, vmax=max_adc)
        ax_v1[i] = draw_hits(pos=get_hits_ch(1, sel), time=get_hits_tdc(1, sel), z= get_hits_adc(1, sel), ax=ax_v1[i],  cmap=cmap_ed, s=marker_size, vmin=0, vmax=max_adc)
        i+=1


    """ doesn't work for some reasons"""
    for ax0,ax1 in zip(ax_v0, ax_v1):
        ax0.yaxis.set_major_locator(plt.MaxNLocator(4))
        ax1.yaxis.set_major_locator(plt.MaxNLocator(4))
        

    """ color bar """
    ax_col.set_title('Hit Max ADC')

    cb = fig.colorbar(ax_v0[0].collections[0], cax=ax_col, orientation='horizontal')
    cb.ax.xaxis.set_ticks_position('top')
    cb.ax.xaxis.set_label_position('top')

    plt.tight_layout()
    plt.subplots_adjust(wspace=0.02)


    if(option):
        option = "_"+option
    else:
        option = ""


    run_nb = str(dc.evt_list[-1].run_nb)
    evt_nb = str(dc.evt_list[-1].evt_nb_glob)


    plt.savefig('ED/hit_crp'+option+'_run_'+run_nb+'_evt_'+evt_nb+'.png')

    if(to_be_shown):
        plt.show()
    plt.close()





def plot_2dview_hits(option=None, to_be_shown=False):

    fig, ax_col, ax_v0, ax_v1 = template_data_view()

    max_adc=50
    ax_v0, ax_v1 = draw_all_hits(ax_v0, ax_v1, adc=True, cmap=cmap_ed, s=marker_size, vmin=0, vmax=max_adc)
    
    """ color bar """
    ax_col.set_title('Hit Max ADC')

    cb = fig.colorbar(ax_v0.collections[0], cax=ax_col, orientation='horizontal')
    cb.ax.xaxis.set_ticks_position('top')
    cb.ax.xaxis.set_label_position('top')


    if(option):
        option = "_"+option
    else:
        option = ""


    run_nb = str(dc.evt_list[-1].run_nb)
    evt_nb = str(dc.evt_list[-1].evt_nb_glob)


    plt.savefig('ED/hit_view'+option+'_run_'+run_nb+'_evt_'+evt_nb+'.png')
    if(to_be_shown):
        plt.show()
    plt.close()



def plot_2dcrp_clusters(option=None, to_be_shown=False):
    fig, ax_col, ax_v0, ax_v1 = template_data_crp()

    """ clustered hits """
    i=0
    for icrp in range(cf.n_CRPUsed):
        if(icrp == 2): continue
        for icl in range(dc.evt_list[-1].nClusters[icrp, 0]):
            sel = 'x.crp=='+str(icrp)+' and x.cluster=='+str(icl)
            ax_v0[i] = draw_hits(pos=get_hits_ch(0, sel), time=get_hits_tdc(0, sel), ax=ax_v0[i], s=marker_size, marker='o')

        for icl in range(dc.evt_list[-1].nClusters[icrp, 1]):
            sel = 'x.crp=='+str(icrp)+' and x.cluster=='+str(icl)        
            ax_v1[i] = draw_hits(pos=get_hits_ch(1, sel), time=get_hits_tdc(1, sel), ax=ax_v1[i], s=marker_size, marker='o')
        

        """ unclustered hits """
        sel = 'x.crp=='+str(icrp)+' and x.cluster==-1'
        ax_v0[i] = draw_hits(pos=get_hits_ch(0, sel), time=get_hits_tdc(0, sel), ax=ax_v0[i], c=color_noise, s=marker_size, marker='o')
        ax_v1[i] = draw_hits(pos=get_hits_ch(1, sel), time=get_hits_tdc(1, sel), ax=ax_v1[i], c=color_noise, s=marker_size, marker='o')

        i+=1


    """ for some reasons, this doesn't work """
    for ax0,ax1 in zip(ax_v0, ax_v1):
        ax0.yaxis.set_major_locator(plt.MaxNLocator(4))
        ax1.yaxis.set_major_locator(plt.MaxNLocator(4))


    ax_col.axis('off')


    plt.tight_layout()
    plt.subplots_adjust(wspace=0.02)

    if(option):
        option = "_"+option
    else:
        option = ""


    run_nb = str(dc.evt_list[-1].run_nb)
    evt_nb = str(dc.evt_list[-1].evt_nb_glob)


    plt.savefig('ED/cluster_2dcrp'+option+'_run_'+run_nb+'_evt_'+evt_nb+'.png')
    if(to_be_shown):
        plt.show()
    plt.close()


def plot_2dview_clusters(option=None, to_be_shown=False):
    
    fig, ax_col, ax_v0, ax_v1 = template_data_view()

    """ clustered hits """
    for icrp in range(cf.n_CRPUsed):
        if(icrp == 2): continue
        for icl in range(dc.evt_list[-1].nClusters[icrp, 0]):
            sel = 'x.crp=='+str(icrp)+' and x.cluster=='+str(icl)
            ax_v0 = draw_hits(pos=get_hits_pos(0, sel), time=get_hits_z(0, sel), ax=ax_v0, s=marker_size, marker='o')

        for icl in range(dc.evt_list[-1].nClusters[icrp, 1]):
            sel = 'x.crp=='+str(icrp)+' and x.cluster=='+str(icl)        
            ax_v1 = draw_hits(pos=get_hits_pos(1, sel), time=get_hits_z(1, sel), ax=ax_v1, s=marker_size, marker='o')



    """ unclustered hits """
    sel = 'x.cluster==-1'
    ax_v0, ax_v1 = draw_all_hits(ax_v0, ax_v1, sel, c=color_noise, s=marker_size, marker='o')

    ax_col.axis('off')

    if(option):
        option = "_"+option
    else:
        option = ""


    run_nb = str(dc.evt_list[-1].run_nb)
    evt_nb = str(dc.evt_list[-1].evt_nb_glob)


    plt.savefig('ED/cluster_2dview'+option+'_run_'+run_nb+'_evt_'+evt_nb+'.png')
    if(to_be_shown):
        plt.show()
    plt.close()






def plot_2dview_2dtracks(option=None, to_be_shown=False):
    fig, ax_leg, ax_v0, ax_v1 = template_data_view()
    
    """ all hits """
    ax_v0, ax_v1 = draw_all_hits(ax_v0, ax_v1, c="#e6e6e6", s=marker_size, marker='o', label='Noise Hits')

    
    
    """ 2D tracks """
    ax_v0, ax_v1 = draw_all_tracks(ax_v0, ax_v1, linewidth=1, legend='2D Track')


    """ legend """
    ax_leg.axis('off')

    
    if(option):
        option = "_"+option
    else:
        option = ""


    run_nb = str(dc.evt_list[-1].run_nb)
    evt_nb = str(dc.evt_list[-1].evt_nb_glob)

    plt.savefig('ED/alltrack2D'+option+'_run_'+run_nb+'_evt_'+evt_nb+'.png')
    if(to_be_shown):
        plt.show()
    plt.close()



def plot_2dview_hits_2dtracks(option=None, to_be_shown=False):
    fig, ax_leg, ax_v0, ax_v1 = template_data_view()
    
    """ unclustered hits """
    sel = 'x.cluster == -1'    
    ax_v0, ax_v1 = draw_all_hits(ax_v0, ax_v1, sel, c=color_noise, s=marker_size, marker='o', label='Noise Hits')

    """ clustered hits """
    sel = 'x.cluster > -1'    
    ax_v0, ax_v1 = draw_all_hits(ax_v0, ax_v1, sel, c=color_clustered, s=marker_size, marker='o', label='Hits Clustered')

    """ delta_ray hits attached to track """
    sel = 'x.matched <0 and x.matched > -9999'
    ax_v0, ax_v1 = draw_all_hits(ax_v0, ax_v1, sel, c=color_matched2, s=marker_size, marker='o', label='Delta Rays')

    """ hits attached to track """
    sel = 'x.matched >= 0'
    ax_v0, ax_v1 = draw_all_hits(ax_v0, ax_v1, sel, c=color_matched1, s=marker_size, marker='o', label='Hits Attached to Track')

    
    """ 2D tracks """
    ax_v0, ax_v1 = draw_all_tracks(ax_v0, ax_v1, legend='2D Track', c=color_track2d, linewidth=1)


    """ legend """
    ax_leg.axis('off')
    
    """ re-arrange the legend (line last), and merge blue and green entries """
    """ might not work anymore if the plotting order is changed """
    h, l = ax_v0.get_legend_handles_labels()

    if(False): #len(h)==5):
        leg = ax_leg.legend([h[1], h[2], (h[3], h[4]), h[0]], [l[1], l[2], 'Hits Attached to Track (1,2)', l[0]], loc='center', ncol=4, markerscale=4, handler_map={tuple: HandlerTuple(ndivide=None)})
    else:
        """otherwise this works well """
        leg = ax_leg.legend(*ax_v0.get_legend_handles_labels(),loc='center', ncol=5, markerscale=4, markerfirst=True)
    
    """ make line in the legend bigger """
    for line in leg.get_lines():
        line.set_linewidth(3)
    
    if(option):
        option = "_"+option
    else:
        option = ""


    run_nb = str(dc.evt_list[-1].run_nb)
    evt_nb = str(dc.evt_list[-1].evt_nb_glob)

    plt.savefig('ED/track2D'+option+'_run_'+run_nb+'_evt_'+evt_nb+'.png')
    if(to_be_shown):
        plt.show()
    plt.close()



def plot_2dview_hits_and_3dtracks(option=None, to_be_shown=False):

    fig, ax_leg, ax_v0, ax_v1 = template_data_view()
    
    """ unclustered hits """
    sel = 'x.cluster == -1'    
    ax_v0, ax_v1 = draw_all_hits(ax_v0, ax_v1, sel, c=color_noise, s=marker_size, marker='o', label='Noise')

    """ clustered hits """
    sel = 'x.cluster > -1'    
    ax_v0, ax_v1 = draw_all_hits(ax_v0, ax_v1, sel, c=color_clustered, s=marker_size, marker='o', label='Clustered')

    """ delta rays attached to track """
    sel = 'x.matched <0 and x.matched > -9999'
    ax_v0, ax_v1 = draw_all_hits(ax_v0, ax_v1, sel, c=color_matched2, s=marker_size, marker='o', label='Delta Rays')

    """ hits attached to track """
    sel = 'x.matched >= 0'
    ax_v0, ax_v1 = draw_all_hits(ax_v0, ax_v1, sel, c=color_matched1, s=marker_size, marker='o', label='Attached to Track')

    
    """ 2D tracks """
    ax_v0, ax_v1 = draw_all_tracks(ax_v0, ax_v1, legend='2D Track', c=color_track2d, linewidth=1)


    """ 3D tracks """
    sel = 't.matched >= 0'
    ax_v0, ax_v1 = draw_all_tracks(ax_v0, ax_v1, sel, c=color_track3d, linewidth=2, legend='3D Track')

    
    """ legend """
    ax_leg.axis('off')

    """ re-arrange the legend (lines last), and merge blue and green entries """
    h, l = ax_v0.get_legend_handles_labels()

    if(False): #len(h)==6):
        leg = ax_leg.legend([h[2], h[3], (h[4], h[5]), h[0], h[1]], [l[2], l[3], 'Hits Attached to Track (1,2)', l[0], l[1]], loc='center', ncol=5, markerscale=4, handler_map={tuple: HandlerTuple(ndivide=None)})
    else:
        leg = ax_leg.legend(*ax_v0.get_legend_handles_labels(),loc='center', ncol=6, markerscale=4, markerfirst=True)

    """ make lines in the legend bigger """
    for line in leg.get_lines():
        line.set_linewidth(3)


    if(option):
        option = "_"+option
    else:
        option = ""


    run_nb = str(dc.evt_list[-1].run_nb)
    evt_nb = str(dc.evt_list[-1].evt_nb_glob)

    plt.savefig('ED/track3D_proj'+option+'_run_'+run_nb+'_evt_'+evt_nb+'.png')
    if(to_be_shown):
        plt.show()
    plt.close()



