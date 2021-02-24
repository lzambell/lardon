import config as cf
import data_containers as dc

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import collections  as mc
import itertools as itr
import math
import colorcet as cc


adc_min = -5
adc_max = 25
cmap_ed = cc.cm.linear_tritanopic_krjcw_5_95_c24_r


def draw(crp, view, ax=None, t_min = -1, t_max = -1, ch_min = -1, ch_max = -1):

    ax = plt.gca() if ax is None else ax
    
    ax.imshow(dc.data[crp, view, :, :].transpose(), 
              origin = 'lower', 
              aspect = 'auto', 
              #interpolation='none',
              cmap   = cmap_ed,
              vmin   = adc_min, 
              vmax   = adc_max)


    if(ch_min > -1):
        ax.set_xbound(lower=ch_min)
    if(ch_max > -1):
        ax.set_xbound(upper=ch_max)

    if(t_min > -1):
        ax.set_ybound(lower=t_min)
    if(t_max > -1):
        ax.set_ybound(upper=t_max)
    
    return ax


def plot_ed_zoom(crp, view, t_min = -1, t_max = -1, ch_min = -1, ch_max = -1, option=None, to_be_shown=False):

    fig = plt.figure( figsize=(6,6))
    gs = gridspec.GridSpec(nrows=2, ncols=1, height_ratios=[1,30])
    
    ax_data = fig.add_subplot(gs[1, 0])
    ax_data = draw(crp, view, ax=ax_data, t_min = t_min, t_max = t_max, ch_min = ch_min, ch_max = ch_max)
    ax_data.set_xlabel('View Channel')
    ax_data.set_ylabel('Time')
    ax_data.yaxis.set_major_locator(plt.MaxNLocator(4))
    ax_data.set_title('CRP '+str(crp)+' - View '+str(view))

    ax_col = fig.add_subplot(gs[0, 0])
    ax_col.set_title('Collected Charge [ADC]')
        
                         
    cb = fig.colorbar(ax_data.images[-1], cax=ax_col, orientation='horizontal')
    cb.ax.xaxis.set_ticks_position('top')
    cb.ax.xaxis.set_label_position('top')

    plt.tight_layout()

    if(option):
        option = "_"+option
    else:
        option = ""

    run_nb = str(dc.evt_list[-1].run_nb)
    evt_nb = str(dc.evt_list[-1].evt_nb_glob)

    plt.savefig(cf.plot_path+'/ed_zoom_crp_'+str(crp)+'_v'+str(view)+option+'_run_'+run_nb+'_evt_'+evt_nb+'.png')
    if(to_be_shown):
        plt.show()
    plt.close()



def plot_ed(crp_list, nameout="", option=None, to_be_shown=False):
    n_crp = len(crp_list)

    fig = plt.figure(figsize=(9,3*n_crp))
    gs = gridspec.GridSpec(nrows=n_crp+1, 
                           ncols=2, 
                           height_ratios=[1 if i == 0 else 10 for i in range(n_crp+1)])
    ax_col = fig.add_subplot(gs[0,:])
    ax_v0 = [fig.add_subplot(gs[i+1, 0]) for i in range(n_crp) ]
    ax_v1 = [fig.add_subplot(gs[i+1, 1], sharey=ax_v0[i]) for i in range(n_crp) ]
    

    """ draw data """
    i = 0
    for icrp in crp_list:

        """ View 0 """
        ax_v0[i] = draw(icrp, 0, ax=ax_v0[i])
        ax_v0[i].set_title('CRP '+str(icrp)+' - View 0')
        ax_v0[i].set_ylabel('Time')

        """ View 1 """
        ax_v1[i] = draw(icrp, 1, ax=ax_v1[i])        
        ax_v1[i].set_title('CRP '+str(icrp)+' - View 1')        
        ax_v1[i].set_ylabel('Time')#, rotation=270)        
        ax_v1[i].yaxis.tick_right()
        ax_v1[i].yaxis.set_label_position("right")

        i += 1

    """ special care for x-labels """        
    for a,b in zip(ax_v0[:-1], ax_v1[:-1]):
        a.tick_params(labelbottom=False)
        b.tick_params(labelbottom=False)
        
    ax_v0[-1].set_xlabel('View Channel')
    ax_v1[-1].set_xlabel('View Channel')


    for ax in ax_v0:
        ax.yaxis.set_major_locator(plt.MaxNLocator(4))

    ax_col.set_title('Collected Charge [ADC]')
        
                         
    cb = fig.colorbar(ax_v0[0].images[-1], cax=ax_col, orientation='horizontal')
    cb.ax.xaxis.set_ticks_position('top')
    cb.ax.xaxis.set_label_position('top')
        
    plt.tight_layout()
    plt.subplots_adjust(wspace=0.02)

    """
    plt.subplots_adjust(wspace=0.02, 
                        hspace=0.25, 
                        top=0.92, 
                        bottom=0.08, 
                        left=0.1, 
                        right=0.9)
    """

    if(option):
        option = "_"+option
    else:
        option = ""

    run_nb = str(dc.evt_list[-1].run_nb)
    evt_nb = str(dc.evt_list[-1].evt_nb_glob)

    plt.savefig(cf.plot_path+'/ed_'+nameout+option+'_run_'+run_nb+'_evt_'+evt_nb+'.png')
    if(to_be_shown):
        plt.show()
    plt.close()


def plot_ed_data(option=None, to_be_shown=False):
    if(cf.n_CRPUsed==4):
        plot_ed(crp_list=[0,1,3], nameout="data", option=option, to_be_shown=to_be_shown)
    else:
        plot_ed_two_crps(crp1=0, crp2=1, option=option, to_be_shown=to_be_shown)

def plot_ed_one_crp(crp, option=None, to_be_shown=False):
    plot_ed(crp_list=[crp], nameout="crp_"+str(crp), option=option, to_be_shown=to_be_shown)

def plot_ed_two_crps(crp1=0, crp2=1, option=None, to_be_shown=False):
    plot_ed(crp_list=[crp1, crp2], nameout="crp_"+str(crp1)+"_"+str(crp2), option=option, to_be_shown=to_be_shown)

def plot_ed_three_crps(crp1=0, crp2=1, crp3=3, option=None, to_be_shown=False):
    plot_ed(crp_list=[crp1, crp2, crp3], nameout="crp_"+str(crp1)+"_"+str(crp2)+"_"+str(crp3), option=option, to_be_shown=to_be_shown)

def plot_ed_all(option=None, to_be_shown=False):
    plot_ed(crp_list=[0, 1, 2, 3], nameout="allcrp", option=option, to_be_shown=to_be_shown)
