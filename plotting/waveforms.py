import config as cf
import data_containers as dc

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import collections  as mc
from matplotlib.legend_handler import HandlerTuple
import itertools as itr
import math
import colorcet as cc




def draw_current_waveform(crp, view, ch, ax=None, **kwargs):    
    
    ax = plt.gca() if ax is None else ax
    
    ax.plot(dc.data[crp, view, ch, :], **kwargs)
    return ax


def plot_wvf_single_current(wvf_list, adc_min=-1, adc_max=-1, option=None, to_be_shown=False):

    """ wvf_list should be a list of tuples as (crp, view, ch)"""

    n_wvf = len(wvf_list)

    fig = plt.figure(figsize=(12, 3*n_wvf))
    gs = gridspec.GridSpec(nrows=n_wvf, ncols=1)
    ax = [fig.add_subplot(gs[i,0]) for i in range(n_wvf)]

    
    for i in range(n_wvf):
        (crp, view, ch) = wvf_list[i]
        legend = "crp "+str(crp)+ " v"+str(view)+" ch"+str(ch)
        ax[i] = draw_current_waveform(crp, view, ch, ax=ax[i], label=legend, c='k')
        ax[i].set_ylabel('ADC')
        ax[i].legend(loc='upper right')
        ax[i].set_xlim([0, 10000])
        if(adc_min > -1):
            ax[i].set_ybound(lower=adc_min)
        if(adc_max > -1):
            ax[i].set_ybound(upper=adc_max)
            
    for a in ax[:-1]:
        a.tick_params(labelbottom=False)
    ax[-1].set_xlabel('Time')
    

    plt.tight_layout()

    if(option):
        option = "_"+option
    else:
        option = ""


    run_nb = str(dc.evt_list[-1].run_nb)
    evt_nb = str(dc.evt_list[-1].evt_nb_glob)

    plt.savefig(cf.plot_path+'/waveforms'+option+'_run_'+run_nb+'_evt_'+evt_nb+'.png')
    if(to_be_shown):
        plt.show()
    plt.close()


def plot_wvf_multi_current(wvf_list, adc_min=-1, adc_max=-1, option=None, to_be_shown=False):
    """ wvf_list should be a list of tuples as (crp, view, ch)"""

    n_wvf = len(wvf_list)

    fig = plt.figure(figsize=(12, 3))
    ax = plt.gca()

    
    for i in range(n_wvf):
        (crp, view, ch) = wvf_list[i]
        legend = "crp "+str(crp)+ " v"+str(view)+" ch"+str(ch)
        ax = draw_current_waveform(crp, view, ch, ax=ax, label=legend)
        
    ax.set_ylabel('ADC')
    ax.legend(loc='upper right')
    ax.set_xlim([0, 10000])
    if(adc_min > -1):
        ax.set_ybound(lower=adc_min)
    if(adc_max > -1):
        ax.set_ybound(upper=adc_max)
    
    ax.set_xlabel('Time')
    

    plt.tight_layout()

    if(option):
        option = "_"+option
    else:
        option = ""


    run_nb = str(dc.evt_list[-1].run_nb)
    evt_nb = str(dc.evt_list[-1].evt_nb_glob)

    plt.savefig(cf.plot_path+'/waveforms_multi'+option+'_run_'+run_nb+'_evt_'+evt_nb+'.png')
    if(to_be_shown):
        plt.show()
    plt.close()



def draw_data(data, ax=None, **kwargs):
    
    ax = plt.gca() if ax is None else ax    
    ax.plot(data, **kwargs)
    return ax
    

def plot_wvf_evo(data, title="", legends=[], adc_min=-1, adc_max=-1, option=None, to_be_shown=False):

    fig = plt.figure(figsize=(12, 3))
    ax = plt.gca()

    for d in range(len(data)):
        ax = draw_data(data[d], ax=ax, label=legends[d])
        
    ax.set_xlabel('Time')
    ax.set_ylabel('ADC')
    ax.legend(loc='upper right')
    ax.set_xlim([0, 10000])
    if(adc_min > -1):
        ax.set_ybound(lower=adc_min)
    if(adc_max > -1):
        ax.set_ybound(upper=adc_max)
    
    ax.set_title(title)
    

    plt.tight_layout()

    if(option):
        option = "_"+option
    else:
        option = ""


    run_nb = str(dc.evt_list[-1].run_nb)
    evt_nb = str(dc.evt_list[-1].evt_nb_glob)

    plt.savefig(cf.plot_path+'/waveform_evo'+option+'_run_'+run_nb+'_evt_'+evt_nb+'.png')
    if(to_be_shown):
        plt.show()
    plt.close()
