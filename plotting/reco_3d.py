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
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.cbook import flatten

color_v0 = 'c'
color_v1 = 'orange'


def plot_3d(option=None, to_be_shown=True):
    
    fig = plt.figure(figsize=(6, 6))
    ax  = fig.add_subplot(111, projection='3d')

    x0, y0, z0 = list(flatten(get_3dtracks_x(0))), list(flatten(get_3dtracks_y(0))), list(flatten(get_3dtracks_z(0)))

    x1, y1, z1 = list(flatten(get_3dtracks_x(1))), list(flatten(get_3dtracks_y(1))), list(flatten(get_3dtracks_z(1)))

    """ shadow on the walls """
    ax.scatter(x0,
               [-100 for i in y0], 
               z0, 
               c="#f2f2f2", s=4)


    ax.scatter([300 for i in x1], 
               y1,
               z1,
               c="#f2f2f2", s=4)

    
    ax.scatter(x0, 
               y0, 
               z0, 
               c=color_v0, s=4)

    ax.scatter(x1, y1, z1,
               c=color_v1, s=4)




    ax.set_xlim3d(-300, 300)
    ax.set_ylim3d(-100, 300)
    ax.set_zlim3d(-30, 300)
    
    ax.set_xlabel('View 0/X [cm]')
    ax.set_ylabel('View 1/Y [cm]')
    ax.set_zlabel('Drift/Z [cm]')


    ax.grid(False)
    ax.xaxis.pane.set_edgecolor('black')
    ax.yaxis.pane.set_edgecolor('black')
    ax.zaxis.pane.set_edgecolor('black')
    ax.xaxis.pane.fill = False
    ax.yaxis.pane.fill = False
    ax.zaxis.pane.fill = False
    ax.view_init(elev=10, azim=135)


    ax.xaxis.set_major_locator(plt.MaxNLocator(5))
    ax.yaxis.set_major_locator(plt.MaxNLocator(5))
    ax.zaxis.set_major_locator(plt.MaxNLocator(7))
    

    ax.xaxis._axinfo['tick']['inward_factor'] = 0.4
    ax.xaxis._axinfo['tick']['outward_factor'] = 0.4
    ax.yaxis._axinfo['tick']['inward_factor'] = 0.4
    ax.yaxis._axinfo['tick']['outward_factor'] = 0.4
    ax.zaxis._axinfo['tick']['outward_factor'] = 0.3
    ax.zaxis._axinfo['tick']['inward_factor'] = 0.3
    

    plt.subplots_adjust(top=0.95,
                        bottom=0.05,
                        left=0.05,
                        right=0.95,
                        hspace=0.2,
                        wspace=0.2)


    if(option):
        option = "_"+option
    else:
        option = ""


    run_nb = str(dc.evt_list[-1].run_nb)
    evt_nb = str(dc.evt_list[-1].evt_nb_glob)


    plt.savefig(cf.plot_path+'/track3D'+option+'_run_'+run_nb+'_evt_'+evt_nb+'.png')
    if(to_be_shown):
        plt.show()
    plt.close()

