import matplotlib as mpl



def set_style():    
    params = {'figure.figsize': (12, 6),

              'legend.fontsize':'x-large',              
              'axes.labelsize': 'x-large',
              'axes.titlesize': 'x-large',
              'xtick.labelsize':'x-large',
              'ytick.labelsize':'x-large'}
              
    mpl.rcParams.update(params)
