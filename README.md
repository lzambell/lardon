#Liquid Argon Reconstruction Done in PythON
![Logo](lardon.png)


## How to run lardon
You need miniconda installed :
https://docs.conda.io/en/latest/miniconda.html#linux-installers

and then get the librairies as stated in lardenv.yml : 
conda env create -f lardenv.yml
[It'll take about 2-3 GB!]

then :
conda activate lardenv

Before running anything, please check config.py :
 - modify the store_path and plot_path to provide repositories on where to store the output and control plots (if any)
 - comment / uncomment run specific parameters (Nb of CRP taking data, Broken Channels, Drift Field, LAr Temperature, ...) -> This part will be automatized at some point

To run lardon on data : 
python reader.py -run <run nb> -sub <subfile name> -n <nb of events (optional)> -out <output file option>
e.g. : python reader.py -run 1415 -sub 1_a -n 10 -out example


To run lardon on MC : 
add -mc <the MC root file to be read>
the simulation will be added to a run of data (preferably a noise run!)
Note that the MC file should be in a specific format - This part is still under development

** CRP/VIEW CONVENTION **
(y)
 1  
 w  |-------------> View 0 (x) 
 e  -------------------
 i  |       |         |
 V  |   1   |    0    |
 ^  |       |         |
 |  --------•----------
 |  |       |         |
 |  |   2   |    3    |
 |  |       |         |
 -  -------------------
 
 origin in the center (at •)
 electrons drift along z axis
 all distance are in cm


By default, no control plots are produced, but you can call the plotting functions in reader.py anywhere in the reconstruction loop
   all plot function have the two options 
   option="extra_output_name_if_you_want" [default is none] 
   to_be_shown=True/False if you want to see the plot live [default is False]

 - The current event display with :
     plot.plot_ed_data()
     take a look at plotting/event_display.py to see all possible plots (just one crp/one view, zoom, ...)
  - Current waveform(s) with :
    plot.plot_wvf_single_current([(crp,view,ch),(crp,view,ch),..]) [one wvf/plot]
    plot.plot_wvf_multi_current([(crp,view,ch),(crp,view,ch),..]) [all wvf superimposed]
 - waveform evolution with : 
    plot.plot_wvf_evo([wvf1, wvf2,...], title="your plot title", legends=['leg1', 'leg2',...]
    where wvf1 would be a copy of the waveform at some point in the reconstruction with :
    wvf1 = dc.data(crp, view, ch, :].copy()
 - hit found / cluster found : 
   plot.plot_2dcrp_clusters() <- shown in CRP/View format
   plot.plot_2dview_clusters() <- shown in x/y format
   [same for plot_2dxxx_hits()]
 - 2D tracks (and hits)
   plot.plot_2dview_2dtracks()
 -3D tracks
   plot.plot_2dview_hits_and_3dtracks() <- two views separated (all hits and 2d tracks also shown)
   plot.plot_3d() <- in 3d 

