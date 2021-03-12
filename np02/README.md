# NP02 RUN CONGIFURATION

All runs marked as *.cosmics* are listed in `runs.yaml`.

The number of CRP usefull for the reconstruction is already set. 
What could be missing is the eventual broken/noisy channels. 

Such not yet documented runs have their `daq_broken_channels = [-1]`. 
Hence, by default all channels will be considered in the reconstruction.
If this is really the case, please update it to an empty list: `[]`.

If some channels should be switched off, please add them using **the corresponding DAQ channels**.
The file `convert_channels_to_daq.py` will help you to get the (crp, view, channel) to DAQ channel correspondance.

Don't worry about un-instrumented channels (CRP 2 and some of CRP 3), this is taken care of in lardon.