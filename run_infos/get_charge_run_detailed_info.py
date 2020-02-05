import os
import sys
import glob
import pandas as pd
import numpy as np

if(len(sys.argv) != 4):
    exit()


run_nb = sys.argv[1]
sub_file_ini = sys.argv[2]
sub_file_fin = sys.argv[3]

print " looking at run ", run_nb, " files ", sub_file_ini, " -> ", sub_file_fin
data_path = "/eos/experiment/neutplatform/protodune/rawdata/np02/rawdata"

subnb_ini  = int(sub_file_ini[:sub_file_ini.find("_")])
sublet_ini = sub_file_ini[sub_file_ini.find("_")+1:]

subnb_fin  = int(sub_file_fin[:sub_file_fin.find("_")])
sublet_fin = sub_file_fin[sub_file_fin.find("_")+1:]

subfiles = []
for isub in range(subnb_ini,subnb_fin+1):
    for ichunk in ['a','b','c','d']:
        subfiles.append(str(isub)+"_"+ichunk)




infos = []
event_times = []
ntot = 0
earliest_time = -1.
latest_time   = 0.

for idx, sub in reversed(list(enumerate(subfiles))): 
    if(sub == sub_file_fin):
        subfiles = subfiles[:idx+1]
        break
        
for ifile in subfiles:
    file_path = data_path+"/"+run_nb+"/"+run_nb+"_"+ifile+".cosmics"
    if(os.path.exists(file_path) is False):
        print " --> file ", ifile, " does not exist !"
        #Add that info on file
        infos.append((ifile, -1, -1, -1, -1, -1))
        break

    
    with open(file_path, "rb") as data:
        run, nb_evt = np.fromfile(data, dtype='<u4', count=2)
        ntot += nb_evt

        sequence = []
        print " run :", run, " file ", ifile, " : nb of event ", nb_evt
        
        n_flag_ok_v0 = 0
        n_flag_ok_v1 = 0
        tfirst = -1.
        tlast = 0.

        for i in range(nb_evt):
            seq  = np.fromfile( data, dtype='<u4', count=4)
            # 4 uint of [event number - event total size with header- event data size - 0]
            sequence.append(seq[1])
        event_pos = []
        event_pos.append( data.tell() )
        for i in range(nb_evt-1):
            data.seek(sequence[i], 1)
            event_pos.append( data.tell() ) #get the byte position of each event
        
        for idx in range(nb_evt):
            data.seek(event_pos[idx], 0)
            data.seek(15, 1)
            time_s, time_ns = np.fromfile(data, dtype='<i8', count=2)
            flag_v0 = np.fromfile(data, dtype='<B', count=1)
            evt_num, lro, cro = np.fromfile(data, dtype='<i4', count=3)
            data.seek(lro+cro+1 + 31, 1)
            flag_v1 = np.fromfile(data, dtype='<B', count=1)        

            if(flag_v0 == 25): n_flag_ok_v0 +=1
            if(flag_v1 == 25): n_flag_ok_v1 +=1
            time = np.float64(time_s) + np.float64(time_ns*1e-9)
            if(time < tfirst or tfirst < 0):
                tfirst = time
            if(time > tlast):
                tlast = time
            event_times.append((int(run), ifile, idx, int(evt_num), "%.4f"%time))
        if(earliest_time < 1 or tfirst < earliest_time):
            earliest_time = tfirst
        if(tlast > latest_time):
            latest_time = tlast
        infos.append((ifile, nb_evt, "%.3f"%tfirst, "%.3f"%tlast, n_flag_ok_v0, n_flag_ok_v1))

col_infos = ['sub_run', 
             'n_events',
             'first_evt_time',
             'last_evt_time',
             'nb_flag_ok_v0',
             'nb_flag_ok_v1',]

df_gen = pd.DataFrame(infos, columns=col_infos)
df_gen.to_csv("run_"+run_nb+"_infos.csv",index=False)

col_times = ['run', 'file', 'file_event', 'run_event','time']
df_times = pd.DataFrame(event_times, columns=col_times)
df_times = df_times.sort(['run_event'])
#df_times.sort_values(by='run_event')
df_times.to_csv("run_"+run_nb+"_times.csv",index=False)
print type(ntot)

len_s = (latest_time - earliest_time)
rate =  ntot/len_s
print "in total ", ntot, " event taken in ", latest_time - earliest_time, "s <-> r = ", rate, " Hz"
print type(ntot)

df_cosmics = pd.read_csv("cosmics_charge_run_list.csv")
df_cosmics.set_index('run', inplace=True)
print df_cosmics.dtypes
df_cosmics.at[int(run_nb),'nb_event'] = ntot
df_cosmics.at[int(run_nb),'tstart'] = "%.3f"%earliest_time
df_cosmics.at[int(run_nb),'tstop'] = "%.3f"%latest_time
df_cosmics.at[int(run_nb),'len_s'] = "%.3f"%len_s
df_cosmics.at[int(run_nb),'rate'] = "%.1f"%rate


df_cosmics.to_csv("cosmics_charge_run_list.csv")

