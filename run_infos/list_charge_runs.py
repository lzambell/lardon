import glob
import pandas as pd
import os



def list_all_runs(folder, last_run):
    d = []
    letter_dict = {'a':1, 'b':2, 'c':3, 'd':4}

    
    for l in folder:

        run_nb = int( l[path_len:])
        if(run_nb <= last_run): continue
        sub_runs = glob.glob(l+"/*")
    
        afile = sub_runs[0]
        datatype = afile[afile.find(".")+1:]
        nb_files = len(sub_runs)
        nevent = nb_files * 30
    
        sub_runs = [ x[path_len + len(str(run_nb)) + 1 + len(str(run_nb)) + 1:] for x in sub_runs]
    
        """there must be a simpler way"""
        sub_runs = [ x[:len(x)-len(datatype)-1] for x in sub_runs]
        sub_list = [( int(x[:x.find("_")]), x[x.find("_")+1:]) for x in sub_runs]
        sub_list.sort()

        first_nb, first_letter = sub_list[0]
        if(len(first_letter)>1):
            print run_nb, " problem ? ", first_letter
            first_letter = first_letter[0]

        last_nb, last_letter = sub_list[-1]
        if(len(last_letter)>1):
            print run_nb, " problem ? ", last_letter
            last_letter = last_letter[0]

        first_sub = str(first_nb)+"_"+first_letter
        last_sub = str(last_nb)+"_"+last_letter
    
    
    
        count_files = (last_nb - first_nb)*4 + letter_dict[last_letter]
        missing = abs(count_files - nb_files)
        if(missing >0):
            print " missing files in run ", run_nb, " nb of files: ", nb_files, " counting ", count_files
        d.append( (run_nb, nb_files, nevent, datatype, first_sub, last_sub, missing) )
    return d




data_path = "/eos/experiment/neutplatform/protodune/rawdata/np02/rawdata"
path_len = len(data_path) + 1

file_all_runs = "all_charge_run_list.csv"
file_cosmics_runs = "cosmics_charge_run_list.csv"
file_pedestal_runs = "pedestal_charge_run_list.csv"


list_all =  glob.glob(data_path+"/*")
list_folder = [x for x in list_all if os.path.isdir(x)]


col = ['run', 
       'n_files',
       'n_event_up',
       'type',
       'first_sub_file',
       'last_sub_file',
       'missing_files']



if(os.path.exists(file_all_runs)):

    current_df = pd.read_csv(file_all_runs)
    n_runs_logged = len(current_df)
    last_run_logged = current_df['run'][n_runs_logged-1]

    print " --> Currently, ", n_runs_logged, " logged up to run ", last_run_logged

    if(n_runs_logged == len(list_folder)):
        print " files are up to date !"
        exit()
    else:
        update = list_all_runs(list_folder, last_run_logged)

        update_df = pd.DataFrame(update, columns=col)        
        update_df.sort(['run'])

        df = [current_df, update_df]
        df = pd.concat([current_df, update_df])
        print " now --> ", len(df), " logged files "

        df.to_csv(file_all_runs,index=False)

        df_cos = df.groupby('type').get_group('cosmics')
        df_cos.to_csv(file_cosmics_runs,index=False)

        df_ped = df.groupby('type').get_group('pedestal')
        df_ped.to_csv(file_pedestal_runs,index=False)
        
else:
    data = list_all_runs(list_folder, 0)

    df = pd.DataFrame(data, columns=col)

    df = df.sort(['run'])
    df.to_csv(file_all_runs,index=False)
    
    df_cos = df.groupby('type').get_group('cosmics')
    df_cos.to_csv(file_cosmics_runs,index=False)
    
    df_ped = df.groupby('type').get_group('pedestal')
    df_ped.to_csv(file_pedestal_runs,index=False)
    
    
