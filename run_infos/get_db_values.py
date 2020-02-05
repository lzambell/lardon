import numpy as np
import urllib, json
from datetime import datetime
import slow_control_dict as sc
import matplotlib.pyplot as plt
import sys


def merge_dicts(dicts):#, y):
    """Given two dicts, merge them into a new dict as a shallow copy."""
    z = dicts[0].copy()
    for idi in range(len(dicts)-1):
        z.update(dicts[idi+1])
    return z

slow_control = "https://np02-slow-control.web.cern.ch/np02-slow-control/app/php-db-conn/histogramrange.conn.php?"

if(len(sys.argv) != 5):
    print " needs <run nb> <nb events> <timestamp start> <timestamp stop>"
    exit()

run_nb = int(sys.argv[1]) #1267 #1323
nb_event = int(sys.argv[2]) #4505 #11464
time_start = float(sys.argv[3]) #1570117084.08 #1574347905.61
time_stop = float(sys.argv[4]) #1570117535.33 #1574349053.5


output = open("run_"+str(run_nb)+"_db.py","w")
output.write("run_parameters={\n")
output.write("\t'run':%i,\n"%run_nb)
output.write("\t'nb_of_events':%i,\n"%nb_event)

#dates are DD-MM-YYYY

dt_start   = datetime.fromtimestamp(time_start)
date_start = dt_start.strftime("%d-%m-%Y")


"""request needs two different days to work"""
if(time_stop - time_start < 86400):
    """artificially adds 24h"""
    dt_stop   = datetime.fromtimestamp(time_stop + 86400.)
    date_stop = dt_stop.strftime("%d-%m-%Y")
else:
    dt_stop   = datetime.fromtimestamp(time_stop)
    date_stop = dt_stop.strftime("%d-%m-%Y")

output.write("\t'time_start':'%s',\n"%dt_start.strftime("%d-%m-%Y %H:%M:%S"))
dt_real_stop = datetime.fromtimestamp(time_stop)
output.write("\t'time_stop':'%s',\n"%dt_real_stop.strftime("%d-%m-%Y %H:%M:%S"))


print " from ", date_start, " to ", date_stop
print " --> length of run ", time_stop - time_start, " s or %.2f min"%((time_stop - time_start)/60.)


tspark_start = []
tspark_max   = []
tspark_end   = []
spark_id = []

all_dict = merge_dicts([sc.voltage, sc.level_meter, sc.cryo, sc.ambiant])

for name, elemId in all_dict.items(): 
    print " --- > > > ", name
    
    url = slow_control+"elemId="+str(elemId)+"&start="+date_start+"&end="+date_stop


    response = urllib.urlopen(url)
    response_read = response.read()
    print type(response), " ", len(response_read)
    if(len(response_read) < 20):
        print "no json measurement available ... "
        print response_read
        output.write("\t'%s':%.2f,\n"%(name, -1.))
        output.write("\t'%s_rms':%.2f,\n"%(name, -1.))
        continue

    json_data = json.loads(response_read)
    data = np.array(json_data['records']) 
    #--> returns a 2d table of [timestamp in ms, value]
    """converts timestamp to s"""
    data[:,0] *= 1.e-3 
    

    time_sel = (data[:,0]>time_start) & (data[:,0] < time_stop)
    """selects only times during data taking"""

    last_meas = np.where(data[:,0]>time_start)
    """last point measured before run"""

    if(len(last_meas[0]) == 0 and np.count_nonzero(time_sel) > 0):
        print "no measurement available ... "
        output.write("\t'%s':%.2f,\n"%(name, -1.))
        output.write("\t'%s_rms':%.2f,\n"%(name, -1.))
        continue
    
    if(len(last_meas[0]) == 0):
        print " no prior measurements available "
        output.write("\t'%s':%.2f,\n"%(name, -1.))
        output.write("\t'%s_rms':%.2f,\n"%(name, -1.))
        continue
    else:
        last = data[ last_meas[0][0]-1]        
        print " --> last point before run :", last[1], " (%.2f mn ago)"% ((time_start - last[0])/60.)
        val = last[1]
        val_rms = 0.

        if(name.find("lem") >= 0):
            output.write("\t'%s':%.2f,\n"%(name, last[1]))
            output.write("\t'%s_rms':%.2f,\n"%(name, 0.))

        if(name.find("grid") >= 0):
            output.write("\t'%s':%.2f,\n"%(name, last[1]))
            output.write("\t'%s_rms':%.2f,\n"%(name, 0.))

        if(np.count_nonzero(time_sel) > 0 and name.find("lem") < 0): 
            print " during run : "
            val = data[ time_sel ].mean(axis=0)[1]
            val_rms  = data[ time_sel].std(axis=0)[1]
            print " --> mean : ", val
            print " --> std : ", val_rms

        if(name.find("lem") <= 0): #not lem
            output.write("\t'%s':%.2f,\n"%(name, val))
            output.write("\t'%s_rms':%.2f,\n"%(name, val_rms))
                
        if(np.count_nonzero(time_sel) > 0 and name.find("lem") >= 0): #lem
            print " during run : "
            val = data[ time_sel ].mean(axis=0)[1]
            val_rms  = data[ time_sel].std(axis=0)[1]
            print " --> mean : ", val
            print " --> std : ", val_rms
            plt.plot(data[:,0], data[:,1],'k')
            plt.axvspan(time_start, time_stop, color='red', alpha=0.5)
            plt.savefig("run_"+str(run_nb)+"_"+name+".png")
            plt.close()


            data_z = data[time_sel]
            deriv = data_z[:,1]-last[1]
            deriv /= last[1]
            diff  = np.diff(data_z[:,1])
                
            idx = 0
                
            """naive spark hunter"""
            while(idx < len(deriv)-1):   
                if(deriv[idx] < -0.2):
                    tspark_start.append(data_z[idx,0])
                    spark_id.append(name)

                    while(idx < len(deriv)-1 and deriv[idx] < -0.2 and diff[idx] < 0):
                        idx += 1
                            
                    tspark_max.append(data_z[idx,0])
                        
                    
                    while(idx < len(deriv)-1 and deriv[idx] < -0.2 and diff[idx] > 0):
                        idx += 1

                    tspark_end.append(data_z[idx,0])

                else:
                    idx += 1

            
            if( (len(tspark_start) != len(tspark_max)) or (len(tspark_start) != len(tspark_end))):
                print "spark algo had issues ... "
                continue




nsparks = len(tspark_start)
output.write("\t'nsparks':%i,\n"%nsparks)
print " found ", nsparks
for isp in range(nsparks):
    print " spark ", isp, " : on ", spark_id[isp], " at ", tspark_start[isp], " - ", tspark_max[isp], " - ", tspark_end[isp]
    output.write("\t'spark_%i_id':'%s',\n"%(isp,spark_id[isp]))
    output.write("\t'spark_%i_time':%.0f,\n"%(isp,tspark_start[isp]))
    
output.write("}\n")
output.close()
