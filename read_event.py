import config as cf
import data_containers as dc


import numpy as np
import numba as nb

evskey = 0xFF
endkey = 0xF0
evdcard0 = 0x19


header_type = np.dtype([
    ('k0','B'),
    ('k1','B'),
    ('run_num', '<u4'), 
    ('run_flag', 'c'),
    ('trig_type', '<B'),
    ('padding', '<3c'),
    ('trig_num', '<u4'),
    ('time_s', '<i8'),
    ('time_ns', '<i8'), 
    ('evt_flag', '<B'), 
    ('evt_num', '<u4'),
    ('lro', '<u4'),
    ('cro', '<u4')
])

def read_event(data, idx, iev):
    data.seek(idx,0)
    header_size = header_type.itemsize

    v0, lro, cro = read_event_header( data.read(header_size) )

    if(cro < 0):
        return -1 #v0, [], []

    """in case one day light readout is also written in the data file"""
    if(lro > 0): 
        data.read(lro, 1)


    shape_and_store( read_evt_uint12_nb( data.read(cro)) , 0)

    data.read(1) #BRUNO BYTE (?)
        
    v1, lro, cro = read_event_header( data.read(header_size))

    if(cro < 0):
        return -1 #v1, [], []

    if(lro > 0):
        data.read(lro, 1)

    shape_and_store( read_evt_uint12_nb( data.read(cro)) , 1)

    check_and_merge_events(v0, v1, iev)
    return 0
    
def read_event_header(data):
    head = np.frombuffer(data, dtype=header_type)
            
    return check_and_build_event(head)


def shape_and_store(data_raw, vread):

    if(len(data_raw)/cf.n_Sample != cf.n_ChanPerView):
        print(" PBM OF Nb of CHANNELS in View", view, "!!! ", len(data_raw)/cf.n_Sample , " vs ", cf.n_ChanPerView)
        return
    
    data_raw = np.split(data_raw, cf.n_ChanPerView)

    shift = 0
    if(vread == 1):
        shift = cf.n_ChanPerView

    """reshape the array and subtract reference pedestal"""    
    for idq in range(cf.n_ChanPerView):

        crp, view, vch = dc.map_ped[idq+shift].get_ana_chan()

        if(view != vread): continue
        if(crp >= cf.n_CRPUsed): continue 
        pedval = dc.map_ped[idq+shift].ref_ped 

        if(crp < 0 or view < 0 or vch < 0):
            print(" ERROR ? ", idq)

        dc.data[crp,view,vch] = data_raw[idq] - pedval

    

@nb.jit
def read_evt_uint12_nb(data):
    tt = np.frombuffer(data, dtype=np.uint8)
    assert np.mod(tt.shape[0],3)==0

    out=np.empty(tt.shape[0]//3*2,dtype=np.uint16)

    for i in nb.prange(tt.shape[0]//3):
        fst_uint8=np.uint16(data[i*3])
        mid_uint8=np.uint16(data[i*3+1])
        lst_uint8=np.uint16(data[i*3+2])

        out[i*2]   = (fst_uint8 << 4) + (mid_uint8 >> 4)
        out[i*2+1] = ((mid_uint8 % 16) << 8) + lst_uint8

    return out

        
def read_evt_uint12(data):
    #solution from 
    #https://stackoverflow.com/questions/44735756/python-reading-12-bit-binary-files

    tt = np.frombuffer(data, dtype=np.uint8)
    
    fst_uint8, mid_uint8, lst_uint8 = np.reshape(tt, (tt.shape[0] // 3, 3)).astype(np.uint16).T
    fst_uint12 = (fst_uint8 << 4) + (mid_uint8 >> 4)
    snd_uint12 = ((mid_uint8 % 16) << 8) + lst_uint8
    return np.reshape(np.concatenate((fst_uint12[:, None], snd_uint12[:, None]), axis=1), 2 * fst_uint12.shape[0])



def check_and_build_event(header):
    if( not((header['k0'][0] & 0xFF)==evskey and (header['k1'][0] & 0xFF)==evskey)):
        print(" problem ")
        return throw_bad_event(), -1, -1

    good_evt = (header['evt_flag'][0] & 0x3F) == evdcard0

    #print("header event flag is ", header['evt_flag'][0])
    if(not good_evt):
        print(" problem, the event is marked as bad ... ")
        #return throw_bad_event(), -1, -1

    ev = dc.event(header['run_num'][0], header['evt_num'][0], header['time_s'][0], header['time_ns'][0], good_evt)

    return ev, header['lro'][0], header['cro'][0]


def check_and_merge_events(v0, v1, iev):
    if(v0 == v1 and v0.evt_flag == True): #in case both are bad events
        v0.evt_nb_loc = iev
        dc.evt_list.append(v0)
        #return v0
    else:
        print("Event headers are different ... ")
        v0.evt_nb_loc = iev
        dc.evt_list.append(v0)

        #return v0



def throw_bad_event():
    return dc.event(-1, -1, -1, -1, False)

def dump_header(header):
    print(" run_nb : ", header.run_nb)
    print(" evt_nb_loc : ", header.evt_nb_loc)
    print(" evt_nb_glob : ", header.evt_nb_glob)
    print(" time_s :", header.time_s)
    print(" time_ns :", header.time_ns)
    print(" evt_flag :", header.evt_flag)
