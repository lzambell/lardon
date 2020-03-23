import config as cf
import struct
import numpy as np

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

def read_event(data, idx):
    data.seek(idx,0)
    header_size = header_type.itemsize

    v0, lro, cro = read_event_header( data.read(header_size) )

    if(cro < 0):
        return v0, [], []

    if(lro > 0): #in case one day light readout is also written in the data file
        data.read(lro, 1)

    npv0 = read_evt_uint32( data.read(cro))

    data.read(1) #BRUNO BYTE (?)
        
    v1, lro, cro = read_event_header( data.read(header_size))

    if(cro < 0):
        return v1, [], []

    if(lro > 0):
        data.read(lro, 1)

    npv1 = read_evt_uint32( data.read(cro))
    v = check_and_merge_events(v0, v1)
    return v, npv0, npv1
    
def read_event_header(data):
    head = np.frombuffer(data, dtype=header_type)
            
    return check_and_build_event(head)

        
def read_evt_uint32(data):
    #solution from 
    #https://stackoverflow.com/questions/44735756/python-reading-12-bit-binary-files
    #could be faster with numba


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
    print(" evt_flag is ", header['evt_flag'][0])
    if(not good_evt):
        print(" problem, the event is marked as bad -->")
        #return throw_bad_event(), -1, -1

    ev = cf.event(header['run_num'][0], header['evt_num'][0], header['time_s'][0], header['time_ns'][0], good_evt)

    return ev, header['lro'][0], header['cro'][0]


def check_and_merge_events(v0, v1):
    if(v0 == v1 and v0.evt_flag == True): #in case both are bad events
        return v0
    else:
        print("Event headers are different ... ")
        """
        print " view 0 header : "
        dump_header(v0)
        print " view 1 header : "
        dump_header(v1)
        #return throw_bad_event()
        """
        return v0
def throw_bad_event():
    return cf.event(-1, -1, -1, -1, False)

def dump_header(header):
    print(" run_nb : ", header.run_nb)
    print(" evt_nb_loc : ", header.evt_nb_loc)
    print(" evt_nb_glob : ", header.evt_nb_glob)
    print(" time_s :", header.time_s)
    print(" time_ns :", header.time_ns)
    print(" evt_flag :", header.evt_flag)



"""
#Slow - kept just in case
def read_event_header(data, idx):
    data.seek(idx,0)
    fmt_control = 'BB'
    size_control = struct.calcsize(fmt_control)
    k0,k1 = struct.unpack('BB', data.read( size_control ))
        
    if( not((k0 & 0xFF)==evskey and (k1 & 0xFF)==evskey)):
        print " problem "
    else:
        print " ok ! "
        
    fmt_runevent = 'Ic'
    size_runevent = struct.calcsize(fmt_runevent)
    runnum, runflag = struct.unpack(fmt_runevent, data.read( size_runevent ))

    fmt_time= '<B3cIqq'
    size_time = struct.calcsize(fmt_time)

    trigtype, tmp, tmp, tmp, trignum, time_s, time_ns = struct.unpack(fmt_time, data.read( size_time ))
    print " trigger type: ", trigtype, " trigger nb ", trignum, " at ", time_s, "s +", time_ns, " ns"
    
    fmt_event = '<BIII'
    size_event = struct.calcsize(fmt_event)
    
    evflag, evnum, lro, cro = struct.unpack(fmt_event, data.read( size_event ))
    good_evt = (evflag & 0x3F)==evdcard0
    print " Good Event ? ", good_evt
    
    print evnum
    print lro
    print cro
    
    ev = cf.event()
    ev.run_nb = runnum
    ev.evt_nb_glob = evnum
    ev.time_s = time_s
    ev.time_ns = time_ns
    ev.evt_flag = good_evt

    return ev, lro, cro
"""
