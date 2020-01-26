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
    #at = 0
    v0, lro, cro = read_event_header_np( data.read(header_size) )
    #at += header_size
    npv0 = read_evt_uint32( data.read(cro))

    data.read(1) #BRUNO BYTE
    
    
    v1, lro, cro = read_event_header_np( data.read(header_size))
    #at += header_size
    npv1 = read_evt_uint32( data.read(cro))
    v = v0
    return v, npv0, npv1
    
def read_event_header_np(data):#, idx):
    #data.seek(idx,0)
    #head = np.frombuffer(data.read(header_type.itemsize), dtype=header_type)
    head = np.frombuffer(data, dtype=header_type)
    
    if( not((head['k0'][0] & 0xFF)==evskey and (head['k1'][0] & 0xFF)==evskey)):
        print " problem " ## AND SHOULD DO SOMETHING
    good_evt = (head['evt_flag'][0] & 0x3F)==evdcard0

    ev = cf.event()
    ev.run_nb = head['run_num'][0]
    ev.evt_nb_glob = head['evt_num'][0]
    ev.time_s = head['time_s'][0]
    ev.time_ns = head['time_ns'][0]
    ev.evt_flag = good_evt

    return ev, head['lro'][0], head['cro'][0]
        
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

def read_evt_uint32(data):
    #solution from 
    #https://stackoverflow.com/questions/44735756/python-reading-12-bit-binary-files
    #could be faster with numba


    tt = np.frombuffer(data, dtype=np.uint8)
    
    fst_uint8, mid_uint8, lst_uint8 = np.reshape(tt, (tt.shape[0] // 3, 3)).astype(np.uint16).T
    fst_uint12 = (fst_uint8 << 4) + (mid_uint8 >> 4)
    snd_uint12 = ((mid_uint8 % 16) << 8) + lst_uint8
    return np.reshape(np.concatenate((fst_uint12[:, None], snd_uint12[:, None]), axis=1), 2 * fst_uint12.shape[0])
