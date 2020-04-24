import config as cf
from tables import *
import numpy as np
import data_containers as dc

class Infos(IsDescription):
    run = UInt16Col()
    subfile = StringCol(8)
    nEvent = UInt8Col()
    process_date = UInt32Col()

class Event(IsDescription):
    evt_nb_global = UInt32Col()
    time_s        = UInt32Col()
    time_ns       = UInt32Col()
    nHits         = UInt32Col(shape=(cf.n_CRP,cf.n_View))
    nClusters     = UInt32Col(shape=(cf.n_CRP,cf.n_View))
    nTracks2D     = UInt32Col(shape=(cf.n_View))
    nTracks3D     = UInt32Col()

class Pedestal(IsDescription):
    crp       = UInt8Col()
    view      = UInt8Col()
    channel   = UInt16Col()
    raw_mean  = Float16Col()
    raw_rms   = Float16Col()
    filt_mean = Float16Col()
    filt_rms  = Float16Col()


class Hits(IsDescription):
    crp     = UInt8Col()
    view    = UInt8Col()
    channel = UInt16Col()
    tdc_max = UInt16Col()
    z       = Float16Col()
    x       = Float16Col()
    dt      = Float16Col()
    adc_max = UInt16Col()
    charge  = Float16Col()
    cluster = Int16Col()


class Tracks2D(IsDescription):
    view    = UInt8Col()
    crp_ini = UInt8Col()
    crp_end = UInt8Col()
    pos_ini = Float16Col()
    pos_end = Float16Col()
    z_ini   = Float16Col()
    z_end   = Float16Col()
    nHits   = UInt16Col()
    chi2    = Float16Col()
    slope_ini = Float16Col()
    slope_end = Float16Col()
    len_straight = Float16Col()
    len_path     = Float16Col()
    total_charge = Float16Col()

class Tracks3D(IsDescription):
    crp_ini = UInt8Col()
    crp_end = UInt8Col()
    x_ini   = Float16Col()
    y_ini   = Float16Col()
    z_ini   = Float16Col()
    x_end   = Float16Col()
    y_end   = Float16Col()
    z_end   = Float16Col()
    chi2    = Float16Col()
    theta_ini = Float16Col()
    theta_end = Float16Col()
    phi_ini   = Float16Col()
    phi_end   = Float16Col()

    nHits        = UInt16Col(shape=(cf.n_View))
    len_straight = Float16Col(shape=(cf.n_View))
    len_path     = Float16Col(shape=(cf.n_View))
    total_charge   = Float16Col(shape=(cf.n_View))
    
def new_event(h5file, event_nb):
    return h5file.create_group("/", 'event_'+str(event_nb), 'Event '+str(event_nb))    

def store_infos(h5file, run, subfile, nevt, time):
    table = h5file.create_table("/", 'infos', Infos, 'Infos')
    inf = table.row
    inf['run'] = run

    inf['subfile'] = subfile.ljust(8)[:8] #so the string is exactly 8 caracters
    print(subfile)
    print(subfile.ljust(8)[:8])
    inf['nEvent'] = nevt
    inf['process_date'] = time
    inf.append()
    table.flush()

def store_event(h5file, group):
    table = h5file.create_table(group, 'event', Event, "Event")

    evt = table.row

    evt['evt_nb_global'] = dc.evt_list[-1].evt_nb_glob
    evt['time_s']        = dc.evt_list[-1].time_s
    evt['time_ns']       = dc.evt_list[-1].time_ns
    evt['nHits']         = dc.evt_list[-1].nHits
    evt['nClusters']     = dc.evt_list[-1].nClusters
    evt['nTracks2D']     = dc.evt_list[-1].nTracks2D
    evt['nTracks3D']     = dc.evt_list[-1].nTracks3D
    evt.append()
    table.flush()


def store_pedestal(h5file, group):
    table = h5file.create_table(group, 'pedestals', Pedestal, 'Pedestals')

    ped = table.row

    for x in dc.map_ped:
        ped['crp']       = x.crp
        ped['view']      = x.view
        ped['channel']   = x.vchan
        ped['raw_mean']  = x.raw_ped
        ped['raw_rms']   = x.raw_rms
        ped['filt_mean'] = x.evt_ped
        ped['filt_rms']  = x.evt_rms
        ped.append()
    table.flush()


def store_hits(h5file, group):
    table = h5file.create_table(group, 'hits', Hits, 'Hits')

    hit = table.row
    for h in dc.hits_list:
        hit['crp']     = h.crp
        hit['view']    = h.view
        hit['channel'] = h.channel
        hit['tdc_max'] = h.max_t
        hit['z']       = h.Z
        hit['x']       = h.X
        hit['dt']      = (h.stop - h.start)*cf.n_Sampling
        hit['adc_max'] = h.max_adc
        hit['charge']  = h.charge
        hit['cluster'] = h.cluster
        hit.append()
    table.flush()


def store_tracks2D(h5file, group):    
    table = h5file.create_table(group, 'tracks2D', Tracks2D, "Tracks 2D")       

    t2d_hits_v0 = h5file.create_group(group, 'tracks2D_v0', 'Tracks 2D View0')
    t2d_hits_v1 = h5file.create_group(group, 'tracks2D_v1', 'Tracks 2D View1')

    t2d = table.row
    i_v0 = 0
    i_v1 = 0
    for t in dc.tracks2D_list:
        t2d['view']      = t.view
        t2d['crp_ini']   = t.ini_crp
        t2d['crp_end']   = t.end_crp
        t2d['pos_ini']   = t.path[0][0]
        t2d['z_ini']     = t.path[0][1]
        t2d['pos_end']   = t.path[-1][0]
        t2d['z_end']     = t.path[-1][1]
        t2d['nHits']     = t.nHits
        t2d['slope_ini'] = t.ini_slope
        t2d['slope_end'] = t.end_slope
        t2d['chi2']      = t.chi2
        t2d['len_straight'] = t.len_straight
        t2d['len_path']     = t.len_path
        t2d['total_charge'] = t.tot_charge
        pts = [[p[0], p[1], q] for p,q in zip(t.path,t.dQ)]
        if(t.view==0):
            h5file.create_array(t2d_hits_v0, 'track_%i'%(i_v0), np.asarray(pts), 'track hits')
            i_v0 += 1
        else:
            h5file.create_array(t2d_hits_v1, 'track_%i'%(i_v1), np.asarray(pts), 'track hits')
            i_v1 += 1

        t2d.append()
    table.flush()



def store_tracks3D(h5file, group):
    table = h5file.create_table(group, 'tracks3D', Tracks3D, "Tracks 3D")       

    t3d_hits_v0 = h5file.create_group(group, 'tracks3D_v0', 'Tracks 3D View0')
    t3d_hits_v1 = h5file.create_group(group, 'tracks3D_v1', 'Tracks 3D View1')

    t3d = table.row
    i = 0

    for t in dc.tracks3D_list:
        t3d['crp_ini']   = t.ini_crp
        t3d['crp_end']   = t.end_crp

        t3d['x_ini']     = t.ini_x
        t3d['y_ini']     = t.ini_y
        t3d['z_ini']     = t.ini_z
        t3d['x_end']     = t.end_x
        t3d['y_end']     = t.end_y
        t3d['z_end']     = t.end_z

        t3d['chi2']      = t.chi2
        t3d['nHits']     = [t.nHits_v0, t.nHits_v1]
        t3d['len_straight'] = [t.len_straight_v0, t.len_straight_v1]
        t3d['len_path']     = [t.len_path_v0, t.len_path_v1]
        t3d['total_charge'] = [t.tot_charge_v0, t.tot_charge_v1]


        t3d['theta_ini'] = t.ini_theta
        t3d['theta_end'] = t.end_theta
        t3d['phi_ini']   = t.ini_phi
        t3d['phi_end']   = t.end_phi

        pts_v0 = [[p[0], p[1], p[2], q] for p,q in zip(t.path_v0,t.dQds_v0)]
        pts_v1 = [[p[0], p[1], p[2], q] for p,q in zip(t.path_v1,t.dQds_v1)]
        

        h5file.create_array(t3d_hits_v0, 'track_%i'%(i), np.asarray(pts_v0), 'track hits')
        h5file.create_array(t3d_hits_v1, 'track_%i'%(i), np.asarray(pts_v1), 'track hits')
        i += 1

        t3d.append()
    table.flush()
