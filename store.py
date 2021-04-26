import config as cf
from tables import *
import numpy as np
import data_containers as dc

class Infos(IsDescription):
    run          = UInt16Col()
    subfile      = StringCol(8)
    nEvent       = UInt8Col()
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
    raw_mean  = Float32Col()
    raw_rms   = Float32Col()
    filt_mean = Float32Col()
    filt_rms  = Float32Col()


class Hits(IsDescription):
    crp     = UInt8Col()
    view    = UInt8Col()
    channel = UInt16Col()
    tdc_max = UInt16Col()
    z       = Float32Col()
    x       = Float32Col()
    dt      = Float32Col()
    adc_max = Float32Col()
    charge  = Float32Col()
    cluster = Int16Col()
    z_start = Float32Col()
    z_stop  = Float32Col()

class Tracks2D(IsDescription):
    view    = UInt8Col()
    crp_ini = UInt8Col()
    crp_end = UInt8Col()
    pos_ini = Float32Col()
    pos_end = Float32Col()
    z_ini   = Float32Col()
    z_end   = Float32Col()
    nHits   = UInt16Col()
    chi2    = Float32Col()

    slope_ini = Float32Col()
    slope_end = Float32Col()

    len_straight = Float32Col()
    len_path     = Float32Col()
    total_charge = Float32Col()

class Tracks3D(IsDescription):
    crp_ini = UInt8Col()
    crp_end = UInt8Col()
    x_ini   = Float32Col()
    y_ini   = Float32Col()
    z_ini   = Float32Col()
    x_end   = Float32Col()
    y_end   = Float32Col()
    z_end   = Float32Col()
    chi2    = Float32Col()

    theta_ini = Float32Col()
    theta_end = Float32Col()
    phi_ini   = Float32Col()
    phi_end   = Float32Col()

    nHits        = UInt16Col(shape=(cf.n_View))
    len_straight = Float32Col(shape=(cf.n_View))
    len_path     = Float32Col(shape=(cf.n_View))
    total_charge   = Float32Col(shape=(cf.n_View))

    len_straight_field_corr = Float32Col(shape=(cf.n_View))
    len_path_field_corr     = Float32Col(shape=(cf.n_View))
    total_charge_field_corr   = Float32Col(shape=(cf.n_View))


    z0_corr = Float32Col()
    t0_corr = Float32Col()
    
def new_event(h5file, event_nb):
    return h5file.create_group("/", 'event_'+str(event_nb), 'Event '+str(event_nb))    

def store_infos(h5file, run, subfile, nevt, time):
    table = h5file.create_table("/", 'infos', Infos, 'Infos')
    inf = table.row
    inf['run'] = run

    inf['subfile'] = subfile.ljust(8)[:8] #so the string is exactly 8 caracters
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
        hit['z_start'] = h.Z_start
        hit['z_stop']  = h.Z_stop
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

        t3d['len_straight_field_corr'] = [t.len_straight_field_corr_v0, t.len_straight_field_corr_v1]
        t3d['len_path_field_corr']     = [t.len_path_field_corr_v0, t.len_path_field_corr_v1]
        t3d['total_charge_field_corr'] = [t.tot_charge_field_corr_v0, t.tot_charge_field_corr_v1]


        t3d['theta_ini'] = t.ini_theta
        t3d['theta_end'] = t.end_theta
        t3d['phi_ini']   = t.ini_phi
        t3d['phi_end']   = t.end_phi

        t3d['z0_corr']    = t.z0_corr
        t3d['t0_corr']    = t.t0_corr

        pts_v0 = [[p[0], p[1], p[2], q[0]/q[1], pc, qc[0]/qc[1]] for p,q,pc,qc in zip(t.path_v0, t.dQds_v0, t.z_field_corr_v0, t.dQds_field_corr_v0)]
        pts_v1 = [[p[0], p[1], p[2], q[0]/q[1], pc, qc[0]/qc[1]] for p,q,pc,qc in zip(t.path_v1, t.dQds_v1, t.z_field_corr_v1, t.dQds_field_corr_v1)]
        

        h5file.create_array(t3d_hits_v0, 'track_%i'%(i), np.asarray(pts_v0), 'track hits and charge (x, y, z, dqdx, zcorr, dqdx corr)')
        h5file.create_array(t3d_hits_v1, 'track_%i'%(i), np.asarray(pts_v1), 'track hits and charge (x, y, z, dqdx, zcorr, dqdx corr)')
        i += 1

        t3d.append()
    table.flush()



def create_lite(h5file):
    table = h5file.create_table("/", 'tracks3D', Tracks3D, "Tracks 3D")           
    t = h5file.create_vlarray("/", 'pathv0', Float32Atom(shape=(6)), "Path V0 (x, y, z, q, zcorr, qcorr)")    
    t = h5file.create_vlarray("/", 'pathv1', Float32Atom(shape=(6)), "Path V1 (x, y, z, q, zcorr, qcorr)")    
       

def store_tracks3D_lite(h5file):    
    t3d = h5file.root.tracks3D.row
    vlv0  = h5file.root.pathv0
    vlv1  = h5file.root.pathv1

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


        t3d['len_straight_field_corr'] = [t.len_straight_field_corr_v0, t.len_straight_field_corr_v1]
        t3d['len_path_field_corr']     = [t.len_path_field_corr_v0, t.len_path_field_corr_v1]
        t3d['total_charge_field_corr'] = [t.tot_charge_field_corr_v0, t.tot_charge_field_corr_v1]


        t3d['theta_ini'] = t.ini_theta
        t3d['theta_end'] = t.end_theta
        t3d['phi_ini']   = t.ini_phi
        t3d['phi_end']   = t.end_phi

        t3d['z0_corr']    = t.z0_corr
        t3d['t0_corr']    = t.t0_corr

        pts_v0 = [[p[0], p[1], p[2], q[0]/q[1], pc, qc[0]/qc[1]] for p,q,pc,qc in zip(t.path_v0, t.dQds_v0, t.z_field_corr_v0, t.dQds_field_corr_v0)]
        pts_v1 = [[p[0], p[1], p[2], q[0]/q[1], pc, qc[0]/qc[1]] for p,q,pc,qc in zip(t.path_v1, t.dQds_v1, t.z_field_corr_v1, t.dQds_field_corr_v1)]

        t3d.append()
        vlv0.append(pts_v0)
        vlv1.append(pts_v1)
    #table.flush()
    
