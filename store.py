from tables import *
import numpy as np
import data_containers as dc

class Tracks(IsDescription):
    view = UInt8Col()
    crp_ini = UInt8Col()
    crp_end = UInt8Col()
    pos_ini = Int32Col()
    pos_end = Int32Col()
    z_ini = Int32Col()
    z_end = Int32Col()
    nHits = UInt16Col()



def store_tracks(h5file, event_nb):
    #h5file = open_file(name+".h5", mode="w", title="Test file")
    group = h5file.create_group("/", 'event_'+str(event_nb), 'Event '+str(event_nb))
    table = h5file.create_table(group, 'tracks', Tracks, "Tracks 2D")       

    t2d = table.row

    for t in dc.tracks2D_list:
        #t2d['ID'] = i
        t2d['view'] = t.view
        t2d['crp_ini'] = t.ini_crp
        t2d['crp_end'] = t.end_crp
        t2d['pos_ini']  = t.path[0][0]
        t2d['z_ini']  = t.path[0][1]
        t2d['pos_end']  = t.path[-1][0]
        t2d['z_end']  = t.path[-1][1]
        t2d['nHits']  = t.nHits
        t2d.append()
    table.flush()
    #h5file.close()
