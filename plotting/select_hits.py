import config as cf
import data_containers as dc


def get_hits_pos(view, selection="True"):
    return [x.X for x in dc.hits_list if x.view == view and eval(selection)]


def get_hits_ch(view, selection="True"):
    return [x.channel for x in dc.hits_list if x.view == view and eval(selection)]

def get_hits_tdc(view, selection="True"):
    return [x.max_t for x in dc.hits_list if x.view == view and eval(selection)]

def get_hits_z(view, selection="True"):
    return [x.Z for x in dc.hits_list if x.view == view and eval(selection)]

def get_hits_cluster(view, selection="True"):
    return [x.cluster for x in dc.hits_list if x.view == view and eval(selection)]

def get_hits_charge(view, selection="True"):
    return [x.charge for x in dc.hits_list if x.view == view and eval(selection)]

def get_hits_adc(view, selection="True"):
    return [x.max_adc for x in dc.hits_list if x.view == view and eval(selection)]


def get_2dtracks_pos(view, selection="True"):
    return [[p[0] for p in t.path] for t in dc.tracks2D_list if t.view==view and eval(selection)]

def get_2dtracks_z(view, selection="True"):
    return [[p[1] for p in t.path] for t in dc.tracks2D_list if t.view==view and eval(selection)]


def get_3dtracks(view, axis, selection="True"):
    if(view==0):
        return [[p[axis] for p in t.path_v0] for t in dc.tracks3D_list if eval(selection)]
    elif(view==1):
        return [[p[axis] for p in t.path_v1] for t in dc.tracks3D_list if eval(selection)] 
    else:
        raise ValueError



def get_3dtracks_x(view, selection="True"):
    return get_3dtracks(view, 0, selection)

def get_3dtracks_y(view, selection="True"):
    return get_3dtracks(view, 1, selection)

def get_3dtracks_z(view, selection="True"):
    return get_3dtracks(view, 2, selection)
