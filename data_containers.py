import config as cf
import numpy as np
import math
import time

map_ped = []
evt_list = []
hits_list = []
tracks2D_list = []
tracks3D_list = []

data = np.zeros((cf.n_CRPUsed, cf.n_View, cf.n_ChanPerCRP, cf.n_Sample)) #crp, view, vchan


"""
the mask will be used to differentiate background (True for noise processing) from signal (False for noise processing)
at first everything is considered background (all at True)
"""
mask = np.ones((cf.n_CRPUsed, cf.n_View, cf.n_ChanPerCRP, cf.n_Sample), dtype=bool)

"""
alive_chan mask intends to not take into account broken channels
True : not broken
False : broken
"""
alive_chan = np.ones((cf.n_CRPUsed,cf.n_View, cf.n_ChanPerCRP, cf.n_Sample), dtype=bool)


ped_rms = np.zeros((cf.n_CRPUsed, cf.n_View, cf.n_ChanPerCRP)) 
ped_mean = np.zeros((cf.n_CRPUsed, cf.n_View, cf.n_ChanPerCRP)) 

def reset_event():
    data[:,:,:,:] = 0.
    mask[:,:,:,:] = True
    ped_rms[:,:,:] = 0.
    ped_mean[:,:,:] = 0.

    hits_list.clear()
    tracks2D_list.clear()
    tracks3D_list.clear()

    [x.reset() for x in map_ped]


class pdmap: 
    def __init__(self, crp, view, vchan):
        self.crp   = crp
        self.view  = view
        self.vchan = vchan
        self.ref_ped   = 0.
        self.ref_rms   = -1.
        self.raw_ped   = 0.
        self.raw_rms   = -1.
        self.evt_ped   = 0.
        self.evt_rms   = -1.
        
    def set_ref_pedestal(self, ped, rms):
        self.ref_ped = ped
        self.ref_rms = rms

    def set_raw_pedestal(self, ped, rms):
        self.raw_ped = ped
        self.raw_rms = rms

    def set_evt_pedestal(self, ped, rms):
        self.evt_ped = ped
        self.evt_rms = rms
        
    def reset(self):
        self.raw_ped = 0.
        self.raw_rms = -1.
        self.evt_ped = 0.
        self.evt_rms = -1.


    def get_ana_chan(self):
        return self.crp, self.view, self.vchan

class event:
    def __init__(self, run_nb, evt_glob, t_s, t_ns, flag):
        self.run_nb      = run_nb
        self.evt_nb_glob = evt_glob
        self.evt_nb_loc  = -1
        self.time_s      = t_s
        self.time_ns     = t_ns
        self.evt_flag    = flag
        self.nHits       = np.zeros((cf.n_CRP, cf.n_View), dtype=int)
        self.nClusters   = np.zeros((cf.n_CRP, cf.n_View), dtype=int)
        self.nTracks2D   = np.zeros((cf.n_View), dtype=int)
        self.nTracks3D   = 0
        
    def __eq__(self, other):
        return (self.run_nb, self.evt_nb_glob, self.evt_nb_loc, self.time_s, self.time_ns, self.evt_flag) == (other.run_nb, other.evt_nb_glob, other.evt_nb_loc, other.time_s, other.time_ns, other.evt_flag)

    def dump(self):
        print("RUN ",self.run_nb, " EVENT ", self.evt_nb_loc, " / ", self.evt_nb_glob,)
        print("Taken at ", time.ctime(self.time_s), " + ", self.time_ns, " ns ")
        print(" TS = ", self.time_s, " +", self.time_ns)
    def dump_reco(self):
        print("\n*-*-*-* Reconstruction Summary *-*-*-*")
        for icrp in range(cf.n_CRPUsed):
            for iv in range(cf.n_View):
                print("CRP ", icrp, " View ", iv, " : ")
                print("\tNb of Hits :", self.nHits[icrp,iv])
                print("\tNb of Clusters :", self.nClusters[icrp,iv])
        print("\n")
        for iv in range(cf.n_View):
            print("View ", iv)
            print("\tNb of 2D tracks :", self.nTracks2D[iv])

        print("\n")
        print("--> Nb of 3D tracks ", self.nTracks3D)
        print("\n")


class hits:
    def __init__(self, crp, view, channel, start, stop, charge, max_t, max_adc):
        self.crp     = crp
        self.view    = view
        self.channel = channel
        self.start   = start
        self.stop    = stop
        self.charge  = charge
        self.max_t   = max_t
        self.max_adc = max_adc
        self.cluster = -1 #cluster
        self.X       = -1
        self.Z       = -1
        self.matched = -1


    def __lt__(self,other):
        #"""sort hits by increasing channel and increasing Z"""
        #return (self.X < other.X) or (self.X== other.X and self.Z < other.Z)

        """ sort hits by decreasing Z and increasing channel """
        return (self.Z > other.Z) or (self.Z == other.Z and self.X < other.X)

    def hit_positions(self, v):
        self.X = self.channel*cf.ChanPitch

        if(self.view == 0):
            if(self.crp == 1 or self.crp == 2):
                self.X -= cf.n_ChanPerCRP * cf.ChanPitch
        
        else:
            if(self.crp == 2 or self.crp == 3):
                self.X -= cf.n_ChanPerCRP * cf.ChanPitch                

        self.Z = cf.Anode_Z - self.max_t*cf.n_Sampling*v*0.1

    def hit_charge(self):
        self.charge *= cf.n_Sampling 
        self.charge /= cf.ADCtofC

    def set_match(self, ID):
        self.matched=ID

class trk2D:
    def __init__(self, ini_crp, view, ini_slope, ini_slope_err, x0, y0, q0, chi2, cluster):
        self.ini_crp = ini_crp
        self.end_crp = ini_crp
        self.view    = view
        self.ini_slope       = ini_slope
        self.ini_slope_err   = ini_slope_err
        self.end_slope       = ini_slope
        self.end_slope_err   = ini_slope_err
        self.nHits   = 1
        self.path    = [(x0,y0)]
        self.dQ      = [q0]
        self.chi2    = chi2

        self.tot_charge = q0
        self.len_straight = 0.
        self.len_path = 0.
        self.matched = -1
        self.cluster = cluster
        
    def __lt__(self,other):
        """ sort tracks by decreasing Z and increasing channel """
        return (self.path[0][1] > other.path[0][1]) or (self.path[0][1] == other.path[0][1] and self.path[0][0] < other.path[0][0])


    def add_hit(self, x, y, q):
        self.nHits += 1
        self.len_path += math.sqrt( pow(self.path[-1][0]-x, 2) + pow(self.path[-1][1]-y,2) )
        #beware to append (x,y) after !
        self.path.append((x,y))
        self.dQ.append(q)
        self.tot_charge += q
        self.len_straight = math.sqrt( pow(self.path[0][0]-self.path[-1][0], 2) + pow(self.path[0][1]-self.path[-1][1], 2) )


    def add_hit_update(self, slope, slope_err, x, y, q, chi2):
        self.end_slope = slope
        self.end_slope_err = slope_err
        self.nHits += 1
        self.len_path += math.sqrt( pow(self.path[-1][0]-x, 2) + pow(self.path[-1][1]-y,2) )

        #beware to append (x,y) after !
        self.path.append((x,y))
        self.dQ.append(q)
        self.chi2 = chi2
        self.tot_charge += q
        self.len_straight = math.sqrt( pow(self.path[0][0]-self.path[-1][0], 2) + pow(self.path[0][1]-self.path[-1][1], 2) )
        

    def dist(self, other):
        return math.sqrt(pow( self.path[-1][0] - other.path[0][0], 2) + pow(self.path[-1][1] - other.path[0][1], 2))

    def slope_comp(self, other, sigcut):
        return (math.fabs( self.end_slope - other.ini_slope) < (sigcut*self.end_slope_err + sigcut*other.ini_slope_err))

    def x_extrapolate(self, other, rcut):
        xa, za = self.path[-1][0], self.path[-1][1]
        xb, zb = other.path[0][0], other.path[0][1]

        return ( math.fabs( xb - (xa+(zb-za)*self.end_slope)) < rcut) and (math.fabs( xa-(xb+(za-zb)*other.ini_slope)) < rcut)

    def z_extrapolate(self, other, rcut):
        xa, za = self.path[-1][0], self.path[-1][1]
        xb, zb = other.path[0][0], other.path[0][1]
        if(self.end_slope == 0 and other.ini_slope == 0) : 
            return True

        if(self.end_slope == 0):
            return (math.fabs(za - zb - (xa-xb)/other.ini_slope) < rcut)
        elif( other.ini_slope == 0):
            return ( math.fabs(zb - za - (xb-xa)/self.end_slope) < rcut)
        else:
            return ( math.fabs(zb - za - (xb-xa)/self.end_slope) < rcut) and (math.fabs(za - zb - (xa-xb)/other.ini_slope) < rcut)


    def joinable(self, other, dcut, sigcut, rcut):
        if(self.view != other.view): 
            return False
        if( self.dist(other) < dcut and self.slope_comp(other, sigcut) == True and self.x_extrapolate(other, rcut) and self.z_extrapolate(other, rcut)):            
            return True


    def merge(self, other):
        self.nHits += other.nHits
        self.chi2 += other.chi2 #should be refiltered though
        self.tot_charge += other.tot_charge
        self.len_path += other.len_path 
        self.len_path += self.dist(other)
        self.matched = -1

        if(self.path[0][1] > other.path[0][1]):
               self.ini_crp = self.ini_crp
               self.end_crp = other.end_crp

               self.ini_slope = self.ini_slope
               self.ini_slope_err = self.ini_slope_err
               self.end_slope = other.end_slope
               self.end_slope_err = other.end_slope_err
               
               self.path.extend(other.path)
               self.dQ.extend(other.dQ)
               self.len_straight = math.sqrt( pow(self.path[0][0]-self.path[-1][0], 2) + pow(self.path[0][1]-self.path[-1][1], 2) )

        else:
               self.ini_crp = other.ini_crp
               self.end_crp = self.end_crp

               self.ini_slope = other.ini_slope
               self.ini_slope_err = other.ini_slope_err
               self.end_slope = self.end_slope
               self.end_slope_err = self.end_slope_err

               self.path = other.path + self.path
               self.dQ = other.dQ + self.dQ
               self.len_straight = math.sqrt( pow(self.path[0][0]-self.path[-1][0], 2) + pow(self.path[0][1]-self.path[-1][1],2) )        

    def mini_dump(self):
        print("[", self.ini_crp, " -> ", self.end_crp,",",self.view,"] from (%.1f,%.1f)"%(self.path[0][0], self.path[0][1]), " to (%.1f, %.1f)"%(self.path[-1][0], self.path[-1][1]), " N = ", self.nHits, " L = %.1f/%.1f"%(self.len_straight, self.len_path), " Q = ", self.tot_charge )
               

class trk3D:
    def __init__(self, tv0, tv1):
        self.ini_crp = tv0.ini_crp
        self.end_crp = tv0.end_crp

        self.chi2    = 0.5*(tv0.chi2 + tv1.chi2)
        self.momentum = -1

        self.nHits_v0   = tv0.nHits
        self.nHits_v1   = tv1.nHits

        self.len_straight_v0 = tv0.len_straight
        self.len_straight_v1 = tv1.len_straight

        self.len_path_v0 = -1 #tv0.len_path
        self.len_path_v1 = -1 #tv1.len_path

        self.tot_charge_v0 = tv0.tot_charge
        self.tot_charge_v1 = tv1.tot_charge

        self.ini_theta = -1
        self.end_theta = -1
        self.ini_phi = -1
        self.end_phi = -1

        self.ini_x = tv0.path[0][0]
        self.ini_y = tv1.path[0][0]
        self.ini_z = 0.5*(tv0.path[0][1] + tv1.path[0][1])

        self.end_x = tv0.path[-1][0]
        self.end_y = tv1.path[-1][0]
        self.end_z = 0.5*(tv0.path[-1][1] + tv1.path[-1][1])

        self.t0 = 0.
        
        self.path_v0 = []
        self.path_v1 = []
        self.dQds_v0 = []
        self.dQds_v1 = []

    def set_view0(self, length, path, dqds):
        self.len_path_v0 = length
        self.path_v0     = path
        self.dQds_v0     = dqds
        
    def set_view1(self, length, path, dqds):
        self.len_path_v1 = length
        self.path_v1     = path
        self.dQds_v1     = dqds

    def matched(self, tv0, tv1):
        tv0.matched = evt_list[-1].nTracks3D
        tv1.matched = evt_list[-1].nTracks3D


    def angles(self, tv0, tv1):

        """ initial angles """
        slope_v0 = tv0.ini_slope #dx/dz
        slope_v1 = tv1.ini_slope #dy/dz
        self.ini_phi = math.degrees(math.atan2(slope_v1, slope_v0))
        self.ini_theta = math.degrees(math.atan2(math.sqrt(pow(slope_v0,2)+pow(slope_v1,2)),1))

        """ end angles """
        slope_v0 = tv0.end_slope #dx/dz
        slope_v1 = tv1.end_slope #dy/dz
        self.end_phi = math.degrees(math.atan2(slope_v1, slope_v0))
        self.end_theta = math.degrees(math.atan2(math.sqrt(pow(slope_v0,2)+pow(slope_v1,2)),1))

        
    def dump(self):
        print(" from (%.2f, %.2f, %.2f) to (%.2f, %.2f, %.2f)"%(self.ini_x, self.ini_y, self.ini_z, self.end_x, self.end_y, self.end_z))
        print(" %.2f ; %.2f"%(self.ini_theta, self.ini_phi), " -> %.2f ; %.2f "%( self.end_theta, self.end_phi), " L = %.2f / %.2f"%(self.len_path_v0, self.len_path_v1))


