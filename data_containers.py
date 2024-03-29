import config as cf
import numpy as np
import math
import time

map_ped = []
evt_list = []
hits_list = []
tracks2D_list = []
tracks3D_list = []


data = np.zeros((cf.n_CRPUsed, cf.n_View, cf.n_ChanPerCRP, cf.n_Sample), dtype=np.float32) #crp, view, vchan


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


ped_rms = np.zeros((cf.n_CRPUsed, cf.n_View, cf.n_ChanPerCRP), dtype=np.float32) 
ped_mean = np.zeros((cf.n_CRPUsed, cf.n_View, cf.n_ChanPerCRP), dtype=np.float32) 


CRP_gap = np.zeros( (cf.n_View, cf.n_CRP, cf.n_CRP), dtype=np.float32)

CRP_gap[0, 0, 1] = cf.crp_01_x/2.
CRP_gap[0, 1, 0] = -cf.crp_01_x/2.
CRP_gap[1, 0, 1] = cf.crp_01_y/2.
CRP_gap[1, 1, 0] = -cf.crp_01_y/2.

CRP_gap[0, 0, 3] = cf.crp_03_x/2.
CRP_gap[0, 3, 0] = -cf.crp_03_x/2.
CRP_gap[1, 0, 3] = cf.crp_03_y/2.
CRP_gap[1, 3, 0] = -cf.crp_03_y/2.


def reset_event():
    data[:,:,:,:] = 0.
    mask[:,:,:,:] = True
    ped_rms[:,:,:] = 0.
    ped_mean[:,:,:] = 0.

    hits_list.clear()
    tracks2D_list.clear()
    tracks3D_list.clear()
    evt_list.clear()

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

    def dump(self, verb):
        if(verb is True):

            print("RUN ",self.run_nb, " EVENT ", self.evt_nb_loc, " / ", self.evt_nb_glob,)
            print("Taken at ", time.ctime(self.time_s), " + ", self.time_ns, " ns ")
            print(" TS = ", self.time_s, " +", self.time_ns)

    def dump_reco(self, verb):
        if(verb is True):

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
        self.idx     = -1
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
        self.Z_start = -1
        self.Z_stop  = -1
        self.matched = -9999
        

    def __lt__(self,other):
        #"""sort hits by increasing channel and increasing Z"""
        #return (self.X < other.X) or (self.X== other.X and self.Z < other.Z)

        """ sort hits by decreasing Z and increasing channel """
        return (self.Z > other.Z) or (self.Z == other.Z and self.X < other.X)


    def set_index(self, idx):
        self.idx = idx


    def hit_positions(self, v):
        self.X = self.channel*cf.ChanPitch

        if(self.view == 0):
            if(self.crp == 1 or self.crp == 2):
                self.X -= cf.n_ChanPerCRP * cf.ChanPitch
        
        else:
            if(self.crp == 2 or self.crp == 3):
                self.X -= cf.n_ChanPerCRP * cf.ChanPitch                

        self.Z = cf.Anode_Z - self.max_t*cf.n_Sampling*v*0.1
        self.Z_start = cf.Anode_Z - self.start*cf.n_Sampling*v*0.1
        self.Z_stop = cf.Anode_Z - self.stop*cf.n_Sampling*v*0.1

    def hit_charge(self):
        self.charge *= cf.n_Sampling 
        self.charge /= cf.ADCtofC

    def set_match(self, ID):
        self.matched = ID

        
    def dump(self):
        print(self.crp, " ", self.view, " ", self.X, " [", self.Z_start, ", ", self.Z, ", ", self.Z_stop, "]")


class trk2D:
    def __init__(self, ID, ini_crp, view, ini_slope, ini_slope_err, x0, y0, t0, q0, chi2, cluster):
        self.trackID = ID
        self.ini_crp = ini_crp
        self.end_crp = ini_crp
        self.view    = view
    
        self.ini_slope       = ini_slope
        self.ini_slope_err   = ini_slope_err
        self.end_slope       = ini_slope
        self.end_slope_err   = ini_slope_err

        self.nHits      = 1
        self.nHits_dray = 0

        self.path    = [(x0,y0)]
        self.dQ      = [q0]

        self.chi2_fwd    = chi2
        self.chi2_bkwd   = chi2

        self.drays   = []
        
        self.tot_charge = q0
        self.dray_charge = 0.

        self.len_straight = 0.
        self.len_path = 0.

        self.matched = -1
        self.cluster = cluster

        self.ini_time = t0
        self.end_time = t0
        
    def __lt__(self,other):
        """ sort tracks by decreasing Z and increasing channel """
        return (self.path[0][1] > other.path[0][1]) or (self.path[0][1] == other.path[0][1] and self.path[0][0] < other.path[0][0])


    def add_drays(self, x, y, q):
        self.drays.append((x,y,q))
        self.dray_charge += q
        self.nHits_dray += 1
        self.remove_hit(x, y, q)


    def remove_hit(self, x, y, q):
        pos = -1
        for p,t in enumerate(self.path):
            if(t[0] == x and t[1] == y and self.dQ[p]==q):
                pos = p
                break

        if(pos >= 0):
            self.path.pop(pos)
            self.dQ.pop(pos)
            self.nHits -= 1
            self.tot_charge -= q
        else:
            print("?! cannot remove hit ", x, " ", y, " ", q, " pos ", pos)
            
    def add_hit(self, x, y, q, t):
        self.nHits += 1
        
        self.len_path += math.sqrt( pow(self.path[-1][0]-x, 2) + pow(self.path[-1][1]-y,2) )
        #beware to append (x,y) after !
        self.path.append((x,y))
        self.dQ.append(q)
        self.tot_charge += q
        self.len_straight = math.sqrt( pow(self.path[0][0]-self.path[-1][0], 2) + pow(self.path[0][1]-self.path[-1][1], 2) )
        self.end_time = t

    def add_hit_update(self, slope, slope_err, x, y, t, q, chi2):
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
        self.end_time = t

    def update_forward(self, chi2, slope, slope_err):
        self.chi2_fwd = chi2
        self.end_slope = slope
        self.end_slope_err = slope_err

    def update_backward(self, chi2, slope, slope_err):
        self.chi2_bkwd = chi2
        self.ini_slope = slope
        self.ini_slope_err = slope_err

    def reset_path(self, path, dQ):
        self.path = path
        self.dQ = dQ
        self.finalize_track()
        

    def finalize_track(self):
        if(self.path[-1][1] > self.path[0][1]):

            self.path.reverse()
            self.dQ.reverse()
            self.ini_crp, self.end_crp = self.end_crp, self.ini_crp
            self.ini_slope, self.end_slope = self.end_slope, self.ini_slope
            self.ini_slope_err, self.end_slope_err = self.end_slope_err, self.ini_slope_err

            self.chi2_fwd, self.chi2_bkwd = self.chi2_bkwd, self.chi2_fwd 
            print(self.trackID, " : wrong order check :", self.path[0][1], " to ", self.path[-1][1])
            self.ini_time, self.end_time = self.end_time, self.ini_time

        self.nHits = len(self.path)
        self.tot_charge = sum(self.dQ)

        self.nHits_dray = len(self.drays)
        self.dray_charge = sum(k for i,j,k in self.drays)

        self.len_straight = math.sqrt( pow(self.path[0][0]-self.path[-1][0], 2) + pow(self.path[0][1]-self.path[-1][1], 2) )        
        self.len_path = 0.
        for i in range(self.nHits-1):
            self.len_path +=  math.sqrt( pow(self.path[i][0]-self.path[i+1][0], 2) + pow(self.path[i][1]-self.path[i+1][1],2) )
            
            

    def dist(self, other, i=-1, j=0):
        return math.sqrt(pow( self.path[i][0] - other.path[j][0], 2) + pow(self.path[i][1] - other.path[j][1], 2))



    def slope_comp(self, other):#, sigcut):
        """ check if both tracks have the same slope direction """
        if(self.end_slope * other.ini_slope < 0.):
            return 9999. #False

        """ if slope error is too low, re-assign it to 5 percent """
        if(self.end_slope_err == 0 or math.fabs(self.end_slope_err/self.end_slope) < 0.05):
            end_err = math.fabs(self.end_slope*0.05)
        else:
            end_err = self.end_slope_err

        if(other.ini_slope_err == 0 or math.fabs(other.ini_slope_err/other.ini_slope) < 0.05):
            ini_err = math.fabs(other.ini_slope*0.05)
        else:
            ini_err = other.ini_slope_err

        #return (math.fabs( self.end_slope - other.ini_slope) < (sigcut*end_err + sigcut*ini_err))

        return math.fabs( self.end_slope - other.ini_slope) / (end_err + ini_err)


    def x_extrapolate(self, other, rcut):
        view = self.view
        a_crp = self.end_crp
        b_crp = other.ini_crp

        xa = self.path[-1][0] + CRP_gap[view,a_crp,b_crp]
        za = self.path[-1][1]
        xb = other.path[0][0] + CRP_gap[view,b_crp,a_crp]
        zb = other.path[0][1]
        return ( math.fabs( xb - (xa+(zb-za)*self.end_slope)) < rcut) and (math.fabs( xa-(xb+(za-zb)*other.ini_slope)) < rcut)

    def z_extrapolate(self, other, rcut):
        view = self.view
        a_crp = self.end_crp
        b_crp = other.ini_crp

        xa = self.path[-1][0] + CRP_gap[view,a_crp,b_crp]
        za = self.path[-1][1]
        xb = other.path[0][0] + CRP_gap[view,b_crp,a_crp]
        zb = other.path[0][1]

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
        if( self.dist(other) < dcut and self.slope_comp(other) <  sigcut and self.x_extrapolate(other, rcut) and self.z_extrapolate(other, rcut)):            
            return True


    def merge(self, other):
        self.nHits += other.nHits
        self.nHits_dray += other.nHits_dray
        self.chi2 += other.chi2 #should be refiltered though
        self.tot_charge += other.tot_charge
        self.dray_charge += other.dray_charge
        self.len_path += other.len_path 
        self.len_path += self.dist(other)
        self.matched = -1
        self.drays.extend(other.drays)

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
               self.end_time = other.end_time

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
               self.ini_time = other.ini_time

    def charge_in_z_interval(self, start, stop):
        return sum([q for q, (x, z) in zip(self.dQ, self.path) if z >= start and z <= stop])

    def mini_dump(self, verb=True):
        if(verb is True):
            if(True):#self.ini_crp != self.end_crp):            
                print("view : ", self.view, " [", self.ini_crp, " -> ", self.end_crp,",",self.view,"] from (%.1f,%.1f)"%(self.path[0][0], self.path[0][1]), " to (%.1f, %.1f)"%(self.path[-1][0], self.path[-1][1]), " N = ", self.nHits, " L = %.1f/%.1f"%(self.len_straight, self.len_path), " Q = ", self.tot_charge, " Dray N = ", self.nHits_dray, " Qdray ", self.dray_charge)
               

class trk3D:
    def __init__(self, tv0, tv1):
        self.ini_crp = tv0.ini_crp
        self.end_crp = tv0.end_crp

        self.chi2    = 0.5*(tv0.chi2 + tv1.chi2)
        self.momentum = -1

        self.nHits_v0   = tv0.nHits
        self.nHits_v1   = tv1.nHits

        self.len_straight_v0 = -1
        self.len_straight_v1 = -1

        self.len_path_v0 = -1 
        self.len_path_v1 = -1 

        self.tot_charge_v0 = -1
        self.tot_charge_v1 = -1
        self.dray_charge_v0 = -1
        self.dray_charge_v1 = -1

        self.ini_theta = -1
        self.end_theta = -1
        self.ini_phi = -1
        self.end_phi = -1

        self.t0_corr = 0.
        self.z0_corr = 0.

        self.ini_time = min(tv0.ini_time, tv1.ini_time)
        self.end_time = max(tv0.end_time, tv1.end_time)

        """ field corrected """
        self.len_straight_field_corr_v0 = -1
        self.len_straight_field_corr_v1 = -1

        self.len_path_field_corr_v0 = -1 
        self.len_path_field_corr_v1 = -1 

        self.tot_charge_field_corr_v0 = -1
        self.tot_charge_field_corr_v1 = -1

        self.ini_theta_field_corr = -1
        self.end_theta_field_corr = -1
        self.ini_phi_field_corr = -1
        self.end_phi_field_corr = -1


        self.t0_field_corr = 0.
        self.z0_field_corr = 0.

        
        ''' track boundaries '''                
        '''
        self.ini_x = tv0.path[0][0]
        self.ini_y = tv1.path[0][0]
        self.ini_z = 0.5*(tv0.path[0][1] + tv1.path[0][1])

        self.end_x = tv0.path[-1][0]
        self.end_y = tv1.path[-1][0]
        self.end_z = 0.5*(tv0.path[-1][1] + tv1.path[-1][1])
        '''

        self.path_v0 = []
        self.path_v1 = []
        self.dQds_v0 = []
        self.dQds_v1 = []

        self.z_field_corr_v0 = []
        self.z_field_corr_v1 = []
        self.dQds_field_corr_v0 = []
        self.dQds_field_corr_v1 = []

    def set_view0(self, path, dqds):
        self.len_straight_v0 = math.sqrt( pow(path[0][0]-path[-1][0], 2) + pow(path[0][1]-path[-1][1],2) + pow(path[0][2]-path[-1][2],2) )     
        self.len_path_v0 = 0.

        for i in range(len(path)-1):
            self.len_path_v0 +=  math.sqrt( pow(path[i][0]-path[i+1][0], 2) + pow(path[i][1]-path[i+1][1],2)+ pow(path[i][2]-path[i+1][2],2) )
            

        self.path_v0     = path
        self.dQds_v0     = dqds
        self.tot_charge_v0 = sum(q for q,s in dqds)

    def set_view1(self, path, dqds):
        self.len_straight_v1 = math.sqrt( pow(path[0][0]-path[-1][0], 2) + pow(path[0][1]-path[-1][1],2) + pow(path[0][2]-path[-1][2],2) )     

        self.len_path_v1 = 0.
        for i in range(len(path)-1):
            self.len_path_v1 +=  math.sqrt( pow(path[i][0]-path[i+1][0], 2) + pow(path[i][1]-path[i+1][1],2)+ pow(path[i][2]-path[i+1][2],2) )

        self.path_v1     = path
        self.dQds_v1     = dqds
        self.tot_charge_v1 = sum(q for q,s in dqds)

    def boundaries(self):
        ''' begining '''
        self.ini_x = self.path_v0[0][0] if self.path_v0[0][2] > self.path_v1[0][2] else self.path_v1[0][0]
        self.ini_y = self.path_v0[0][1] if self.path_v0[0][2] > self.path_v1[0][2] else self.path_v1[0][1]
        self.ini_z = self.path_v0[0][2] if self.path_v0[0][2] > self.path_v1[0][2] else self.path_v1[0][2]
        self.ini_z_overlap = min(self.path_v0[0][2], self.path_v1[0][2])

        ''' end '''
        self.end_x = self.path_v0[-1][0] if self.path_v0[-1][2] < self.path_v1[-1][2] else self.path_v1[-1][0]
        self.end_y = self.path_v0[-1][1] if self.path_v0[-1][2] < self.path_v1[-1][2] else self.path_v1[-1][1]
        self.end_z = self.path_v0[-1][2] if self.path_v0[-1][2] < self.path_v1[-1][2] else self.path_v1[-1][2]
        self.end_z_overlap = max(self.path_v0[-1][2], self.path_v1[-1][2])
        


    def matched(self, tv0, tv1):
        tv0.matched = evt_list[-1].nTracks3D
        tv1.matched = evt_list[-1].nTracks3D

    def set_t0_z0_corr(self, t0, z0):
        self.t0_corr = t0
        self.z0_corr = z0

    def set_field_correction(self, t0, z0, z_path_v0, z_path_v1, dqds_v0, dqds_v1):
        self.t0_field_corr = t0
        self.z0_field_corr = z0

        self.z_field_corr_v0 = z_path_v0
        self.z_field_corr_v1 = z_path_v1

        self.ini_field_corr_z = z_path_v0[0] if z_path_v0[0] > z_path_v1[0] else z_path_v1[0]
        self.end_field_corr_z = z_path_v0[-1] if z_path_v0[-1] < z_path_v1[-1] else z_path_v1[-1]

        self.ini_field_corr_z_overlap = min(z_path_v0[0], z_path_v1[0])
        self.end_field_corr_z_overlap = max(z_path_v0[-1], z_path_v1[-1])

        self.dQds_field_corr_v0 = dqds_v0
        self.dQds_field_corr_v1 = dqds_v1

        self.tot_charge_field_corr_v0 = sum(q for q,s in dqds_v0)
        self.tot_charge_field_corr_v1 = sum(q for q,s in dqds_v1)

        self.len_straight_field_corr_v0 = math.sqrt( pow(self.path_v0[0][0]-self.path_v0[-1][0], 2) + pow(self.path_v0[0][1]-self.path_v0[-1][1],2) + pow(z_path_v0[0]-z_path_v0[-1],2) )     

        self.len_straight_field_corr_v1 = math.sqrt( pow(self.path_v1[0][0]-self.path_v1[-1][0], 2) + pow(self.path_v1[0][1]-self.path_v1[-1][1],2) + pow(z_path_v1[0]-z_path_v1[-1],2) )     

        self.len_path_field_corr_v0 = 0.
        for i in range(self.nHits_v0-1):
            self.len_path_field_corr_v0 +=  math.sqrt( pow(self.path_v0[i][0]-self.path_v0[i+1][0], 2) + pow(self.path_v0[i][1]-self.path_v0[i+1][1],2)+ pow(z_path_v0[i]-z_path_v0[i+1],2) )

        self.len_path_field_corr_v1 = 0.
        for i in range(self.nHits_v1-1):
            self.len_path_field_corr_v1 +=  math.sqrt( pow(self.path_v1[i][0]-self.path_v1[i+1][0], 2) + pow(self.path_v1[i][1]-self.path_v1[i+1][1],2)+ pow(z_path_v1[i]-z_path_v1[i+1],2) )

    
        
    def angles(self, tv0, tv1):

        """ initial angles """
        slope_v0 = tv0.ini_slope #dx/dz
        slope_v1 = tv1.ini_slope #dy/dz
        self.ini_phi = math.degrees(math.atan2(slope_v1, slope_v0))
        self.ini_theta = math.degrees(math.atan2(math.sqrt(pow(slope_v0,2)+pow(slope_v1,2)),-1.))

        """ end angles """
        slope_v0 = tv0.end_slope #dx/dz
        slope_v1 = tv1.end_slope #dy/dz
        self.end_phi = math.degrees(math.atan2(slope_v1, slope_v0))
        self.end_theta = math.degrees(math.atan2(math.sqrt(pow(slope_v0,2)+pow(slope_v1,2)),-1.))


    def angles_field_corr(self, slope_v0_ini, slope_v0_end, slope_v1_ini, slope_v1_end):
        self.ini_phi_field_corr = math.degrees(math.atan2(slope_v1_ini, slope_v0_ini))
        self.ini_theta_field_corr = math.degrees(math.atan2(math.sqrt(pow(slope_v0_ini,2)+pow(slope_v1_ini,2)),-1.))

        self.end_phi_field_corr = math.degrees(math.atan2(slope_v1_end, slope_v0_end))
        self.end_theta_field_corr = math.degrees(math.atan2(math.sqrt(pow(slope_v0_end,2)+pow(slope_v1_end,2)),-1.))
  
        
    def dump(self, verb=True):
        if(verb is True):
            print('----')
            print(" from (%.2f, %.2f, %.2f) to (%.2f, %.2f, %.2f)"%(self.ini_x, self.ini_y, self.ini_z, self.end_x, self.end_y, self.end_z))
            print('z-overlap ', self.ini_z_overlap, ' to ', self.end_z_overlap)
            print(" theta, phi: [ini] %.2f ; %.2f"%(self.ini_theta, self.ini_phi), " -> [end] %.2f ; %.2f "%( self.end_theta, self.end_phi), " L = (P) %.2f / %.2f ; (S) %.2f / %.2f"%(self.len_path_v0, self.len_path_v1, self.len_straight_v0, self.len_straight_v1))

            print(" corr : z: %.2f cm / t: %.2f mus"%(self.z0_corr, self.t0_corr))
            
            print(" charge V0: %.2f, V1: %.2f A= %.3f"%(self.tot_charge_v0, self.tot_charge_v1, (self.tot_charge_v0-self.tot_charge_v1)/(self.tot_charge_v0+self.tot_charge_v1)))
            print("field corrected ")
            print(' track z bound ', self.ini_field_corr_z, ' to ', self.end_field_corr_z)
            print(' track z overlap ', self.ini_field_corr_z_overlap, ' to ', self.end_field_corr_z_overlap)

            print(" theta, phi: [ini] %.2f ; %.2f"%(self.ini_theta_field_corr, self.ini_phi_field_corr), " -> [end] %.2f ; %.2f "%( self.end_theta_field_corr, self.end_phi_field_corr), " L = (P) %.2f / %.2f ; (S) %.2f / %.2f"%(self.len_path_field_corr_v0, self.len_path_field_corr_v1, self.len_straight_field_corr_v0, self.len_straight_field_corr_v1))
            print(" corr : z %.2f cm / t %.2f mus"%(self.z0_field_corr, self.t0_field_corr))
            print(" charge V0: %.2f, V1: %.2f A= %.3f"%(self.tot_charge_field_corr_v0, self.tot_charge_field_corr_v1, (self.tot_charge_field_corr_v0-self.tot_charge_field_corr_v1)/(self.tot_charge_field_corr_v0+self.tot_charge_field_corr_v1)))
