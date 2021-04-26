import config as cf
import data_containers as dc
import math
import numpy as np
from operator import itemgetter
#import sklearn.neighbors as skn
import scipy.stats as stat
import scipy.spatial as spatial
import scipy.sparse.csgraph as csgr

import R_tree as myrtree

from scipy.interpolate import UnivariateSpline

import pierre_filter as pf


def linear_reg(X, Y, cut):
    N = len(X)
    x0, y0 = X[0], Y[0] 
    fits = []
    for i in range(1, N):
        for j in range(i+1, N):
            x = [x0, X[i], X[j]]
            y = [y0, Y[i], Y[j]]
            slope, intercept, r, p, s = stat.linregress(x, y)
            if(r**2 > cut):
                fits.append( (r**2, slope, intercept, i, j) )

    if(len(fits)==0):
        return (9999., 0, 0, -1, -1)
    elif(len(fits)==1):
        return fits[0]
    else:
        """ sort by correlation factor """
        fits = sorted(fits, key=itemgetter(1))        
        return fits[0]

def get_path(Pr, i, j):
    path = [j]
    k = j
    while Pr[i, k] != -9999:
        path.append(Pr[i, k])
        k = Pr[i, k]
    return path[::-1]


def fit_slopes(n, t):
    #track begining
    pts = np.array([[p[1],p[0]] for p in  t.path[:n]])
    slope, interc, r, prob, see = stat.linregress(pts)
    mx = pts.mean(axis=0)[0]
    sx2 = ((pts[:,0]-mx)**2).sum()

    slope_err = see * np.sqrt( 1./sx2) if sx2 != 0 else 0.

    t.ini_slope = slope
    t.ini_slope_err = slope_err


    #track ending
    pts = np.array([[p[1],p[0]] for p in  t.path[-n:]])
    slope, interc, r, prob, see = stat.linregress(pts)
    mx = pts.mean(axis=0)[0]
    sx2 = ((pts[:,0]-mx)**2).sum()
    slope_err = see * np.sqrt( 1./sx2) if sx2 != 0 else 0.

    t.end_slope = slope
    t.end_slope_err = slope_err



def dump_track(idx):
    t = dc.tracks2D_list[idx]

    print("Track ", idx)
    print("CRP : ", t.ini_crp, " to ", t.end_crp)
    print("VIEW : ", t.view)
    print("NB of hits attached ", t.nHits)
    print("X : ", t.path[0][0], " to ", t.path[0][1])
    print("Z : ", t.path[-1][0], " to ", t.path[-1][1])
    print("slope : %.2f [%.2f] to %.2f [%.2f]"%(t.ini_slope, t.ini_slope_err, t.end_slope, t.end_slope_err))
    print("Final Chi2 %.2f"%(t.chi2))
    print(" ")
    

def find_tracks(min_hits, rcut, chicut, y_err, slope_err, pbeta):

    """error on y axis, error on slope, pbeta hyp"""
    filt = pf.PFilter(y_err, slope_err, pbeta)

    """track ID starts at 1 """
    trackID = len(dc.tracks2D_list)+1

    for icrp in range(cf.n_CRPUsed):
        for iview in range(cf.n_View):
            for icl in range(dc.evt_list[-1].nClusters[icrp,iview]):
                hits = [x for x in dc.hits_list if x.crp==icrp and x.view==iview and x.matched == -9999 and x.cluster==icl]
                
                """sort by decreasing Z and increasing channel """
                hits.sort()

                nHits = len(hits)
                if(nHits < min_hits):
                    continue
                visited = np.zeros((nHits),dtype=bool)
                X = np.asarray([[x.X,x.Z] for x in hits])                

                """build the NN Tree"""
                tree = spatial.cKDTree(X)

                seeded = False
                
                while(np.sum(visited) < nHits):
                    
                    idx_list = []

                    """get the first not yet visited hit in the list"""
                    idx = np.argmax(visited==0)
                    idx_list.append(idx)

                    """ in case everything has already been visited """
                    if(idx == 0 and visited[idx] is True):
                        break
                    
                    visited[idx] = True
                    
                    
                    """gets NN indices within rcut, 
                    first nearest point is itself"""
                    nn = tree.query_ball_point(X[idx], rcut, return_sorted=True)#[1:]

                    nn = [k for k in nn if (k != idx) and (visited[k] == 0)]



                    """double rcut in case nothing found on first trial"""
                    #if(len(nn)==0):
                        #nn = tree.query_ball_point(X[idx], 2.*rcut, return_sorted=True)#[1:]
 

                    #nn = [k for k in nn if (k != idx) and (visited[k] == 0)]
                   
                    """give up if still nothing"""
                    if(len(nn)==0):
                        seeded = False
                        continue

                    """start the filter"""
                    if(seeded is False):
                        nn_idx = nn[0]#np.argmax(visited[nn]==0)]

                        x0, x1 = hits[idx].Z, hits[nn_idx].Z
                        y0, y1 = hits[idx].X, hits[nn_idx].X
                        
                        if(x1 == x0):
                            seeded = False
                            continue

                        slope = (y1-y0)/(x1-x0)
                        intercept = y1 - slope * x1

                        ystart = slope * x0 + intercept

                        t0, t1 = hits[idx].max_t, hits[nn_idx].max_t
                        filt.initiate(ystart, slope)
                        track = dc.trk2D(trackID, icrp,iview,slope, slope_err, y0, x0, t0, hits[idx].charge, filt.getChi2(), icl)
                        track.add_hit_update(slope, filt.getSlopeErr(), y1, x1, t1, hits[nn_idx].charge, filt.getChi2())
                        idx_list.append(nn_idx)
                        
                        seeded = True
                        visited[nn_idx] = True

                    """update the track from nearby hits"""                   
                    finished = False
                    
                    while(seeded is True and finished is False and np.sum(visited) < nHits):
                        idx = nn_idx
                        x0, y0 = x1, y1
                        nn = tree.query_ball_point(X[idx], rcut, return_sorted=True)#[1:]
                        nn = [k for k in nn if (k != idx) and (visited[k] == 0)]                        


                        if(len(nn)==0):
                            finished = True
                            seeded = False
                            if(track.nHits >= min_hits):
                                dc.tracks2D_list.append(track)
                                [hits[i].set_match(trackID) for i in idx_list]
                                refilter_and_find_drays(trackID,
                                                        y_err, slope_err, pbeta)
                                dc.evt_list[-1].nTracks2D[iview] += 1
                                trackID += 1

                            continue

                        if(np.sum(visited[nn]) == len(nn)):
                            finished=True
                            seeded=False
                            if(track.nHits >= min_hits):
                                dc.tracks2D_list.append(track)
                                [hits[i].set_match(trackID) for i in idx_list]
                                refilter_and_find_drays(trackID,
                                                        y_err, slope_err, pbeta)
                                dc.evt_list[-1].nTracks2D[iview] += 1
                                trackID += 1
                            continue
                            
                        
                        updated = False

                        ok_idx = []
                        best_idx = -1
                        best_chi2 = 99999.
                        
                        for j in nn: 
                            nn_idx = j
                            
                            x1, y1 = hits[nn_idx].Z, hits[nn_idx].X

                            if(x1 >= x0):
                                continue
                        
                            yp = filt.predict(x1-x0)
                            chi2m = filt.computeChi2(y1, x1-x0)
    
                            if(chi2m < chicut):
                                ok_idx.append(nn_idx)
                                
                                if(chi2m < best_chi2):
                                    best_idx = nn_idx
                                    best_chi2 = chi2m

                        if(len(ok_idx) > 0):
                            for j in ok_idx:
                                nn_idx = j
                                x1, y1 ,t1 = hits[nn_idx].Z, hits[nn_idx].X, hits[nn_idx].max_t
                                idx_list.append(nn_idx)
                                
                                if(j != best_idx):
                                    track.add_hit(y1, x1, hits[nn_idx].charge)
                                else:
                                    chi2_up = filt.update(y1, x1-x0)
                                    tot_chi2 = filt.getChi2()
                                    track.add_hit_update(filt.getSlope(), filt.getSlopeErr(), y1, x1, t1, hits[nn_idx].charge, tot_chi2)
                                visited[nn_idx]=True
                                updated = True

                                if(j==best_idx):
                                    break

                        if(updated is False or np.sum(visited) == nHits):
                            finished = True
                            seeded = False
                            if(track.nHits >=  min_hits):
                                
                                dc.tracks2D_list.append(track)
                                [hits[i].set_match(trackID) for i in idx_list]
                                refilter_and_find_drays(trackID,
                                                        y_err, slope_err, pbeta)

                                dc.evt_list[-1].nTracks2D[iview] += 1
                                trackID += 1
                            continue

    return
    




def refilter_and_find_drays(idtrk, y_err, slope_err, pbeta):


    """error on y axis, error on slope, pbeta hyp"""
    filt = pf.PFilter(y_err, slope_err, pbeta)
    n_NN = 4


    track = [x for x in dc.tracks2D_list if x.trackID == idtrk]
    if(len(track) == 0 or len(track) > 1):
        print(" THERE IS AN ID PROBLEM !!")
    else:
        track = track[0]

    hits = [x for x in dc.hits_list if x.matched==idtrk]

    """sort by decreasing Z and increasing channel """
    hits.sort()

    """ sort by decreasing z and increasing x """
    coord = [(x.X, x.Z) for x in hits]
    charge = [x.charge for x in hits]    
    #print(coord)

    track.reset_path(coord, charge)
    
    coord = np.asarray(coord)
    """ compute the distance between each points"""
    graph = spatial.distance.cdist(coord, coord, 'euclidean')
    """ keep only the two closest points """
    graph = graph * (graph < np.sort(graph, axis=-1)[:,[n_NN]])
    """ keep only short edges """
    #graph[graph > rcut] = 0.
        
    """ compute the MST from this graph """
    T = csgr.minimum_spanning_tree(csgraph=graph)

    """ get the number of disconnected graphs """
    #n_components, labels = csgr.connected_components(csgraph=T, directed=False, return_labels=True)
    T = T.toarray()
    n_elem = np.count_nonzero(T, axis=0) + np.count_nonzero(T, axis=1)
    #solo = np.nonzero(n_elem==0)[0]
    borders = np.nonzero(n_elem==1)[0]
    vertex  = np.nonzero(n_elem>2)[0]




    """ identify potential delta rays from MST"""
    drays = []
    D, Pr = csgr.shortest_path(T, directed=False, method='FW', return_predecessors=True)
    for v in vertex:
        d_min = 9999.
        for b in borders:
            if(D[v, b] < d_min):
                p = get_path(Pr,v,b)
                d_min = D[v,b]
        for i in p[1:]:
            drays.append(i)

    #print("FOUND ", len(drays), " DELTA RAYS")
    for l in drays:
        h = hits[l]
        if(h.matched >= 0):
            track.add_drays(h.X, h.Z, h.charge)
            h.set_match(-1*h.matched)
            #print(h.X, " ", h.Z)

    drays.clear()

    #reversed because spline wants an increasing x only
    x_r = [x.X for x in reversed(hits) if x.matched > 0]
    z_r = [x.Z for x in reversed(hits) if x.matched > 0]

    """ spline needs unique 'x' points to work --> remove duplicate """
    z_r_unique, idx = np.unique(z_r, return_index=True)
    x_r = np.asarray(x_r)
    x_r_unique = x_r[idx]


    """at least 3 points for the spline """
    if(len(z_r_unique) < 4):
        print("PROBBLEM")
        track.update_forward(9999., 9999., 9999.)
        track.update_backward(9999., 9999., 9999.)

    else:
        spline = UnivariateSpline(z_r_unique, x_r_unique)

        deriv = spline.derivative()
        '''
        print("-----")
        print(" nb of hits unique ", len(x_r_unique), " vs ", len(x_r))
        print(" deriv at (top) z = %.2f = %.2f"%(z_r[-1],  deriv(z_r[-1])))
        print(" deriv at (bottom) z = %.2f = %.2f"%(z_r[0],  deriv(z_r[0])))
        print(" --> error ", spline.get_residual())
        '''

        for i, (ix, iz) in enumerate(zip(x_r, z_r)):
            if(np.fabs(spline(iz) - ix) > 1.):
                drays.append(i)

        
        #print("FOUND ", len(drays), " DELTA RAYS (far from spline)")
        for l in drays:
            h = hits[l]
            if(h.matched >= 0):
                track.add_drays(h.X, h.Z, h.charge)
                h.set_match(-1*h.matched)
                #print(h.X, " ", h.Z)



        track.update_forward(spline.get_residual(), deriv(z_r[0]), deriv(z_r[0])*0.05)
        track.update_backward(spline.get_residual(), deriv(z_r[-1]), deriv(z_r[-1])*0.05)

        

        
    """ Forward filter the track with w/o drays """
    """
    tot_fwd_chi2 = -1
    i=0
    while(i in drays):
        i+=1
    x0, y0 = coord[i][1], coord[i][0]
    x1 = x0


    while(x1 == x0):
        i+=1
        while(i in drays):
            i+=1
        x1, y1 = coord[i][1], coord[i][0]


    slope = (y1-y0)/(x1-x0)
    intercept = y1 - slope * x1
    ystart = slope * x0 + intercept
        
    filt.initiate(ystart, slope)

    x0, y0 = x1, y1
    maxchi = -1
    for label, ih in enumerate(coord[i+1:]):
        label += i+1
        if(label in drays):
            continue
        x1, y1 = ih[1], ih[0]
        chi2_up = filt.update(y1, x1-x0)
        tot_fwd_chi2 = filt.getChi2()
        x0, y0 = x1, y1

    track.update_forward(tot_fwd_chi2, filt.getSlope(), filt.getSlopeErr())
    print(" from fwd filter slope ", filt.getSlope(), " +/- ", filt.getSlopeErr())      
    """

    """ Backward filter the track with w/o drays """
    
    """
    tot_bkd_chi2 = -1
    i=len(coord)-1
    while(i in drays):
        i-=1
    x0, y0 = coord[i][1], coord[i][0]
    x1 = x0


    while(x1 == x0):
        i-=1
        while(i in drays):
            i-=1
        x1, y1 = coord[i][1], coord[i][0]


    slope = (y1-y0)/(x1-x0)
    intercept = y1 - slope * x1
    ystart = slope * x0 + intercept
        
    filt.initiate(ystart, slope)

    x0, y0 = x1, y1
    maxchi = -1
    for label, ih in enumerate(coord[i-1::-1]):
        label = i-1-label
        if(label in drays):
            continue
        x1, y1 = ih[1], ih[0]
        chi2_up = filt.update(y1, x1-x0)
        tot_bkd_chi2 = filt.getChi2()
        x0, y0 = x1, y1

    track.update_backward(tot_bkd_chi2, filt.getSlope(), filt.getSlopeErr())

    print(" from bkw filter slope ", filt.getSlope(), " +/- ", filt.getSlopeErr())      
    """
    track.finalize_track()

    


def stitch_tracks(dist_min, slope_err_tol, r_extrapol_min, y_err, slope_err, pbeta):


    dc.tracks2D_list.sort()

    i = 0

    while(i < len(dc.tracks2D_list)):
        ti = dc.tracks2D_list[i]
        if(ti.chi2_fwd == 9999. and ti.chi2_bkwd == 9999.):
            i += 1
            continue
        j = 0 
        join = []
        while( j < len(dc.tracks2D_list) ):
            if(i==j): 
                j += 1
                if(j >= len(dc.tracks2D_list) ):
                    break

            tj = dc.tracks2D_list[j]
            if(tj.chi2_fwd == 9999. and tj.chi2_bkwd == 9999.):
                j += 1
                continue
            
            if(ti.path[0][1] > tj.path[0][1]): 
                if(ti.joinable(tj, dist_min, slope_err_tol, r_extrapol_min)):
                    join.append( (tj, j, ti.slope_comp(tj)))
            j += 1

        if(len(join)==0):
            i += 1
            continue

        join = sorted(join, key=itemgetter(2))
        for tm, m, sm in join:

            """ check if there is no better match for tm """
            k = 0
            better_option = False
            while( k < len(dc.tracks2D_list) ):
                if(k == i or k==m): 
                    k += 1
                    if(k >= len(dc.tracks2D_list) ):
                        break

                tk = dc.tracks2D_list[k]
            
                if(tk.path[0][1] > tm.path[0][1]):
                    if(tk.joinable(tm, dist_min, slope_err_tol, r_extrapol_min)):
                        if(tk.slope_comp(tm) < sm):
                            better_option = True
                            break
                k += 1

            if(better_option == False):
                tmID = tm.trackID
                tiID = ti.trackID

                ti.merge(tm)
                hits = [h for h in dc.hits_list if math.fabs(h.matched)==tmID]
            
                for h in hits:
                    h.set_match( tiID if h.matched > 0 else -tiID )

                refilter_and_find_drays(tiID, y_err, slope_err, pbeta) 
                del dc.tracks2D_list[m]
                dc.evt_list[-1].nTracks2D[ti.view] -= 1
                i = 0
                break
        else:    
            i += 1








def find_tracks_rtree(min_hits, rcut, chicut, y_err, slope_err, pbeta):

    """error on y axis, error on slope, pbeta hyp"""
    filt = pf.PFilter(y_err, slope_err, pbeta)

    """track ID starts at 1 """
    trackID = len(dc.tracks2D_list)+1

    """ initialize the R-tree """
    tt = myrtree.R_tree(rcut)

    for icrp in range(cf.n_CRPUsed):
        for iview in range(cf.n_View):
            icl=0
            hits = [x for x in dc.hits_list if x.crp==icrp and x.view==iview and x.matched == -9999]

            nHits = len(hits)
            #print(nHits)
            if(nHits < min_hits):
                continue

            """sort by decreasing Z and increasing channel """
            hits.sort()            
            
            """ build the R-tree """
            tt.create_index()

            [tt.insert_hit(h, i) for i,h in enumerate(hits)]

            visited = np.zeros((nHits),dtype=bool)
            
            seeded = False  

            while(np.sum(visited) < nHits):
                #print(" new attempt, ", np.sum(visited), " hits left")
                idx_list = []

                """get the first not yet visited hit in the list"""
                idx = np.argmax(visited==0)
                idx_list.append(idx)
                #print(" First not yet visited ", idx)

                """ in case everything has already been visited """
                if(idx == 0 and visited[idx] is True):
                    break
                    
                visited[idx] = True
                tt.remove_hit(hits[idx], idx)

                """ get the N nearest hits """
                nn_id = tt.nearest_id(hits[idx], 5)                
                """ sort them by closest distance """
                nn = [(i,tt.distance(hits[idx], hits[i])) for i in nn_id if tt.close_enough(hits[idx], hits[i])]
                nn.sort(key=itemgetter(1))

                #hits[idx].dump()
                #print("NN: ", nn)

                """start the filter"""
                if(seeded is False and len(nn)>2):
                    
                    """ fit lines with all NN combination and idx """
                    X = [hits[idx].Z]
                    Y = [hits[idx].X]
                    X.extend([hits[i].Z for i,d in nn])
                    Y.extend([hits[i].X for i,d in nn])

                    """ returns the sorted fit results """
                    fit = linear_reg(X, Y, 0.9)
                    

                    #print("SEEDING with:")
                    #print([(x,z) for x,z in zip(X,Y)])
                    #print("BEST is ")
                    #print(fit)

                    x0, y0, t0 = hits[idx].Z, hits[idx].X, hits[idx].max_t

                    if(fit[0] != 9999):
                        seeded = True
                        slope  = fit[1]
                        intercept = fit[2]
                        ystart  = slope*x0 + intercept
                        filt.initiate(ystart, slope)
                        track = dc.trk2D(trackID, icrp, iview, slope, slope_err, y0, x0, t0, hits[idx].charge, filt.getChi2(), icl)
                        #print(" track seeded slope ", slope, " intcpt ", intercept, " ystart ", ystart)
                        
                        for i in [fit[3], fit[4]]:
                            """ add the seeding hits to the filter and remove them from the index """
                            nn_idx = nn[i-1][0]
                            x1, y1, t1 = hits[nn_idx].Z, hits[nn_idx].X, hits[nn_idx].max_t
                            chi2_up = filt.update(y1, x1-x0)
                            track.add_hit_update(slope, filt.getSlopeErr(), y1, x1, t1, hits[nn_idx].charge, filt.getChi2())
                            idx_list.append(nn_idx)
                    
                            visited[nn_idx] = True
                            tt.remove_hit(hits[nn_idx], nn_idx)
                            #print(" --> adding ", nn_idx, " chi2_up ", chi2_up)


                """update the track from nearby hits"""                   
                finished = False
                    
                while(seeded is True and finished is False and np.sum(visited) < nHits):
                    idx = nn_idx
                    x0, y0 = x1, y1
                    #print(" searching to extend track from ", idx, " at ", x0, ", ", y0)
                    #print("Filter slope ", filt.getSlope(), " ys ", filt.getY())

                    """ get the N nearest hits """
                    nn_id = tt.nearest_id(hits[idx], 15)                
                    #nn = [(i,tt.distance(hits[idx], hits[i])) for i in nn_id if tt.close_enough(hits[idx], hits[i])]
                    """ sort the NN hits by |pred-meas| distance if they are close enough """
                    nn = [(i,filt.delta_y(hits[i].X, hits[i].Z-x0)) for i in nn_id if tt.close_enough(hits[idx], hits[i])]
                    nn.sort(key=itemgetter(1))
                    #print(" NN: ", nn)

                    if(len(nn)==0):
                        finished = True
                        seeded = False
                        #print("no more close by hits")
                        if(track.nHits >= min_hits):
                            dc.tracks2D_list.append(track)
                            [hits[i].set_match(trackID) for i in idx_list]
                            #print(" storing track, Nidx ", len(idx_list), " chi tot ", filt.getChi2())
                            #dc.tracks2D_list[-1].mini_dump()
                            refilter_and_find_drays(trackID,
                                                    y_err, slope_err, pbeta)
                            #print("refilter and shit")
                            dc.tracks2D_list[-1].mini_dump()
                            dc.evt_list[-1].nTracks2D[iview] += 1
                            trackID += 1
                            #print("\n")
                        continue

                    updated = False

                    ok_idx = []
                    best_idx = -1
                    best_chi2 = 99999.
                        
                    for j,d in nn: 
                        nn_idx = j
                        if(d > 2.): continue

                        x1, y1 = hits[nn_idx].Z, hits[nn_idx].X

                        if(x1 >= x0):
                            if(tt.overlap_in_time(hits[idx], hits[nn_idx])==False):
                                continue
                        
                        yp = filt.predict(x1-x0)
                        chi2m = filt.chi2_if_update(y1, x1-x0)
                        slope = (y1-y0)/(x1-x0) if x1!=x0 else 0
                        prod_slope = slope * filt.getSlope()
                        #print("hit ", nn_idx, " at ", x1, ", ", y1, " d: ", d," chi2m ", chi2m, " slope ", slope, " -> ",slope*filt.getSlope())
                        if(chi2m < chicut and prod_slope >= 0.):
                            ok_idx.append( (nn_idx, tt.distance(hits[idx], hits[nn_idx])) )
                                
                            if(chi2m < best_chi2):
                                best_idx = nn_idx
                                best_chi2 = chi2m
                    #print(" --> best hit is ", best_idx, " at ", best_chi2)


                    """ sort accepted hits by distance to idx """ 
                    ok_idx.sort(key=itemgetter(1))
                    #print(ok_idx)

                    if(len(ok_idx) > 0):
                        for j,t in ok_idx:
                            
                            nn_idx = j
                            x1, y1, t1 = hits[nn_idx].Z, hits[nn_idx].X, hits[nn_idx].max_t
                            #print(" update track with ", nn_idx, " at ", x1, ", ", y1, " at a distance ", t)

                            d = filt.delta_y(y1, x1-x0)
                            #print("new pred distance : ", d)
                            if(d>2.):
                                continue

                            chi2m = filt.chi2_if_update(y1, x1-x0)
                            #print(" chi2 --> ", chi2m)
                            if(chi2m < chicut):

                                chi2_up = filt.update(y1, x1-x0)
                                track.add_hit_update(filt.getSlope(), filt.getSlopeErr(), y1, x1, t1, hits[nn_idx].charge, filt.getChi2())
                                idx_list.append(nn_idx)
                                #print("DONE ... slope ", filt.getSlope(), " ys ", filt.getY())
                                
                                visited[nn_idx]=True
                                tt.remove_hit(hits[nn_idx], nn_idx)
                                updated = True
                                x0, y0 = x1, y1
                            else:
                                #print("NOPE")
                                continue


                            if(j==best_idx):
                                break
                            '''
                            if(j != best_idx):
                                """ hits in between idx and best are just added to the track list """
                                #track.add_hit(y1, x1, hits[nn_idx].charge)
                                print(" just added to the track")
                            else:
                                """ best hit updates the filter """
                                chi2_up = filt.update(y1, x1-x0)

                                track.add_hit_update(filt.getSlope(), filt.getSlopeErr(), y1, x1, hits[nn_idx].charge, filt.getChi2())
                                print(" updates the filter ", chi2_up, " tot ", filt.getChi2())
                                break
                            '''
                    if(updated is False or np.sum(visited) == nHits):
                        finished = True
                        seeded = False
                        #print(" couldn't update the track further")
                        if(track.nHits >=  min_hits):
                            dc.tracks2D_list.append(track)
                            [hits[i].set_match(trackID) for i in idx_list]
                            #print(" storing track, Nidx ", len(idx_list), " chi tot ", filt.getChi2())
                            #dc.tracks2D_list[-1].mini_dump()
                            refilter_and_find_drays(trackID,
                                                    y_err, slope_err, pbeta)
                            #print(" refilter and shit")
                            dc.tracks2D_list[-1].mini_dump()
                            dc.evt_list[-1].nTracks2D[iview] += 1
                            trackID += 1
                            #print("\n")
                        continue

    return
