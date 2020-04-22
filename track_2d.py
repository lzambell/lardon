import config as cf
import data_containers as dc
import math
import numpy as np
import sklearn.neighbors as skn
import scipy.stats as stat
import scipy.spatial as spatial

import pierre_filter as pf



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


    for icrp in range(cf.n_CRPUsed):
        for iview in range(cf.n_View):
            for icl in range(dc.evt_list[-1].nClusters[icrp,iview]):
                hits = [x for x in dc.hits_list if x.crp==icrp and x.view==iview and x.cluster==icl]
                
                """sort by decreasing Z and increasing channel """
                hits.sort()

                nHits = len(hits)

                visited = np.zeros((nHits),dtype=bool)
                X = np.asarray([[x.X,x.Z] for x in hits])

                

                """build the NN Tree"""
                tree = spatial.cKDTree(X)

                seeded = False
                
                while(np.sum(visited) < nHits):
                    
                    """get the first not yet visited hit in the list"""
                    idx = np.argmax(visited==0)
                    visited[idx] = True
                    
                    
                    """gets NN indices within rcut, 
                    first nearest point is itself"""
                    nn = tree.query_ball_point(X[idx], rcut,return_sorted=True)[1:]

                    """double rcut in case nothing found on first trial"""
                    if(len(nn)==0):
                        nn = tree.query_ball_point(X[idx], 2.*rcut, return_sorted=True)[1:]
                    
                    """give up if still nothing"""
                    if(len(nn)==0):
                        seeded = False
                        continue

                    """start the filter"""
                    if(seeded is False):
                        nn_idx = nn[np.argmax(visited[nn]==0)]

                        x0, x1 = hits[idx].Z, hits[nn_idx].Z
                        y0, y1 = hits[idx].X, hits[nn_idx].X
                        
                        if(x1 == x0):
                            seeded = False
                            continue

                        slope = (y1-y0)/(x1-x0)
                        intercept = y1 - slope * x1

                        ystart = slope * x0 + intercept
                        
                        filt.initiate(ystart, slope)
                        track = dc.trk2D(icrp,iview,slope, slope_err, y0, x0, hits[idx].charge, filt.getChi2())
                        track.add_hit(slope, filt.getSlopeErr(), y1, x1, hits[nn_idx].charge, filt.getChi2())
                        seeded = True
                        visited[nn_idx] = True

                    """update the track from nearby hits"""                   
                    finished = False
                        
                    while(seeded is True and finished is False and np.sum(visited) < nHits):
                        idx = nn_idx
                        x0, y0 = x1, y1
                        nn = tree.query_ball_point(X[idx], rcut, return_sorted=True)[1:]

                        if(len(nn)==0):
                            finished = True
                            seeded = False
                            if(track.nHits >= min_hits):
                                dc.tracks2D_list.append(track)
                                dc.evt_list[-1].nTracks2D[iview] += 1
                            continue

                        if(np.sum(visited[nn]) == len(nn)):
                            finished=True
                            seeded=False
                            if(track.nHits >= min_hits):
                                dc.tracks2D_list.append(track)
                                dc.evt_list[-1].nTracks2D[iview] += 1
                            continue
                            
                        
                        updated = False
                        best_idx = -1
                        best_chi2 = 99999.
                        
                        for j in nn :
                            first_not_vis = np.argmax(visited[nn]==0)

                            if(first_not_vis == 0 and visited[nn[0]] is False):
                                continue
                            nn_idx = nn[first_not_vis]
                            
                            x1, y1 = hits[nn_idx].Z, hits[nn_idx].X

                            if(x1 >= x0):
                                continue
                        
                            yp = filt.predict(x1-x0)
                            chi2m = filt.computeChi2(y1, x1-x0)
    
                        
                            if(chi2m < chicut):
                                if(chi2m < best_chi2):
                                    best_idx = nn_idx
                                    best_chi2 = chi2m

                        if(best_idx >= 0):
                            nn_idx = best_idx

                            x1, y1 = hits[nn_idx].Z, hits[nn_idx].X
                            chi2_up = filt.update(y1, x1-x0)
                            tot_chi2 = filt.getChi2()

                            track.add_hit(filt.getSlope(), filt.getSlopeErr(), y1, x1, hits[nn_idx].charge, tot_chi2)
                            visited[nn_idx]=True
                            updated = True
                            

                        if(updated is False or np.sum(visited) == nHits):
                            finished = True
                            seeded = False
                            if(track.nHits >=  min_hits):
                                dc.tracks2D_list.append(track)
                                dc.evt_list[-1].nTracks2D[iview] += 1
                            continue



    print("nb of 2D tracks ", len(dc.tracks2D_list))
    return
    

def stitch_tracks(dist_min, slope_err, r_extrapol_min):

    dc.tracks2D_list.sort()

    i = 0

    while(i < len(dc.tracks2D_list)):
        ti = dc.tracks2D_list[i]
        
        j = i+1
        while( j < len(dc.tracks2D_list) ):
            tj = dc.tracks2D_list[j]

            if(ti.joinable(tj, dist_min, slope_err, r_extrapol_min)):
                ti.merge(tj)
                del dc.tracks2D_list[j]

                dc.evt_list[-1].nTracks2D[ti.view] -= 1
                j = i+1
            else:
                j += 1
        i += 1
