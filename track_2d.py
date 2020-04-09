import config as cf
import math
import numpy as np
import sklearn.neighbors as skn
import scipy.stats as stat
import scipy.spatial as spatial

import pierre_filter as pf

def dump_track(idx, minHits):
    t = cf.tracks2D_list[idx]
    """
    if(t.nHits < minHits):
        return
    """
    print("Track ", idx)
    print("CRP : ", t.ini_crp, " to ", t.end_crp)
    print("VIEW : ", t.view)
    print("NB of hits attached ", t.nHits)
    print("X : ", t.path[0][0], " to ", t.path[0][1])
    print("Z : ", t.path[-1][0], " to ", t.path[-1][1])
    print("slope : %.2f [%.2f] to %.2f [%.2f]"%(t.ini_slope, t.ini_slope_err, t.end_slope, t.end_slope_err))
    if(t.nHits <= 3):
        print(t.path)
    print("Final Chi2 %.2f"%(t.chi2))
    print(" ")
    

def FindTracks(ncl, rcut, chicut, y_err, slope_err, pbeta):
    """error on y axis, error on slope, pbeta hyp"""
    #filt = pf.PFilter(0.15, 0.05, 3.)

    """error on y axis, error on slope, pbeta hyp"""
    filt = pf.PFilter(y_err, slope_err, pbeta)


    for icrp in range(2):
        for iview in range(2):
            for icl in range(ncl[icrp,iview]):
                #if(icl >= 1):
                    #return
                hits = [x for x in cf.hits_list if x.crp==icrp and x.view==iview and x.cluster==icl]
                
                """
                print("before sorting")
                for ii in range(len(hits)):
                    print(ii, " ", hits[ii].Z, " ", hits[ii].X)
                """
                #"""sort by increasing channel and increasing time"""
                """sort by decreasing Z and increasing channel """
                hits.sort()
                """
                print("\nafter sorting")
                for ii in range(len(hits)):
                    print(ii, " ", hits[ii].Z, " ", hits[ii].X)
                """

                nHits = len(hits)

                visited = np.zeros((nHits),dtype=bool)
                #X = np.asarray([[x.X,x.Z] for x in hits])
                X = np.asarray([[x.X,x.Z] for x in hits])

                
                #print("*-*-*-*-*-*-")
                #print("*-*-*-*-*-*-")
                #print("CRP ", icrp, " VIEW ", iview, " CLUSTER ",icl, " Npts ", nHits, " current nb of tracks ", len(cf.tracks2D_list))
                #print(" seeder : ",X[0])
                

                """build the NN Tree"""
                tree = spatial.cKDTree(X)

                seeded = False
                
                while(np.sum(visited) < nHits):
                    
                    """get the first not yet visited hit in the list"""
                    idx = np.argmax(visited==0)
                    visited[idx] = True
                    
                    #print(" while loop at idx : ", idx, " n visited ", np.sum(visited))
                    

                    """gets NN indices within rcut, 
                    first nearest point is itself"""
                    nn = tree.query_ball_point(X[idx], rcut,return_sorted=True)[1:]

                    """double rcut in case nothing found on first trial"""
                    if(len(nn)==0):
                        nn = tree.query_ball_point(X[idx], 2.*rcut, return_sorted=True)[1:]
                        #print("....trying with twice rcut")
                    
                    """give up if still nothing"""
                    if(len(nn)==0):
                        #print("still nothing ! ")
                        if(seeded is True):
                            #if(track.nHits == 2):
                                #visited[
                            #cf.tracks2D_list.append(track)
                            #track.reset()
                            seeded = False
                        continue

                    #print(len(nn), " NN : ", nn)

                    """start the filter"""
                    if(seeded is False):
                        nn_idx = nn[np.argmax(visited[nn]==0)]
                        #print("seeding with idx ", nn_idx)
                        
                        #x0, x1 = hits[idx].X, hits[nn_idx].X
                        x0, x1 = hits[idx].Z, hits[nn_idx].Z
                        #y0, y1 = hits[idx].Z, hits[nn_idx].Z
                        y0, y1 = hits[idx].X, hits[nn_idx].X
                        
                        if(x1 == x0):
                            #print("vertical track ... not handled atm")
                            seeded = False
                            continue
                        #print(" x0,y0 : ", x0, ", ", y0)
                        #print(" x1,y1 : ", x1, ", ", y1)

                        #!!!!!!!!!!!!!!!!!!!!!!!!!
                        slope = (y1-y0)/(x1-x0)
                        intercept = y1 - slope * x1

                        ystart = slope * x0 + intercept
                        #print(" ystart : ", ystart, " x ", x0, " slope ", slope, " intercept ", intercept)
                        
                        filt.initiate(ystart, slope)
                        #track = cf.trk2D(icrp,iview,slope, 0.05, x0, y0, hits[idx].charge, filt.getChi2())
                        track = cf.trk2D(icrp,iview,slope, slope_err, y0, x0, hits[idx].charge, filt.getChi2())
                        track.add_hit(slope, filt.getSlopeErr(), y1, x1, hits[nn_idx].charge, filt.getChi2())
                        seeded = True
                        visited[nn_idx] = True
                        #print("seeding finished! idx at ", idx, " nn_idx at ", nn_idx)
                        #continue

                    """update the track from nearby hits"""                   
                    finished = False
                        
                    while(seeded is True and finished is False and np.sum(visited) < nHits):
                        #print(" -> for loop at ",j)
                        idx = nn_idx
                        x0, y0 = x1, y1
                        nn = tree.query_ball_point(X[idx], rcut, return_sorted=True)[1:]
                        #print(" update with ", idx, " nb nn : ", len(nn), " n vis : ", np.sum(visited[nn]))

                        if(len(nn)==0):
                            #print(" ... cannot find any new neighbours")
                            finished = True
                            seeded = False
                            if(track.nHits > 3):
                                #visited[idx] = False
                                #track.reset()
                            #else:
                                cf.tracks2D_list.append(track)
                                #track.reset()
                            continue

                        if(np.sum(visited[nn]) == len(nn)):
                            #print(".... all NN have been visited")
                            finished=True
                            seeded=False
                            if(track.nHits > 3):
                                #visited[idx] = False
                                #track.reset()
                            #else:
                                cf.tracks2D_list.append(track)
                                #track.reset()
                            continue
                            
                        
                        updated = False
                        best_idx = -1
                        best_chi2 = 99999.
                        
                        for j in nn :
                            first_not_vis = np.argmax(visited[nn]==0)
                            #print("for loop at ", first_not_vis, " -> ", visited[first_not_vis])
                            if(first_not_vis == 0 and visited[nn[0]] is False):
                                continue
                            nn_idx = nn[first_not_vis]
                            
                            #x1, y1 = hits[nn_idx].X, hits[nn_idx].Z
                            x1, y1 = hits[nn_idx].Z, hits[nn_idx].X
                            #if(x1 <= x0):
                            if(x1 >= x0):
                                continue
                            #print(" ... trying with ", nn_idx)
                            #print(" x1, y1 : ", x1, ", ", y1)
                        
                            yp = filt.predict(x1-x0)
                            chi2m = filt.computeChi2(y1, x1-x0)
                        
                            #print("--> PREDICTION yp ", yp, " chi2 :", chi2m)
                        
                            if(chi2m < chicut):
                                if(chi2m < best_chi2):
                                    best_idx = nn_idx
                                    best_chi2 = chi2m
                        if(best_idx >= 0):
                            nn_idx = best_idx
                            #x1, y1 = hits[nn_idx].X, hits[nn_idx].Z
                            x1, y1 = hits[nn_idx].Z, hits[nn_idx].X
                            chi2_up = filt.update(y1, x1-x0)
                            tot_chi2 = filt.getChi2()

                            track.add_hit(filt.getSlope(), filt.getSlopeErr(), y1, x1, hits[nn_idx].charge, tot_chi2)
                            visited[nn_idx]=True
                            updated = True
                            #print(" OK NEXT! :: tot Chi2 : ", np.sum(visited),"/",nHits)
                            #break
                            

                        if(updated is False or np.sum(visited) == nHits):
                            #print(" ... all hits were visited !")
                            finished = True
                            seeded = False
                            if(track.nHits > 3):
                                #visited[idx] = False
                                #track.reset()
                            #else:
                                cf.tracks2D_list.append(track)
                                #track.reset()
                            continue



    print("nb tracks ", len(cf.tracks2D_list))
    for t in range(len(cf.tracks2D_list)):
        dump_track(t, 10)
    return
    
