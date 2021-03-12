import config as cf
import data_containers as dc
import math
import numpy as np
from operator import itemgetter
#import sklearn.neighbors as skn
import scipy.stats as stat
import scipy.spatial as spatial
import scipy.sparse.csgraph as csgr



import pierre_filter as pf


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
                        
                        filt.initiate(ystart, slope)
                        track = dc.trk2D(trackID, icrp,iview,slope, slope_err, y0, x0, hits[idx].charge, filt.getChi2(), icl)
                        track.add_hit_update(slope, filt.getSlopeErr(), y1, x1, hits[nn_idx].charge, filt.getChi2())
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
                                x1, y1 = hits[nn_idx].Z, hits[nn_idx].X
                                idx_list.append(nn_idx)
                                
                                if(j != best_idx):
                                    track.add_hit(y1, x1, hits[nn_idx].charge)
                                else:
                                    chi2_up = filt.update(y1, x1-x0)
                                    tot_chi2 = filt.getChi2()
                                    track.add_hit_update(filt.getSlope(), filt.getSlopeErr(), y1, x1, hits[nn_idx].charge, tot_chi2)
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
    n_NN = 3


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




        
    """ Forward filter the track with w/o drays """
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



    """ Backward filter the track with w/o drays """
        
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

    for l in drays:
        h = hits[l]
        if(h.matched >= 0):
            track.add_drays(h.X, h.Z, h.charge)
            h.set_match(-1*h.matched)

    track.finalize_track()

    


def stitch_tracks(dist_min, slope_err_tol, r_extrapol_min, y_err, slope_err, pbeta):


    dc.tracks2D_list.sort()
    i = 0

    while(i < len(dc.tracks2D_list)):
        ti = dc.tracks2D_list[i]

        j = 0 
        join = []
        while( j < len(dc.tracks2D_list) ):
            if(i==j): 
                j += 1
                if(j >= len(dc.tracks2D_list) ):
                    break

            tj = dc.tracks2D_list[j]
            
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




