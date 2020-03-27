import config as cf

import numpy as np
import sklearn.cluster as skc



def rebin(data, chan_rebin, tdc_rebin):
    """get current shape"""
    crp, view, chan, tdc = data.shape

    """make sure rebinning is possible"""
    if( (chan%chan_rebin) > 0): 
        "rebinning along channels not possible: ", chan," is not dividable by ", chan_rebin
        return data

    if( (tdc%tdc_rebin) > 0): 
        "rebinning along time not possible: ", tdc, " is not dividable by ", tdc_rebin
        return data

    new_n_chan = int(chan/chan_rebin)
    new_n_tdc = int(tdc/tdc_rebin)
    
    data = data.reshape(2, cf.n_View, new_n_chan, chan_rebin, new_n_tdc, tdc_rebin)
    data = data.sum(axis=5).sum(axis=3)
    return data


def dbscan(ncl, eps, min_samp, y_squeez):    
    for icrp in range(2):
        for iview in range(2):
            """try to cluster un-clustered hits only"""
            hits = [x for x in cf.hits_list if x.crp==icrp and x.view==iview and x.cluster == -1]

            """ squeeze y axis instead of rebinning or defining a new metric """
            data = [[x.channel,x.max_t*y_squeez] for x in hits]
            X = np.asarray(data)
            db = skc.DBSCAN(eps=eps,min_samples=min_samp).fit(X)
            labels = db.labels_

            for h,cluster in zip(hits, labels):
                h.cluster = cluster if cluster==-1 else cluster + ncl[icrp,iview]

            n_clusters = len(set(labels)) - (1 if -1 in labels else 0)
            #n_noise = list(labels).count(-1)
            #unique_labels = set(labels)
            ncl[icrp,iview] += n_clusters #len(unique_labels)-1

            #print("CRP ", icrp, " View ", iview)
            #print("number of clusters : %d" % n_clusters)
            #print('Number of noise points: %d' % n_noise)
    return ncl
    
