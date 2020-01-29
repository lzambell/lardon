import config as cf

import numpy as np
import sklearn.cluster as skc
import itertools as itr

"""to be removed"""
import matplotlib as plt



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



    
def dbscan(data, eps, min_samp):
    #res = np.where(ROI[1,0,:,:]>0)
    X = np.asarray(data).T
    
    db = skc.DBSCAN(eps=eps, min_samples=min_samp).fit(X)
    labels = db.labels_
    n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
    n_noise_ = list(labels).count(-1)


    print("number of estimated clusters : %d" % n_clusters_)
    print('Estimated number of noise points: %d' % n_noise_)

    """
    X = X[labels>=0]
    db = cluster.DBSCAN(eps=15, min_samples=30).fit(X)    
    labels = db.labels_
    n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
    n_noise_ = list(labels).count(-1)
    print("number of estimated clusters : %d" % n_clusters_)
    print('Estimated number of noise points: %d' % n_noise_)
    """

    unique_labels = set(labels)

    """temporary here"""
    colors = np.array(list(islice(cycle(['#377eb8', '#ff7f00', '#4daf4a',
                                         '#f781bf', '#a65628', '#984ea3',
                                         '#999999', '#e41a1c', '#dede00']),
                                  int(len(unique_labels) + 1))))
    # add black color for outliers (if any)
    colors = np.append(colors, ["#000000"])

    plt.scatter(X[:, 0], X[:, 1], color=colors[labels])
    plt.show()
    #plt.savefig("ED/clustering.png")
    plt.close()    
