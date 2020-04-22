import config as cf
import data_containers as dc

n_Crates = 12
n_CardPerCrate = 10
n_ChPerCard = 64
n_ChPerCrate = n_CardPerCrate * n_ChPerCard 
n_ChPerConnector = 8
HalfCrate = int(n_Crates/2)
QuartCrate = int(HalfCrate/2)
HalfCard = int(n_CardPerCrate/2)
HalfChPerCrate = int(n_ChPerCrate/2)


def check():
    if(len(dc.map_ped) > 0):
        del dc.map_ped[:]

def ChannelMapper():    
    check()
    for idaq in range(cf.n_ChanTot):
        crp, view, vchan = DAQToCRP(idaq)
        dc.map_ped.append(dc.pdmap(crp,view,vchan))
        """
        ev = cf.pdmap()
        ev.view = view
        ev.crp = crp
        ev.vchan = vchan
        cf.map_ped.append(ev)
        """

def DAQToCRP(daqch):
    
    crate = int(daqch/n_ChPerCrate)
    card  = int((daqch - crate * n_ChPerCrate)/n_ChPerCard)
    chcard = int(daqch - crate * n_ChPerCrate - card * n_ChPerCard)
    conn   = int(chcard/n_ChPerConnector)
    view = 0 if crate < HalfCrate else 1

    crp  = -1
    
    if(view==0):
        if(crate < QuartCrate):
            if(card < HalfCard):
                crp = 1
            else:
                crp = 2
        else:
            if(card < HalfCard):
                crp = 0
            else:
                crp = 3
    else:
        if((crate-HalfCrate) < QuartCrate):
            if(card < HalfCard):
                crp = 0
            else:
                crp = 1
        else:
            if(card < HalfCard):
                crp = 3
            else:
                crp = 2
            
    vchan = -1
    if(view==0):
        vchan = abs(chcard-conn*n_ChPerConnector - 7) + conn*n_ChPerConnector + (crate%QuartCrate)*HalfChPerCrate + (card if card < HalfCard else card-HalfCard)*n_ChPerCard
    else:
        vchan = abs(chcard-conn*n_ChPerConnector - 7) + conn*n_ChPerConnector + abs((crate%QuartCrate)-2)*HalfChPerCrate + abs((card if card < HalfCard else card-HalfCard)-4)*n_ChPerCard
    return int(crp), int(view), int(vchan)
