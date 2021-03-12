import config as cf


n_Crates = 12
n_CardPerCrate = 10
n_ChPerCard = 64
n_ChPerCrate = n_CardPerCrate * n_ChPerCard 
n_ChPerConnector = 8
HalfCrate = int(n_Crates/2)
QuartCrate = int(HalfCrate/2)
HalfCard = int(n_CardPerCrate/2)
HalfChPerCrate = int(n_ChPerCrate/2)


def ChannelMapper():    
    import data_containers as dc

    if(len(dc.map_ped) > 0):
        del dc.map_ped[:]

    for idaq in range(cf.n_ChanTot):
        crp, view, vchan = DAQToCRP(idaq)
        dc.map_ped.append(dc.pdmap(crp,view,vchan))

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



def CRPToDAQ(crp, view, vchannel):
    """ please don't ask how it works """

    #get the crate number
    if(view == 0):
        crate = int(vchannel/HalfChPerCrate)
        if(crp==0 or crp==3):
            crate += 3
    else:
        crate = int(abs(vchannel-959)/HalfChPerCrate)
        if(crp==0 or crp==1):
            crate += 6
        else : 
            crate += 9

    #get the card number
    if(view==0):
        card = int((vchannel - crate%QuartCrate*HalfChPerCrate)/n_ChPerCard)
        if(crp == 3 or crp == 2):
            card += 5
    else:
        card = int(abs(vchannel-abs(crate%QuartCrate-2)*HalfChPerCrate-319)/n_ChPerCard)
        if(crp == 1 or crp == 2):
            card += 5


    #get the connector
    if(view==0):
        conn = int((vchannel - crate%QuartCrate*HalfChPerCrate-(card if card < 5 else card-5)*n_ChPerCard)/n_ChPerConnector)
    else:
        conn = int(abs(vchannel-(abs(crate%QuartCrate-2)*HalfChPerCrate)-abs(((card if card < 5 else card-5)-4)*n_ChPerCard))/n_ChPerConnector)




    #get the channel card
    if(view==0):
        chcard = abs(int(vchannel - crate%QuartCrate*HalfChPerCrate-(card if card < 5 else card-5)*n_ChPerCard-conn*n_ChPerConnector)-7)+conn*n_ChPerConnector
    else:
        chcard = abs(int(vchannel - (abs(crate%QuartCrate-2)*HalfChPerCrate)-abs((card if card < 5 else card-5)-4)*n_ChPerCard-conn*n_ChPerConnector)-7)+conn*n_ChPerConnector


    #finally the daq channel
    if(view==0):        
        daqch = crate*n_ChPerCrate + card*n_ChPerCard + chcard
    else:
        daqch = n_ChPerCrate*HalfCrate + (crate-HalfCrate)*n_ChPerCrate + card*n_ChPerCard + chcard
    return daqch
