import numpy as np



def field_moy(x, y):
    """ simple average field parametrization """
    """ to be improved """
    
    d = np.sqrt(pow(x,2)+pow(y,2))

    param = [155.416, -0.020118, 0.001116]
    E = 0.
    for i in range(len(param)):
        E += param[i] * pow(d,i)
    return E*1e-3
