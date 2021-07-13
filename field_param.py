import numpy as np



def field_moy(x, y):
    """ using squircles -- in agreement with comsol within 5% """
    a = 0.026
    b = 155.236
    E = pow(pow(a*x, 10) + pow(a*y, 10), 0.25) + b
    return E*1e-3
               

    """ simple average field parametrization """
    '''
    d = np.sqrt(pow(x,2)+pow(y,2))
    
    param = [155.416, -0.020118, 0.001116]
    E = 0.
    for i in range(len(param)):
        E += param[i] * pow(d,i)
    return E*1e-3
    '''
