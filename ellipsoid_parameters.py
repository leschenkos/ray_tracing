import numpy as np
Pi=np.pi

def ellipsoid_parameters(p,q,alfa):
    """p: input arm, q: output arm, alfa: angle of incidence in degree"""
    a = (p+q)/2
    b = np.sqrt(p*q*(1+np.cos(2*alfa/180*Pi))/2)
    e=(1-b**2/a**2)**0.5
    Y00=p*q*np.sin(2*alfa/180*Pi)/2/a/e
    X00=a*(1-Y00**2/b**2)**0.5
    return a, b, X00

P=ellipsoid_parameters(1000,250,85)
print( 'ellipsoid: ' 'a: ',P[0],' b: ', P[1], ' off-axis position xm ', P[2])