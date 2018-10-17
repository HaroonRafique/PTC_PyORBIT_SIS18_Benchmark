import numpy as np
from orbit.utils.consts import speed_of_light

class DoubleRF(object):
   
    def __init__(self, R, eta, beta, E, phi, h, V):
        cos = np.cos
        sin = np.sin
        pi = np.pi
        sqrt = np.sqrt
        c = speed_of_light

        phi = np.atleast_1d(phi)
        h = np.atleast_1d(h)
        V = np.atleast_1d(V)
        self.beta = beta
        self.eta = eta
        self.T = lambda dp: 0.5*beta*c*eta*dp**2
        
        def v(V, h, phi):
            return lambda z: V/h*(cos(h*z/R+phi) - cos(phi))     
        self.V = lambda z: 1./(E*beta/c)/2/pi * reduce(lambda x, y: x+y, [v(*i)(z) for i in zip(V, h, phi)])
    
    def get_H(self,z,dp):
        return self.T(dp) + self.V(z)
    
    def get_dp(self,z,H):
    	c = speed_of_light
        return np.abs(np.sqrt((H-self.V(z))/(0.5*self.beta*c*self.eta)))