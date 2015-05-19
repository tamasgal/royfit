# coding=utf-8
# Filename: minimiser.py
"""
Description.

"""
from __future__ import division, absolute_import, print_function

__author__ = 'Tamas Gal'
__email__ = 'tamas.gal@physik.uni-erlangen.de'

import numpy as np
from km3pipe import constants


n = 1.3797
c = constants.c / 1e9

class SingleStringParameters(object):
    def __init__(self, t, z, c, sigma_t=10, d0=50, d1=5):
        self.t = np.array(t) # measured hit times
        self.z = np.array(z) # z-component of hit PMT
        self.c = np.array(c) # charges or pmt_hit_counts (within a given deltat)
        self.c_mean = sum(self.c) / len(self.c)
        self.sigma_t = sigma_t # time error in ns
        self.d0 = d0 # distance to photon which induces 1pe
        self.d1 = d1 # minimum distance

    def D_gamma(self, uz, zc, dc):
        """Travel path"""
        return (n/np.sqrt(n**2 - 1))*np.sqrt(dc**2 + ((self.z-zc)**2)*(1-uz**2))

    def T_gamma(self, uz, zc, dc, tc):
        """Arrival time of a Cherenkov photon"""
        return tc+((self.z-zc)*uz+(n**2-1)*self.D_gamma(uz,zc,dc)/n)/c

    def Cos_theta(self, uz, zc, dc):
        """Inclination with respect to (0, 0, 1)"""
        return (1-uz**2)*(self.z-zc)/self.D_gamma(uz, zc, dc) + uz/n


class QualityFunction(SingleStringParameters):
    """Creates the quality function for the minimiser.
    """
    def __call__(self, uz, zc, dc, tc):
        # no weighting
        #return sum((self.T_gamma(uz, zc, dc, tc) - self.t)**2 / self.sigma_t**2 + self.c / self.c_mean)

        # first weighting
        #return sum((self.T_gamma(uz, zc, dc, tc) - self.t)**2 / self.sigma_t**2 + (self.c * np.sqrt(self.d1**2 + self.D_gamma(uz, zc, dc)**2))/(self.c_mean * self.d0))

        aip = 2.*self.c/(self.Cos_theta(uz, zc, dc) + 1.)
        avg_aip = sum(aip)/len(aip)
        D  = np.sqrt(self.d1**2 + self.D_gamma(uz, zc, dc)**2)
        return sum((self.T_gamma( uz, zc, dc, tc) - self.t)**2 / self.sigma_t**2 + aip*D/(avg_aip*self.d0))
