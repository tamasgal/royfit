# coding=utf-8
# Filename: minimiser.py
"""
Description.

"""
from __future__ import division

__author__ = 'Tamas Gal'
__email__ = 'tamas.gal@physik.uni-erlangen.de'

import numpy as np
from icecube import icetray, dataclasses

n = dataclasses.I3Constants.n_water_antares
c = dataclasses.I3Constants.c

class SingleStringParameters(object):
    def __init__(self, t, z):
        self.t = np.array(t) # measured hit times
        self.z = np.array(z) # z-component of hit PMT

    def D_gamma(self, uz, zc, dc):
        return (n/np.sqrt(n**2 - 1))*np.sqrt(dc**2 + ((self.z-zc)**2)*(1-uz**2))

    def T_gamma(self, uz, zc, dc, tc):
        return tc+((self.z-zc)*uz+(n**2-1)*self.D_gamma(uz,zc,dc)/n)/c

    def Cos_theta(self, uz, zc, dc):
        return (1-uz**2)*(self.z-zc)/self.D_gamma(uz, zc, dc) + uz/n


class QualityFunction(SingleStringParameters):
    """Creates the quality function for the minimiser.
    """
    def __call__(self, uz, zc, dc, tc):
        return sum((self.T_gamma(uz, zc, dc, tc) - self.t)**2)
