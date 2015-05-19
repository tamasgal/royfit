# coding=utf-8
# Filename: __init__.py
"""
Restless Oyster Fit: Fast muon reconstruction for KM3NeT.

"""
from __future__ import division, absolute_import, print_function
from royfit.__version__ import version, version_info

__author__ = "Tamas Gal"
__copyright__ = ("Copyright 2015, Tamas Gal and the KM3NeT collaboration "
                 "(http://km3net.org)")
__credits__ = []
__license__ = "MIT"
__version__ = version
__maintainer__ = "Tamas Gal"
__email__ = "tgal@km3net.de"
__status__ = "Development"


from royfit.core import (ZTPlotter, T3HitSelector, TOTFilter, OMRawHitMerger,
                         FirstOMHitFilter, ROyFitter)

from royfit.minimiser import QualityFunction
