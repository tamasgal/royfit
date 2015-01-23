#!/usr/bin/env python

import os

try:
    import km3pipe
except ImportError:
    print("km3pipe could not be imported. Falling back to dev version.")
    import sys
    sys.path.insert(1, '/Users/tamasgal/Dev/km3pipe')
    import km3pipe

from km3pipe import Geometry
from km3pipe.pumps import EvtPump
from km3modules import StatusBar
from royfit import OMRawHitMerger

DATA_PATH='/Users/tamasgal/Data/KM3NeT'
EVT_FILE='mu_single_line/km3net_jul13_90m_muatm50T99_single_line.km3.jte.evt'
GEO_FILE='Detector/km3net_single_line.detx'

def main():
    pipe = km3pipe.Pipeline()
    pipe.attach(StatusBar)
    pipe.attach(EvtPump, filename=os.path.join(DATA_PATH, EVT_FILE))
    pipe.attach(Geometry, filename=os.path.join(DATA_PATH, GEO_FILE))
    pipe.attach(OMRawHitMerger, time_window=10)
    pipe.drain(4)


if __name__ == '__main__':
    main()
