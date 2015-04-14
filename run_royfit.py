#!/usr/bin/env python

import os

try:
    import km3pipe
except ImportError:
    print("km3pipe could not be imported. Falling back to dev version.")
    import sys
    sys.path.insert(1, '/sps/km3net/users/jreubelt/dev/km3pipe')
    import km3pipe

from km3pipe import Geometry
from km3pipe.pumps import EvtPump
from km3modules import StatusBar
from royfit import (OMRawHitMerger, TOTFilter, T3HitSelector, ZTPlotter,
                    FirstOMHitFilter, ROyFitter)

DATA_PATH='/sps/km3net/users/jreubelt/dev/'
#EVT_FILE='mu_single_line/km3net_jul13_90m_muatm50T99_single_line.km3.jte.evt'
EVT_FILE='playground/single_line_prod_jte/old_jte/km3net_jul13_90m_muatm10T'
#EVT_FILE='playground/single_line_prod_jte/old_jte/km3net_jul13_90m_muatm10T97.km3.jte'
GEO_FILE='detector/km3net_single_line.detx'

def main():
    pipe = km3pipe.Pipeline()
        #    pipe.attach(EvtPump, filename=os.path.join(DATA_PATH, EVT_FILE))
    pipe.attach(EvtPump, basename=os.path.join(DATA_PATH, EVT_FILE), index_start=10, index_stop=100)
    pipe.attach(StatusBar)
    pipe.attach(Geometry, filename=os.path.join(DATA_PATH, GEO_FILE))
    pipe.attach(OMRawHitMerger, time_window=20)
    pipe.attach(TOTFilter,
                input_hits='MergedEvtRawHits',
                output_hits='LongToTHits',
                min_tot=20)
    pipe.attach(FirstOMHitFilter,
                input_hits='LongToTHits',
                output_hits='FirstOMHits')
    pipe.attach(T3HitSelector, adjacent_t=200, next_to_adjacent_t=400,
                input_hits='FirstOMHits',
                candidate_hits='EvtRawHits',
                output_hits='T3Hits')
    pipe.attach(FirstOMHitFilter,
                input_hits='T3Hits',
                output_hits='FirstT3Hits')
#    pipe.attach(AdditionalHitSelector,
#                input_hits='FirstT3Hits',
#                hit_pool='EvtRawHits',
#                output_hits='ROyFitHits')
    pipe.attach(ROyFitter,
                input_hits='FirstT3Hits',
                output_track='ROyMuonTrack')
                #  pipe.attach(ZTPlotter)
    pipe.drain(10000)


if __name__ == '__main__':
    main()
