import operator

import numpy as np
import matplotlib.pyplot as plt

from km3pipe import Module

class ZTPlotter(Module):
    def process(self, blob):
        self.scatter(blob['EvtRawHits'],
                     size=2, alpha=1.0, marker='.',
                     label="RawHits")
        self.scatter(blob['MergedEvtRawHits'],
                     size=80, alpha=0.2, color='red', marker='o',
                     label='MergedRawHits')
        self.scatter(blob['LongToTHits'],
                     size=60, alpha=1.0, color='green', marker='x',
                     label='Long ToT Hits')
        plt.xlabel("t [ns]", fontsize=12)
        plt.ylabel("z [m]", fontsize=12)
        plt.legend()
        plt.show()
        return blob

    def scatter(self, hits, color='blue', size=10, alpha=1.0, marker='o', label=None):
        times, zs = self.get_zt_points(hits)
        plt.scatter(times, zs, s=size, c=color, alpha=alpha, marker=marker, label=label)

    def get_zt_points(self, hits):
        times = [hit.time for hit in hits]
        zs = []
        for hit in hits:
            pmt = self.detector.pmt_with_id(hit.pmt_id)
            zs.append(pmt.pos.z)
        return times, zs


class T3HitSelector(Module):
    """Creates a list of hits which contribute to a T3 definition.

    T3 is defined by the coincidence of two hits on adjacent or
    next-to-adjacent OMs within a given time window.
    
    """
    def __init__(self, **context):
        super(self.__class__, self).__init__(**context)
        self.adjacent_t = self.get('adjacent_t') or 80
        self.next_to_adjacent_t = self.get('next_to_adjacent') or 160
        self.input_hits = self.get('input_hits') or 'LongToTHits'
        self.output_hits = self.get('output_hits') or 'T3Hits'

    def process(self, blob):
        om_hits = sort_hits_by_om(self.input_hits, self.detector)


class TOTFilter(Module):
    """Only keeps hits with at least min_tot"""
    def __init__(self, **context):
        super(self.__class__, self).__init__(**context)
        self.min_tot = self.get('min_tot') or 20 
        self.input_hits = self.get('input_hits') or 'EvtRawHits'
        self.output_hits = self.get('output_hits') or 'LongToTHits'

    def process(self, blob):
        hits = blob[self.input_hits]
        filtered_hits = [hit for hit in hits if hit.tot >= self.min_tot]
        blob[self.output_hits] = filtered_hits
        print("Number of long tot hits: " + str(len(filtered_hits)))
        return blob


class OMRawHitMerger(Module):
    """Merges hits on the same OM within a given time window.

    Each merged hit gets the id, pmt_id and time of the earliest contributed
    hit. The time-over-threshold values are added.
    
    """
    def __init__(self, **context):
        super(self.__class__, self).__init__(**context)
        self.time_window = self.get('time_window') or 10
        self.input_hits = self.get('input_hits') or 'EvtRawHits'
        self.output_hits = self.get('output_hits') or 'MergedEvtRawHits'

    def process(self, blob):
        hits = blob[self.input_hits]
        hits.sort(key=lambda x: x.time)

        print("Number of raw hits: " + str(len(hits)))

        om_hits = sort_hits_by_om(hits, self.detector)
        merged_hits = self.merge_hits(om_hits)
        print("Number of merged hits: " + str(len(merged_hits)))

        blob[self.output_hits] = merged_hits
        return blob

    def merge_hits(self, om_hits):
        """Merge hits on each om and return a flat list"""
        merged_hits = []
        for om, hits in om_hits.iteritems():
            hits_to_merge = None
            for hit in hits:
                if not hits_to_merge:
                    hits_to_merge = [hit]
                    continue
                else:
                    if hit.time - hits_to_merge[-1].time <= self.time_window:
                        hits_to_merge.append(hit)
                    else:
                        if len(hits_to_merge) > 1:
                            merged_hit = reduce(operator.add, hits_to_merge)
                            merged_hits.append(merged_hit)
                        hits_to_merge = [hit]
        return merged_hits


class ROyFitter(Module):
    pass


def sort_hits_by_om(hits, detector):
    om_hits = {}
    for hit in hits:
        omkey = detector.pmtid2omkey(hit.pmt_id)
        om_hits.setdefault((omkey[0], omkey[1]), []).append(hit)
    return om_hits

