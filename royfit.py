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
                     size=80, alpha=0.2, color='blue', marker='o',
                     label='MergedRawHits')
        self.scatter(blob['LongToTHits'],
                     size=60, alpha=1.0, color='green', marker='.',
                     label='Long ToT Hits')
        self.scatter(blob['FirstOMHits'],
                     size=60, alpha=1.0, color='red', marker='x',
                     label='First OM Hits')
        self.scatter(blob['T3Hits'],
                     size=60, alpha=1.0, color='yellow', marker='.',
                     label='T3 Hits')
        self.scatter(blob['FirstT3Hits'],
                     size=160, alpha=1.0, color='magenta', marker='+',
                     label='First T3 Hits')
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
        self.adjacent_t = self.get('adjacent_t') or 200
        self.next_to_adjacent_t = self.get('next_to_adjacent') or 400
        self.input_hits = self.get('input_hits') or 'FirstOMHits'
        self.candidate_hits = self.get('candidate_hits') or 'EvtRawHits'
        self.output_hits = self.get('output_hits') or 'T3Hits'

    def process(self, blob):
        om_hits = sort_hits_by_om(blob[self.input_hits], self.detector)
        candidate_hits = blob[self.candidate_hits]
        om_candidate_hits = sort_hits_by_om(candidate_hits, self.detector)
        selected_hits = []
        skip_pmtids = []
        for om, hits in om_hits.iteritems():
            hit = hits[0]

            hits_above = om_candidate_hits.get((om[0], om[1]+1))
            hits_below = om_candidate_hits.get((om[0], om[1]-1))
            adjacent_hits = []
            if hits_above: adjacent_hits += hits_above
            if hits_below: adjacent_hits += hits_below

            hits_next_to_above = om_candidate_hits.get((om[0], om[1]+2))
            hits_next_to_below = om_candidate_hits.get((om[0], om[1]-2))
            next_to_adjacent_hits = []
            if hits_next_to_above: next_to_adjacent_hits += hits_next_to_above
            if hits_next_to_below: next_to_adjacent_hits += hits_next_to_below

            for adjacent_hit in adjacent_hits:
                if 0 <= adjacent_hit.time - hit.time <= self.adjacent_t:
                    if not adjacent_hit.pmt_id in skip_pmtids:
                        selected_hits.append(adjacent_hit)
                        skip_pmtids.append(adjacent_hit.pmt_id)
                    if not hit.pmt_id in skip_pmtids:
                        selected_hits.append(hit)
                        skip_pmtids.append(hit.pmt_id)

            for adjacent_hit in next_to_adjacent_hits:
                if 0 <= adjacent_hit.time - hit.time <= self.next_to_adjacent_t:
                    if not adjacent_hit.pmt_id in skip_pmtids:
                        selected_hits.append(adjacent_hit)
                        skip_pmtids.append(adjacent_hit.pmt_id)
                    if not hit.pmt_id in skip_pmtids:
                        selected_hits.append(hit)
                        skip_pmtids.append(hit.pmt_id)


        blob[self.output_hits] = selected_hits
        return blob



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


class FirstOMHitFilter(Module):
    """Keeps only the first hit on each OM"""
    def __init__(self, **context):
        super(self.__class__, self).__init__(**context)
        self.input_hits = self.get('input_hits') or 'LongToTHits'
        self.output_hits = self.get('output_hits') or 'FirstOMHits'

    def process(self, blob):
        hits = sorted(blob[self.input_hits], key=lambda x: x.time)
        first_hits = []
        skip_oms = []
        for hit in hits:
            om = self.detector.pmtid2omkey(hit.pmt_id)[:2]
            if om in skip_oms:
                continue
            first_hits.append(hit)
            skip_oms.append(om)
        blob[self.output_hits] = first_hits
        return blob


class ROyFitter(Module):
    pass


def sort_hits_by_om(hits, detector):
    om_hits = {}
    for hit in hits:
        omkey = detector.pmtid2omkey(hit.pmt_id)
        om_hits.setdefault((omkey[0], omkey[1]), []).append(hit)
    return om_hits

