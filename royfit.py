import operator

from km3pipe import Module


class OMRawHitMerger(Module):
    """Merges hits on the same OM within a given time window.

    Each merged hit gets the id, pmt_id and time of the earliest contributed
    hit. The time-over-threshold values are added.
    
    """
    def __init__(self, **context):
        super(self.__class__, self).__init__(**context)
        self.time_window = self.get('time_window') or 10
        self.input_hits_key = self.get('input_hits_key') or 'EvtRawHits'
        self.output_hits_key = self.get('output_hits_key') or 'MergedEvtRawHits'

    def process(self, blob):
        hits = blob[self.input_hits_key]
        hits.sort(key=lambda x: x.time)

        print("Number of raw hits: " + str(len(hits)))

        om_hits = self.sort_hits_by_om(hits)
        merged_hits = self.merge_hits(om_hits)
        print("Number of merged hits: " + str(len(merged_hits)))
        
        blob[self.output_hits_key] = merged_hits
        return blob

    def sort_hits_by_om(self, hits):
        om_hits = {}
        for hit in hits:
            omkey = self.detector.pmtid2omkey(hit.pmt_id)
            om_hits.setdefault((omkey[0], omkey[1]), []).append(hit)
        return om_hits

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

