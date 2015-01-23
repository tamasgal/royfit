import operator

from km3pipe import Module

class OMRawHitMerger(Module):
    def process(self, blob):
        detector = blob['Detector']
        hits = blob['EvtRawHits']
        hits.sort(key=lambda x: x.time)

        print("Number of raw hits: " + str(len(hits)))

        om_hits = {}
        for hit in hits:
            omkey = detector.pmtid2omkey(hit.pmt_id)
            om_hits.setdefault((omkey[0], omkey[1]), []).append(hit)

        merged_om_hits = {}
        for om, hits in om_hits.iteritems():
            print "raw hits:", len(hits)
            hits_to_merge = None
            for hit in hits:
                print hit
                if not hits_to_merge:
                    hits_to_merge = [hit]
                    continue
                else:
                    if hit.time - hits_to_merge[-1].time <= 10:
                        hits_to_merge.append(hit)
                    else:
                        if len(hits_to_merge) > 1:
                            merged_hit = reduce(operator.add, hits_to_merge)
                            merged_om_hits.setdefault(om, []).append(merged_hit)
                        hits_to_merge = [hit]

            if merged_om_hits.has_key(om):
                print "merged hits:", len(merged_om_hits[om])
                for hit in merged_om_hits[om]:
                    print hit
        
        return blob
    

class ROyFitter(Module):
    pass

