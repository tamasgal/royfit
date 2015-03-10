import operator
import pickle

import numpy as np
import matplotlib.pyplot as plt

from km3pipe import Module
from km3pipe.tools import unit_vector, angle_between

from minimiser import QualityFunction

import iminuit as minuit
#import minuit2 as minuit


class ZTPlotter(Module):
    """A z-t-plotter"""
    #TODO: rewrite me, I'm hard coded!
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
    """Keeps only the first hit on each OM"""
    def __init__(self, **context):
        super(self.__class__, self).__init__(**context)
        self.input_hits = self.get('input_hits') or 'LongToTHits'
        self.output_hits = self.get('output_hits') or 'FirstOMHits'
        self.zeniths = []
        self.all_zeniths = []
        self.quality_parameters = []
        self.stats = []
        self.processed_events = 0
        self.tried_events = 0

    def process(self, blob):
        self.processed_events += 1
        mc_track = blob['TrackIns'][0]
        zenith = mc_track.dir.zenith * 180.0 / np.pi

        self.all_zeniths.append(zenith)
        

        hits = blob[self.input_hits]
        if len(hits) < 5:
            return blob
        self.tried_events += 1

        hit_times = [hit.time for hit in hits]

        z_coordinates = []
        for hit in hits:
            pmt = self.detector.pmt_with_id(hit.pmt_id)
            z_coordinates.append(pmt.pos.z)


        zc_ini = (min(z_coordinates) + max(z_coordinates)) / 2#sum(z_coordinates)/len(z_coordinates)
        tc_ini = min(hit_times)#sum(hit_times)/len(hit_times)
        dc_ini = 20.
        uz_ini = -0.75 #np.pi - np.cos(zenith)

        dc_lowlimit = 2.
        dc_highlimit = 100.

        quality_function = QualityFunction(hit_times, z_coordinates, 8)
        fitter = minuit.Minuit(quality_function,
                               zc=zc_ini,
                               tc=tc_ini,
                               dc=dc_ini,
                               uz=uz_ini,
                               error_dc=1.0,
                               error_uz=0.01,
                               error_tc=1.0,
                               error_zc=1.0,
                               limit_zc=(min(z_coordinates),
                                         max(z_coordinates)),
                               limit_uz = (-1.0, 1.0),
                               limit_dc = (dc_lowlimit,dc_highlimit))

        fitter.tol = 1
        #fitter.up = 1
        #fitter.maxcalls = 500

        #fitter.printMode = 1
        try:
            fitter.migrad()
        except minuit.MinuitError:
            print("Fitting error!!!")
        else:
            quality_parameter = fitter.fval / 4
            print("Q/4: {0}".format(quality_parameter))
            print("Values:")
            print(fitter.values)
            print("Errors:")
            print(fitter.errors)
            print("MC zenith: {0}".format(zenith))
            reco_zenith = 180 - (np.arccos(fitter.values["uz"]) / (np.pi/180.0))
            #reco_zenith = fitter.values["uz"]
            print("Reconstructed zenith: {0}".format(reco_zenith))

            if fitter.get_fmin().is_valid:
                self.quality_parameters.append(quality_parameter)
                self.zeniths.append((zenith, reco_zenith))
                self.stats.append((zenith, reco_zenith, quality_parameter))

#        x, y = fitter.profile('zc', subtract_min=True)
#        plt.plot(x, y)
#        plt.show()

        return blob

    def finish(self):
        print("Processed {0} events".format(self.processed_events))
        print("Tried fit on {0} events".format(self.tried_events))

        mc_zenith = [zenith for zenith, reco_zenith in self.zeniths]
        reco_zenith = [reco_zenith for zenith, reco_zenith in self.zeniths]

        plt.hist2d(mc_zenith, reco_zenith, bins=30)
        plt.title("ROyFit on {0} MUPAGE events with {1} valid fits" \
                  .format(self.processed_events, len(self.zeniths)),
                  y=1.04)
        plt.xlabel('MC zenith [degree]')
        plt.ylabel('reco zenith [degree]')
        plt.colorbar()
        plt.show()

        plt.hist(self.quality_parameters, bins=40)
        plt.title("ROyFit on {0} MUPAGE events with {1} valid fits - Q/4" \
                  .format(self.processed_events, len(self.zeniths)),
                  y=1.04)
        plt.xlabel("Q/4")
        plt.ylabel("count")
        plt.show()


        plt.hist([zen_mc - zen_reco for zen_mc, zen_reco in self.zeniths])
        plt.show()

        with open('reco_stats.pickle', 'w') as file:
            pickle.dump(self.stats, file)

        with open('all_zeniths.pickle', 'w') as file:
            pickle.dump(self.all_zeniths, file)

        with open('reco_zeniths.pickle', 'w') as file:
            pickle.dump(self.zeniths, file)



def sort_hits_by_om(hits, detector):
    om_hits = {}
    for hit in hits:
        omkey = detector.pmtid2omkey(hit.pmt_id)
        om_hits.setdefault((omkey[0], omkey[1]), []).append(hit)
    return om_hits


