# coding=utf-8
# Filename: royfit.py
"""Main data structures and classes for the ROyFit module."""

__author__ = "Tamas Gal"
__maintainer__ = "Tamas Gal"
__email__ = "tamas.gal@physik.uni-erlangen.de"
__all__ = ('ROyFitter', 'PulseAssistant', 'Detector')

import numpy as np

from icecube import icetray, dataclasses
from i3kit.utilities import PulseAssistant
from i3kit.geometry import Detector
from i3kit.mctools import icy_muon_from_tree

class ROyFitter(icetray.I3Module):

    """The main fitter."""

    def __init__(self, context): # pylint: disable=E1002
        super(self.__class__, self).__init__(context)
        self.AddOutBox("OutBox")
        self.AddParameter("pulse_map_name",
                          "The pulses map",
                          'MultiOMRecoPulseSeries')
        self.detector = None
        self.pulse_map_name = None

    def Geometry(self, frame): # pylint: disable=C0103,C0111
        geometry = frame['I3Geometry']
        self.detector = Detector(geometry)
        self.PushFrame(frame)

    def Configure(self): # pylint: disable=C0103,C0111
        self.pulse_map_name = self.GetParameter("pulse_map_name")
        self.zeniths = []

    def Physics(self, frame): # pylint: disable=C0103,C0111
        pulse_map = frame[self.pulse_map_name]
        pulse_assistant = PulseAssistant(pulse_map)

        merged_pulses = pulse_assistant.get_merged_pulses(time_window=60)
        print "Merged pulses:", len(merged_pulses)
        frame.Put("ROyMergedPulses", dataclasses.I3Double(len(merged_pulses)))

        mctree = frame["I3MCTree"]
        muon = icy_muon_from_tree(mctree)
        zenith = muon.GetZenith() / icetray.I3Units.deg
        print "Muon zenith [deg]: ", zenith
        frame.Put("ROyMuonZenith_in_deg", dataclasses.I3Double(zenith))
        print "Muon energy: ", muon.GetEnergy()
        frame.Put("ROyMuonEnergy_in_GeV", dataclasses.I3Double(muon.GetEnergy()))

        print(69*"#")
        print("Let's try to fit!")

        from minimiser import QualityFunction
        import minuit

        for line in self.detector.lines:
            pulses = pulse_assistant.get_merged_pulses_for_line(line)

            used_doms = []
            selected_pulses = []
            for omkey, pulse in pulses:
                if omkey[1] not in used_doms:
                    used_doms.append(omkey[1])
                    selected_pulses.append((omkey, pulse))

            if len(selected_pulses) > 4:
                print "Trying to perform single track fit on line:", line

                pulse_times = [p.GetTime() for _, p in selected_pulses]
                #print "Pulse times:", pulse_times
                z_coordinates = [self.detector.pmt_pos(omkey).Z for omkey, _ in selected_pulses]
                #print "z-coordinates:", z_coordinates


                zc_ini = z_coordinates[0]
                tc_ini = pulse_times[0]
                dc_ini = 20.

                dc_lowlimit = 2.
                dc_highlimit = 80.

                quality_function = QualityFunction(pulse_times, z_coordinates)
                fitter = minuit.Minuit(quality_function, zc=zc_ini, tc=tc_ini, dc=dc_ini)

                fitter.limits["uz"] = (-0.999,0.999)
                fitter.limits["dc"] = (dc_lowlimit,dc_highlimit)

                fitter.errors["uz"] = 0.01#0.02
                fitter.errors["zc"] = 1#0.5
                fitter.errors["dc"] = 1#0.5
                fitter.errors["tc"] = 3.#5.

                fitter.tol = 1
                fitter.up = 1
                fitter.maxcalls = 500

                #fitter.printMode = 1
                try:
                    fitter.migrad()
                except minuit.MinuitError:
                    print("Fitting error!!!")
                else:
                    print fitter.values
                    print fitter.errors
                    reco_zenith = 180 - (np.arccos(fitter.values["uz"]) / icetray.I3Units.deg)
                    print "Reconstructed zenith:", reco_zenith
                    self.zeniths.append((zenith, reco_zenith))

        self.PushFrame(frame)

    def Finish(self):

        import numpy as np
        import matplotlib as mpl
        import matplotlib.pyplot as plt

        x = [item[0] for item in self.zeniths]
        y = [item[1] for item in self.zeniths]
        H, xedges, yedges = np.histogram2d(y, x, bins=90)
        fig = plt.figure(figsize=(7, 7))
        ax = fig.add_subplot(111)
        im = plt.imshow(H, interpolation=None, origin='low',
                        extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]])
        plt.xlabel("MC zenith [deg]")
        plt.ylabel("reco zenith [deg]")
        plt.show()



if __name__ == '__main__':
    pass
