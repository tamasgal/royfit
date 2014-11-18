# coding=utf-8
# Filename: misc_modules.py
"""
Some SeaTray modules to visualise data and other stuff.

"""
__author__ = "Tamas Gal"
__maintainer__ = "Tamas Gal"
__email__ = "tamas.gal@physik.uni-erlangen.de"
__all__ = ('ROyMonitor', 'ZTPlot', 'PrintFrameIndex', 'PrintLineHits',
           'Sleeper')

import datetime as dt
import socket
import json

from time import sleep

from icecube import icetray, dataclasses

import royfit

from i3kit.geometry import Detector
from i3kit.utilities import PulseAssistant

class ROyMonitor(icetray.I3Module):
    """Find and send any ROy parameters to ROyWeb via UDP."""

    def __init__(self, context): # pylint: disable=E1002
        super(self.__class__, self).__init__(context)
        self.AddParameter("ip",
                          "IP of the ROyWeb server",
                          "131.188.167.62")
        self.AddParameter("port",
                          "UDP Port of the ROyWeb server",
                          9999)
        self.AddOutBox("OutBox")

    def Configure(self): # pylint: disable=C0103,C0111
        self.ip = self.GetParameter("ip")
        self.port = self.GetParameter("port")
        self.event_number = 0

    def Physics(self, frame): # pylint: disable=C0103,C0111
        self.event_number += 1

        for frame_key in frame.keys():
            if frame_key.startswith("ROy"):
                self.send_parameter(frame, frame_key)
        self.PushFrame(frame)

    def send_parameter(self, frame, frame_key):
        parameter_name = frame_key[3:]
        print("Sending: {0}".format(frame_key))
        message = json.dumps(
            {'kind': 'parameter',
             'type': parameter_name,
             'description': '',
             'value': frame[frame_key].value,
             'time': '',
             'event_number': self.event_number})
        sock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
        sock.sendto(message, (self.ip, self.port))


class ConePlot(icetray.I3Module):
    """Plots the summed merged pulses PMT dirs. Just a test..."""
    def __init__(self, context): # pylint: disable=E1002
        super(self.__class__, self).__init__(context)
        self.AddOutBox("OutBox")
        self.AddParameter("pulse_map_name",
                          "The pulses map",
                          'MultiOMRecoPulseSeries')
        self.pulse_map_name = None
        self.detector = None

    def Configure(self): # pylint: disable=C0103,C0111
        self.pulse_map_name = self.GetParameter("pulse_map_name")

    def Geometry(self, frame): # pylint: disable=C0103,C0111
        geometry = frame['I3Geometry']
        self.detector = Detector(geometry)
        self.PushFrame(frame)

    def Physics(self, frame): # pylint: disable=C0103,C0111
        pulse_map = frame[self.pulse_map_name]
        from i3kit.utilities import PulseAssistant
        pulse_assistant = PulseAssistant(pulse_map)

        merged_pulses = pulse_assistant.get_merged_pulses(10)

        from mpl_toolkits.mplot3d import axes3d
        import matplotlib.pyplot as plt

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        summed_x = 0
        summed_y = 0
        summed_z = 0
        for omkey, pulse in merged_pulses:
            dir = self.detector.pmt_dir(omkey)
            x = dir.GetX()
            y = dir.GetY()
            z = dir.GetZ()
            summed_x += x
            summed_y += y
            summed_z += z
            ax.plot_wireframe([0, x], [0, y], [0, z], color='green')
        ax.plot_wireframe([0, summed_x], [0, summed_y], [0, summed_z], color='red')

        mctree = frame["AntMCTree"]
        muon = royfit.icy_muon_from_tree(mctree)
        muon_dir = muon.GetDir()
        x = muon_dir.GetX()
        y = muon_dir.GetY()
        z = muon_dir.GetZ()
        ax.plot_wireframe([0, x], [0, y], [0, z], color='blue')

        plt.show()

        self.PushFrame(frame)



class ZTPlot(icetray.I3Module):
    """Show a z-t-plot for a given line."""

    def __init__(self, context): # pylint: disable=E1002
        super(self.__class__, self).__init__(context)
        self.AddOutBox("OutBox")
        self.AddParameter("pulse_map_name",
                          "The pulses map",
                          'MultiOMRecoPulseSeries')
        self.pulse_map_name = None
        self.detector = None

    def Configure(self): # pylint: disable=C0103,C0111
        self.pulse_map_name = self.GetParameter("pulse_map_name")

    def Geometry(self, frame): # pylint: disable=C0103,C0111
        geometry = frame['I3Geometry']
        self.detector = Detector(geometry)
        self.PushFrame(frame)

    def scatter(self, pulse_info):
        ts = []
        zs = []
        for omkey, pulse in pulse_info:
            pos_z = self.detector.pmt_pos(omkey).Z
            time = pulse.GetTime()
            ts.append(time)
            zs.append(pos_z)
        return ts, zs

    def prepare_plot(self, title):
        import matplotlib.pyplot as plt
        fig = plt.figure(figsize=(6, 6))
        ax = fig.add_subplot(111)
        ax.set_title(title, fontsize=14)
        ax.set_xlabel("t [ns]", fontsize=12)
        ax.set_ylabel("z [m]", fontsize=12)
        return ax, plt

    def Physics(self, frame): # pylint: disable=C0103,C0111
        pulse_map = frame[self.pulse_map_name]
        from i3kit.utilities import PulseAssistant
        pulse_assistant = PulseAssistant(pulse_map)

        merged_pulses = pulse_assistant.get_merged_pulses(100)
        for line in self.detector.lines:
            pulse_info = pulse_assistant.get_pulses_for_line(line)




            if len(pulse_info) > 2:
                merged_line_pulses = pulse_assistant.get_merged_pulses_for_line(line)

                used_doms = []
                selected_pulses = []
                for omkey, pulse in merged_line_pulses:
                    if omkey[1] not in used_doms:
                        used_doms.append(omkey[1])
                        selected_pulses.append((omkey, pulse))



                if len(merged_line_pulses) > 8:
                    title = "z-t-Plot for line {0}".format(line)

                    ts, zs = self.scatter(pulse_info)
                    ts_merged, zs_merged = self.scatter(merged_line_pulses)
                    ts_selected, zs_selected = self.scatter(selected_pulses)

                    ax, plt = self.prepare_plot(title)
                    ax.scatter(ts_selected, zs_selected, s=150, marker="x",
                               alpha=0.4, c='red', label='selected pulses')
                    ax.scatter(ts_merged, zs_merged, s=150,
                               alpha=0.2, c='green', label='merged pulses')
                    ax.scatter(ts, zs, s=10, c='blue', label='PMT hits')

                    plt.legend()
                    plt.show()

        self.PushFrame(frame)






class PrintFrameIndex(icetray.I3Module):
    """Print frame index and time needed to process an event."""

    def __init__(self, context): # pylint: disable=E1002
        super(self.__class__, self).__init__(context)
        self.AddOutBox("OutBox")
        self.timestamp = dt.datetime.now()
        self.frame_index = 1

    def Configure(self): # pylint: disable=C0103,C0111
        pass

    def Process(self): # pylint: disable=C0103,C0111
        new_timestamp = dt.datetime.now()
        delta_t = (new_timestamp - self.timestamp).microseconds / 1e3
        print ("-- Frame:%7u ------ %9.3f ms ----------" %
               (self.frame_index, delta_t))
        self.timestamp = new_timestamp

        frame = self.PopFrame()

        tag = frame.Stop.id
        if tag is "I":
            self.PushFrame(frame)
        elif tag is 'D':
            self.DetectorStatus(frame)
        elif tag is 'G':
            self.Geometry(frame)
        elif tag is 'C':
            self.Calibration(frame)
        elif tag in 'PQ':
            self.PushFrame(frame)
        else:
            self.PushFrame(frame)

        self.frame_index += 1


class PrintLineHits(icetray.I3Module):
    """Print PMT hits for each line of the detector."""

    def __init__(self, context): # pylint: disable=E1002
        super(self.__class__, self).__init__(context)
        self.AddOutBox("OutBox")
        self.AddParameter("pulse_map_name",
                          "The pulses map",
                          'MultiOMRecoPulseSeries_L1')
        self.detector = None
        self.pulse_map_name = None

    def Geometry(self, frame): # pylint: disable=C0103,C0111
        geometry = frame['I3Geometry']
        self.detector = Detector(geometry)
        print self.detector
        self.PushFrame(frame)

    def Configure(self): # pylint: disable=C0103,C0111
        self.pulse_map_name = self.GetParameter("pulse_map_name")

    def Physics(self, frame): # pylint: disable=C0103,C0111
        pulse_map = frame[self.pulse_map_name]

        om_pulses = {}
        for omkey, pulses in pulse_map:
            line = omkey[0]
            om_ix = omkey[1]
            print pulses
            for pulse in pulses:
                charge = pulse.GetCharge()
                time = pulse.GetTime()
                om_pulses.setdefault((line, om_ix), []).append((time, charge))

        for line in self.detector.lines:
            output_line = "{0} ".format(line)
            output_line += str(self.detector.oms)
            print output_line

        selected_om_pulses = dict((k, v) for k, v in om_pulses.iteritems() \
                                  if len(v) > 1)

        print len(om_pulses)
        print len(selected_om_pulses)

        import matplotlib.pyplot as plt
        x_line = [i for i, _ in om_pulses.keys()]
        y_om = [i for _, i in om_pulses.keys()]
        # this one only uses the first time hit
        z_charge = [om_pulses[pos][0][0] for pos in zip(x_line, y_om)]
        for line in x_line:
            plt.axvline(x=x_line, ymin=0., linewidth=1, color='b')
        plt.scatter(x_line, y_om, c=z_charge, s=100,
                    cmap=plt.cm.cool, edgecolors='None', alpha=0.75)
        plt.colorbar()
        plt.show()

class Sleeper(icetray.I3Module):
    """Sleep for a given time in the Physics frame."""

    def __init__(self, context): # pylint: disable=E1002
        super(self.__class__, self).__init__(context)
        self.AddOutBox("OutBox")
        self.AddParameter("milliseconds",
                          "Milliseconds to sleep. Yawn...",
                          1000)

    def Configure(self): # pylint: disable=C0103,C0111
        self.milliseconds = self.GetParameter("milliseconds")

    def Physics(self, frame): # pylint: disable=C0103,C0111
        sleep(self.milliseconds / 1000.)
        self.PushFrame(frame)

if __name__ == '__main__':
    pass
