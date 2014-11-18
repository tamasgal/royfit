#!/usr/bin/env python
# coding=utf-8
# Filename: ice_rofit.py
"""
Yet another BBFit implementation for KM3NeT data.
Based on an idea published by ANTARES http://arxiv.org/abs/1105.4116

"""
__author__ = "Tamas Gal"
__copyright__ = ("Copyright 2014, Tamas Gal and the KM3NeT collaboration "
                 "(http://km3net.org)")
__credits__ = []
__license__ = "GPL"
__version__ = "0.0.1"
__maintainer__ = "Tamas Gal"
__email__ = "tamas.gal@physik.uni-erlangen.de"
__status__ = "Development"


def main():
    """Run a basic reconstruction on a test file."""
    from I3Tray import I3Tray, load
    from icecube import dataclasses, dataio, phys_services

#    inputfile = "/Users/tamasgal/Data/KM3NeT/MC/modk40_v4_numuCC_23.i3"
    inputfile = "/Users/tamasgal/Data/KM3NeT/v4/modk40_v4_numuCC_69.i3"
#    inputfile = "/Users/tamasgal/Data/KM3NeT/2khits/hex115_3inch31pm90_1836_muatm10T13.km3_v5r1.modk40h_Nhit_2k.i3"
    outputfile = "out.i3"

    tray = I3Tray()

    tray.AddModule("I3Reader", "reader", Filename=inputfile)



    from tgseatray.event_selectors import NeutrinoEnergyFilter
    tray.AddModule(NeutrinoEnergyFilter, "energyfilter",
                   MinEnergy=3000, MaxEnergy=1000000)

    from royfit.misc_modules import PrintFrameIndex
    tray.AddModule(PrintFrameIndex, "PrintFrameIndex")



    #from royfit.misc_modules import Sleeper
    #tray.AddModule(Sleeper, "Sleeper")

    from royfit.misc_modules import ZTPlot
    tray.AddModule(ZTPlot, "ZTPlot")

    from royfit.royfit import ROyFitter
    tray.AddModule(ROyFitter, "ROyFitter")

    #from royfit.misc_modules import ROyMonitor
    #tray.AddModule(ROyMonitor, "ROyMonitor")

    #from royfit.misc_modules import ConePlot
    #tray.AddModule(ConePlot, "ConePlot")

    tray.AddModule("I3Writer", "writer")(("filename", outputfile))
    tray.AddModule("TrashCan", "the can")

    tray.Execute()
    tray.Finish()


if __name__ == '__main__':
    main()
