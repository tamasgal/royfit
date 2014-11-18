# coding=utf-8
# Filename: test_royfit.py
# pylint: disable=C0103,C0111,R0201,R0903
"""
Unittests for the royfit module.

"""

__author__ = "Tamas Gal"
__maintainer__ = "Tamas Gal"
__email__ = "tamas.gal@physik.uni-erlangen.de"


class TestSeaTrayFramework(object):
    def test_create_particle(self):
        from icecube import icetray, dataclasses
        particle = dataclasses.I3Particle()
        particle.SetSpeed(1)
        assert particle.GetSpeed() == 1





if __name__ == '__main__':
    pass
