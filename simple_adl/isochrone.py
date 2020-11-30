#!/usr/bin/env python
"""
Generic python script.

This module is an attempt to port over a small amount of the ugali isochrone functionality
to a smaller piece of code with fewer dependencies.
See https://github.com/DarkEnergySurvey/ugali/blob/master/ugali/isochrone/model.py
and https://github.com/DarkEnergySurvey/ugali/blob/master/ugali/isochrone/parsec.py
"""
__author__ = "Sidney Mau"

import numpy as np
import scipy.interpolate
from collections import OrderedDict as odict
import os

import simple_adl

#------------------------------------------------------------------------------

class Isochrone():
    def __init__(self, age=12.0, metallicity=0.00010, distance_modulus=16.0, survey='des', band_1='g', band_2='r', filename=os.path.join(os.path.dirname(simple_adl.__file__),'isochrones/iso_a12.0_z0.00010.dat')):
        self.age = age
        self.metallicity = metallicity
        self.distance_modulus = distance_modulus
        self.survey = survey
        self.band_1 = band_1
        self.band_1_detection = True
        self.band_2 = band_2
        #self.imf_type = Chabrier2003
        self.hb_stage = 4
        self.hb_spread = 0.1
        self.filename = filename
        self._parse(self.filename)

    columns = dict(
        des = odict([
                (3, ('mass_init',float)),
                (4, ('mass_act',float)),
                (5, ('log_lum',float)),
                (10,('g',float)),
                (11,('r',float)),
                (12,('i',float)),
                (13,('z',float)),
                (14,('Y',float)),
                (16,('stage',int)),
                ]),
        sdss = odict([
                (3, ('mass_init',float)),
                (4, ('mass_act',float)),
                (5, ('log_lum',float)),
                (9, ('u',float)),
                (10,('g',float)),
                (11,('r',float)),
                (12,('i',float)),
                (13,('z',float)),
                (15,('stage',int)),
                ]),
        ps1 = odict([
                (3, ('mass_init',float)),
                (4, ('mass_act',float)),
                (5, ('log_lum',float)),
                (9, ('g',float)),
                (10,('r',float)),
                (11,('i',float)),
                (12,('z',float)),
                (13,('y',float)),
                (16,('stage',int)),
                ]),
        lsst = odict([
                (3, ('mass_init',float)),
                (4, ('mass_act',float)),
                (5, ('log_lum',float)),
                (9, ('u',float)),
                (10,('g',float)),
                (11,('r',float)),
                (12,('i',float)),
                (13,('z',float)),
                (14,('Y',float)),
                (16,('stage',float))
                ]),
        )


    def _parse(self,filename):
        """
        Adapted from https://github.com/DarkEnergySurvey/ugali/blob/master/ugali/isochrone/parsec.py
        """
        try:
            columns = self.columns[self.survey.lower()]
        except KeyError as e:
            logger.warning('Unrecognized survey: %s'%(survey))
            raise(e)

        # delimiter='\t' is used to be compatible with OldPadova...
        # ADW: This should be updated, but be careful of column numbering
        kwargs = dict(delimiter='\t',usecols=list(columns.keys()),
                      dtype=list(columns.values()))
        self.data = np.genfromtxt(filename,**kwargs)

        self.mass_init = self.data['mass_init']
        self.mass_act  = self.data['mass_act']
        self.luminosity = 10**self.data['log_lum']
        self.mag_1 = self.data[self.band_1]
        self.mag_2 = self.data[self.band_2]
        self.stage = self.data['stage']

        self.mass_init_upper_bound = np.max(self.mass_init)
        self.index = len(self.mass_init)

        self.mag = self.mag_1 if self.band_1_detection else self.mag_2
        self.color = self.mag_1 - self.mag_2

    def separation(self, mag_1, mag_2):
        """
        Adapted from https://github.com/DarkEnergySurvey/ugali/blob/master/ugali/isochrone/model.py
        Parameters:
        -----------
        mag_1 : The magnitude of the test points in the first band
        mag_2 : The magnitude of the test points in the second band
        Returns:
        --------
        sep : Minimum separation between test points and isochrone interpolation
        """

        iso_mag_1 = self.mag_1 + self.distance_modulus
        iso_mag_2 = self.mag_2 + self.distance_modulus

        def interp_iso(iso_mag_1,iso_mag_2,mag_1,mag_2):
            interp_1 = scipy.interpolate.interp1d(iso_mag_1,iso_mag_2,bounds_error=False)
            interp_2 = scipy.interpolate.interp1d(iso_mag_2,iso_mag_1,bounds_error=False)

            dy = interp_1(mag_1) - mag_2
            dx = interp_2(mag_2) - mag_1

            dmag_1 = np.fabs(dx*dy) / (dx**2 + dy**2) * dy
            dmag_2 = np.fabs(dx*dy) / (dx**2 + dy**2) * dx

            return dmag_1, dmag_2

        # Separate the various stellar evolution stages
        if np.issubdtype(self.stage.dtype,np.number):
            sel = (self.stage < self.hb_stage)
        else:
            sel = (self.stage != self.hb_stage)

        # First do the MS/RGB
        rgb_mag_1 = iso_mag_1[sel]
        rgb_mag_2 = iso_mag_2[sel]
        dmag_1,dmag_2 = interp_iso(rgb_mag_1,rgb_mag_2,mag_1,mag_2)

        # Then do the HB (if it exists)
        if not np.all(sel):
            hb_mag_1 = iso_mag_1[~sel]
            hb_mag_2 = iso_mag_2[~sel]

            hb_dmag_1,hb_dmag_2 = interp_iso(hb_mag_1,hb_mag_2,mag_1,mag_2)

            dmag_1 = np.nanmin([dmag_1,hb_dmag_1],axis=0)
            dmag_2 = np.nanmin([dmag_2,hb_dmag_2],axis=0)

        #return dmag_1,dmag_2
        return(np.sqrt(dmag_1**2 + dmag_2**2))
