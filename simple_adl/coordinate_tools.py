#!/usr/bin/env python
"""
Generic python script.
New file for replacing ugali dependencies
"""
__author__ = "William Cerny"

import numpy as np

def distanceToDistanceModulus(distance):
    """
    Return distance modulus for a given distance (kpc).
    """
    return 5. * (np.log10(np.array(distance) * 1.e3) - 1.)

dist2mod = distanceToDistanceModulus

def distanceModulusToDistance(distance_modulus):
    """
    Return distance (kpc) for a given distance modulus.
    """
    return 10**((0.2 * np.array(distance_modulus)) - 2.)

mod2dist = distanceModulusToDistance

def angsep(lon1,lat1,lon2,lat2):
    """
    Angular separation (deg) between two sky coordinates.
    Borrowed from astropy (www.astropy.org)
    Notes
    -----
    The angular separation is calculated using the Vincenty formula [1],
    which is slighly more complex and computationally expensive than
    some alternatives, but is stable at at all distances, including the
    poles and antipodes.
    [1] http://en.wikipedia.org/wiki/Great-circle_distance
    """
    lon1,lat1 = np.radians([lon1,lat1])
    lon2,lat2 = np.radians([lon2,lat2])
    
    sdlon = np.sin(lon2 - lon1)
    cdlon = np.cos(lon2 - lon1)
    slat1 = np.sin(lat1)
    slat2 = np.sin(lat2)
    clat1 = np.cos(lat1)
    clat2 = np.cos(lat2)

    num1 = clat2 * sdlon
    num2 = clat1 * slat2 - slat1 * clat2 * cdlon
    denominator = slat1 * slat2 + clat1 * clat2 * cdlon

    return np.degrees(np.arctan2(np.hypot(num1,num2), denominator))
