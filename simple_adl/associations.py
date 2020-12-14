#!/usr/bin/env python
"""
Generic python script.
"""
__author__ = "William Cerny and Sidney Mau"

import numpy as np
import astropy
from astropy.coordinates import SkyCoord
import astropy.units as u
#-------------------------------------------------------------------------------


def get_cat_dir():
    catdir = './external_catalogs'
    if not os.path.exists(catdir):
        msg = "Catalog directory not found:\n%s"%catdir

    return catdir




Class ExternalCatalog(object): 
  
