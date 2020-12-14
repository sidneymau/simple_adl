#!/usr/bin/env python
"""
Generic python script.
"""
__author__ = "William Cerny and Sidney Mau (based loosely on ugali's associate.py)"

import numpy as np
import astropy
from astropy.coordinates import SkyCoord
import astropy.units as u
import projector
from projector import gal2cel, cel2gal
import os,sys
from os.path import join,abspath,split
import inspect
from collections import OrderedDict as odict
import pandas as pd 
import fitsio

#-------------------------------------------------------------------------------


def get_cat_dir():
    catdir = './external_catalogs'
    if not os.path.exists(catdir):
        msg = "Catalog directory not found:\n%s"%catdir
    return catdir





class ExternalCatalog(object):
    DATADIR=get_cat_dir()
 
    def __init__(self, filename=None):
        columns = [('ObjectName',object),
                   ('RA',float),
                   ('DEC',float),
                   ('glon',float),
                   ('glat',float)]
        self.data = np.recarray(0,dtype=columns)
        self._load(filename)
        if np.isnan([self.data['glon'],self.data['glat']]).any():
            raise ValueError("Incompatible values")
 
    def __getitem__(self, key):
        """ 
        Support indexing, slicing and direct access.
        """
        try:
            return self.data[key]
        except ValueError as message:
            if key in self.data['name']:
                return self.data[self.data['name'] == key]
            else:
                raise ValueError(message)
 
    def __add__(self, other):
        ret = SourceCatalog()
        ret.data = np.concatenate([self.data,other.data])
        return ret
        
    def __len__(self):
        """ Return the length of the collection.
        """
        return len(self.data)
 
    def _load(self,filename):
        pass
 
    def match(self,lon,lat,coord='gal',tol=0.1,nnearest=1):
        if coord.lower() == 'cel':
            glon, glat = cel2gal(lon,lat)
        else:
            glon,glat = lon, lat
        return projector.match(glon,glat,self['glon'],self['glat'],tol,nnearest)
   
class McConnachie12(SourceCatalog):
    """
    Catalog of nearby dwarf spheroidal galaxies. 
    http://arxiv.org/abs/1204.1562
    https://www.astrosci.ca/users/alan/Nearby_Dwarfs_Database_files/NearbyGalaxies.dat
    """
 
    def _load(self,filename):
        if filename is None: 
            filename = os.path.join(self.DATADIR,"J_AJ_144_4/NearbyGalaxies2012.dat")
        self.filename = filename
 
        raw = np.genfromtxt(filename,delimiter=[19,3,3,5,3,3,3],usecols=range(7),dtype=['|S19']+6*[float],skip_header=36)
 
        self.data.resize(len(raw))
        self.data['name'] = np.char.strip(raw['f0'])
 
        ra = raw[['f1','f2','f3']].view(float).reshape(len(raw),-1)
        dec = raw[['f4','f5','f6']].view(float).reshape(len(raw),-1)
        self.data['ra'] = ugali.utils.projector.hms2dec(ra)
        self.data['dec'] = ugali.utils.projector.dms2dec(dec)
        
        glon,glat = cel2gal(self.data['ra'],self.data['dec'])
        self.data['glon'],self.data['glat'] = glon,glat

class McConnachie15(SourceCatalog):
    """
    Catalog of nearby dwarf spheroidal galaxies. Updated September 2015.
    http://arxiv.org/abs/1204.1562
    http://www.astro.uvic.ca/~alan/Nearby_Dwarf_Database_files/NearbyGalaxies.dat
    """
 
    def _load(self,filename):
        if filename is None: 
            filename = os.path.join(self.DATADIR,"J_AJ_144_4/NearbyGalaxies.dat")
        self.filename = filename
 
        raw = np.genfromtxt(filename,delimiter=[19,3,3,5,3,3,3],usecols=list(range(7)),dtype=['|S19']+6*[float],skip_header=36)

        self.data.resize(len(raw))
        self.data['name'] = np.char.lstrip(np.char.strip(raw['f0']),'*')

        ra = raw[['f1','f2','f3']].view(float).reshape(len(raw),-1)
        dec = raw[['f4','f5','f6']].view(float).reshape(len(raw),-1)
        self.data['ra'] = ugali.utils.projector.hms2dec(ra)
        self.data['dec'] = ugali.utils.projector.dms2dec(dec)
        
        glon,glat = cel2gal(self.data['ra'],self.data['dec'])
        self.data['glon'],self.data['glat'] = glon,glat




















