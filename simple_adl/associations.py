#!/usr/bin/env python
"""
Generic python script.
"""
__author__ = "William Cerny and Sidney Mau (based very loosely on ugali's associate.py)"

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
    DATADIR= get_cat_dir()
    
    
    def __init__(self, filename=None):
        columns = [('ObjectName',object),
                   ('RA',float),
                   ('DEC',float),
                   ('Classification',object)]
        
        self.data = np.recarray(0,dtype=columns)
        self._load(filename)
        if np.isnan([self.data['glon'],self.data['glat']]).any():
            raise ValueError("Incompatible values")
 
#     def __getitem__(self, key):
#         """ 
#         Support indexing, slicing and direct access.
#         """
#         try:
#             return self.data[key]
#         except ValueError as message:
#             if key in self.data['name']:
#                 return self.data[self.data['name'] == key]
#             else:
#                 raise ValueError(message)
 
    def __add__(self, other):
        ret = ExternalCatalog()
        ret.data = np.concatenate([self.data,other.data])
        return ret

 
    def _load(self,filename):
        pass
 
    def match(self,lon,lat,coord='gal',tol=0.1,nnearest=1):
        if coord.lower() == 'cel':
            glon, glat = cel2gal(lon,lat)
        else:
            glon,glat = lon, lat
        return ugali.utils.projector.match(glon,glat,self['glon'],self['glat'],tol,nnearest)
    
