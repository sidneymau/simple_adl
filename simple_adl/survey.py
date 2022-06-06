#!/usr/bin/env python
"""
Generic python script.
"""
__author__ = "Sidney Mau"

import glob
import os

import numpy as np
import pandas as pd
import healpy as hp
import scipy.stats
import scipy.interpolate
import scipy.ndimage

import fitsio as fits

from lsst.rsp import get_tap_service

import simple_adl.query_TAP
import simple_adl.projector

#-------------------------------------------------------------------------------

# From https://github.com/DarkEnergySurvey/ugali/blob/master/ugali/utils/healpix.py

def superpixel(subpix, nside_subpix, nside_superpix):
    """
    Return the indices of the super-pixels which contain each of the sub-pixels.
    """
    if nside_subpix==nside_superpix: return subpix
    theta, phi = hp.pix2ang(nside_subpix, subpix)
    return(hp.ang2pix(nside_superpix, theta, phi))


def subpixel(superpix, nside_superpix, nside_subpix):
    """
    Return the indices of sub-pixels (resolution nside_subpix) within
    the super-pixel with (resolution nside_superpix).

    ADW: It would be better to convert to next and do this explicitly
    """
    if nside_superpix==nside_subpix: return superpix
    vec = hp.pix2vec(nside_superpix, superpix)
    radius = np.degrees(2. * hp.max_pixrad(nside_superpix))
    subpix = hp.query_disc(nside_subpix, vec, np.radians(radius))
    pix_for_subpix = superpixel(subpix,nside_subpix,nside_superpix)
    # Might be able to speed up array indexing...
    return(subpix[pix_for_subpix == superpix])

def get_catalog_file(catalog_dir, mc_source_id):
    """
    Inputs:
        catalog_dir = string corresponding to directory containing the stellar catalog infiles
        mc_source_id = integer corresponding the target MC_SOURCE_ID value
    Outputs:
        catalog_infile = string corresponding to filename of stellar catalog containing mc_source_id
    """
    catalog_infiles = sorted(glob.glob(catalog_dir + '/*catalog*.fits'))
    mc_source_id_array = []
    catalog_infile_index_array = []
    for ii, catalog_infile in enumerate(catalog_infiles):
        mc_source_id_min = int(os.path.basename(catalog_infile).split('.')[0].split('mc_source_id_')[-1].split('-')[0])
        mc_source_id_max = int(os.path.basename(catalog_infile).split('.')[0].split('mc_source_id_')[-1].split('-')[1])
        assert (mc_source_id_max > mc_source_id_min) & (mc_source_id_min >= 1), 'Found invalue MC_SOURCE_ID values in filenames'
        mc_source_id_array.append(np.arange(mc_source_id_min, mc_source_id_max + 1))
        catalog_infile_index_array.append(np.tile(ii, 1 + (mc_source_id_max - mc_source_id_min)))

    mc_source_id_array = np.concatenate(mc_source_id_array)
    catalog_infile_index_array = np.concatenate(catalog_infile_index_array)

    assert len(mc_source_id_array) == len(np.unique(mc_source_id_array)), 'Found non-unique MC_SOURCE_ID values in filenames'
    assert np.in1d(mc_source_id, mc_source_id_array), 'Requested MC_SOURCE_ID value not among files'
    mc_source_id_index = np.nonzero(mc_source_id == mc_source_id_array)[0][0] # second [0] added by smau 7/23/18 to fix incompatiable type bug
    return catalog_infiles[catalog_infile_index_array[mc_source_id_index]]

#-------------------------------------------------------------------------------

class Survey():
    """
    Class to handle survey-specific parameters.
    """
    def __init__(self, iterable=(), **kwargs):
        self.__dict__.update(iterable, **kwargs)

        self.mag_1 = self.catalog['mag'].format(self.band_1)
        self.mag_2 = self.catalog['mag'].format(self.band_2)
        self.mag_dered_1 = self.catalog['mag_dered'].format(self.band_1)
        self.mag_dered_2 = self.catalog['mag_dered'].format(self.band_2)
        self.mag_err_1 = self.catalog['mag_err'].format(self.band_1)
        self.mag_err_2 = self.catalog['mag_err'].format(self.band_2)
        self.sim_dir = self.survey['sim_dir']

        self.load_fracdet

    @property
    def load_fracdet(self):
        """
        Load-in the fracdet map if it exists.
        """
        #if self.survey['fracdet']:
        #    print('Reading fracdet map {} ...'.format(self.survey['fracdet']))
        #    fracdet = ugali.utils.healpix.read_map(self.survey['fracdet'])
        #else:
        #    print('No fracdet map specified ...')
        #    fracdet = None
        ##return(fracdet)
        #self.fracdet = fracdet
        # SM: Commenting this out until I have a fracdet map to debug with
        self.fracdet = None

#-------------------------------------------------------------------------------

class Region():
    """
    Class to handle regions.
    """
    def __init__(self, survey, ra, dec):
        self.survey = survey
        self.nside = self.survey.catalog['nside']
        self.fracdet = self.survey.fracdet

        self.ra = ra
        self.dec = dec
        self.proj = simple_adl.projector.Projector(self.ra, self.dec)
        self.pix_center = hp.ang2pix(self.nside, self.ra, self.dec, lonlat=True)

    def load_data(self, stars=True, galaxies=False):
        # SM: to query the equivalent of hp.get_all_neighbors() for nside=32,
        #     choose a radius of 3 deg:
        #>>> np.sqrt((1/np.pi)*8*hp.nside2pixarea(nside=32, degrees=True))
        #2.9238630046262855

        # Get an instance of the TAP service
        service = get_tap_service()
        assert service is not None
        assert service.baseurl == "https://data.lsst.cloud/api/tap"
        data = simple_adl.query_TAP.query(service, self.ra, self.dec, radius=1.0, gmax=self.survey.catalog['mag_max'], stars=stars, galaxies=galaxies)
        # self.data = data
        return data
    
    def load_satellite_sim(self, mc_source_id):
        """
        Load info for injecting satellite sims
        """
        # self.population_file = population_file
        # self.catalog_file = catalog_file
        
        data_array = []
        cat_file = get_catalog_file(self.survey.sim_dir, mc_source_id)
        cat_data = fits.read(cat_file, ext=1)
        # pix = hp.ang2pix(nside,cat_data[basis_1], cat_data[basis_2], lonlat=True)
        # pix_data = cat_data[np.in1d(pix, pix_nside_neighbors)]
        # data_array.append(pix_data)
        # data = np.concatenate(data_array)
        
        return cat_data
    
    def inject_satellite_sim(self, mc_source_id):
        """
        Inject satellite simulation
        """
        
        sim_data = self.load_satellite_sim(mc_source_id)
        # self.sim_data = sim_data
        # self.combined_data = np.concatenate([self.data, self.sim_data])
        real_data = self.load_data(stars=True, galaxies=False)
        sim_data = self.load_satellite_sim(mc_source_id)
        merged_data = real_data[real_data.columns[:-1]].append(pd.DataFrame(sim_data[real_data.columns[:-1]]))  # pd.Dataframe
        
        return merged_data

    def characteristic_density(self, iso_sel):
        """
        Compute the characteristic density of a region
        Convlve the field and find overdensity peaks
        """

        x, y = self.proj.sphereToImage(self.data[self.survey.catalog['basis_1']][iso_sel], self.data[self.survey.catalog['basis_2']][iso_sel]) # Trimmed magnitude range for hotspot finding
        #x_full, y_full = proj.sphereToImage(data[basis_1], data[basis_2]) # If we want to use full magnitude range for significance evaluation
        delta_x = 0.01
        area = delta_x**2
        smoothing = 2. / 60. # Was 3 arcmin
        bins = np.arange(-8., 8. + 1.e-10, delta_x)
        centers = 0.5 * (bins[0: -1] + bins[1:])
        yy, xx = np.meshgrid(centers, centers)
    
        h = np.histogram2d(x, y, bins=[bins, bins])[0]
    
        h_g = scipy.ndimage.filters.gaussian_filter(h, smoothing / delta_x)
    
        #cut_goodcoverage = (data['NEPOCHS_G'][cut_magnitude_threshold] >= 2) & (data['NEPOCHS_R'][cut_magnitude_threshold] >= 2)
        # expect NEPOCHS to be good in DES data
    
        delta_x_coverage = 0.1
        area_coverage = (delta_x_coverage)**2
        bins_coverage = np.arange(-5., 5. + 1.e-10, delta_x_coverage)
        h_coverage = np.histogram2d(x, y, bins=[bins_coverage, bins_coverage])[0]
        #h_goodcoverage = np.histogram2d(x[cut_goodcoverage], y[cut_goodcoverage], bins=[bins_coverage, bins_coverage])[0]
        h_goodcoverage = np.histogram2d(x, y, bins=[bins_coverage, bins_coverage])[0]
    
        n_goodcoverage = h_coverage[h_goodcoverage > 0].flatten()
    
        #characteristic_density = np.mean(n_goodcoverage) / area_coverage # per square degree
        characteristic_density = np.median(n_goodcoverage) / area_coverage # per square degree
        print('Characteristic density = {:0.1f} deg^-2'.format(characteristic_density))
    
        # Use pixels with fracdet ~1.0 to estimate the characteristic density
        if self.fracdet is not None:
            fracdet_zero = np.tile(0., len(self.fracdet))
            cut = (self.fracdet != hp.UNSEEN)
            fracdet_zero[cut] = self.fracdet[cut]
    
            nside_fracdet = hp.npix2nside(len(self.fracdet))
            
            subpix_region_array = []
            for pix in np.unique(hp.ang2pix(self.nside,
                                            self.data[self.survey.catalog['basis_1']][iso_sel],
                                            self.data[self.survey.catalog['basis_2']][iso_sel],
                                            lonlat=True)):
                subpix_region_array.append(subpixel(self.pix_center, self.nside, nside_fracdet))
            subpix_region_array = np.concatenate(subpix_region_array)
    
            # Compute mean fracdet in the region so that this is available as a correction factor
            cut = (self.fracdet[subpix_region_array] != hp.UNSEEN)
            mean_fracdet = np.mean(self.fracdet[subpix_region_array[cut]])
    
            # Correct the characteristic density by the mean fracdet value
            characteristic_density_raw = 1. * characteristic_density
            characteristic_density /= mean_fracdet 
            print('Characteristic density (fracdet corrected) = {:0.1f} deg^-2'.format(characteristic_density))
    
        return(characteristic_density)
    
    def characteristic_density_local(self, iso_sel, x_peak, y_peak, angsep_peak):
        """
        Compute the local characteristic density of a region
        """
    
        #characteristic_density = self.characteristic_density(iso_sel)
        characteristic_density = self.density
    
        x, y = self.proj.sphereToImage(self.data[self.survey.catalog['basis_1']][iso_sel], self.data[self.survey.catalog['basis_2']][iso_sel]) # Trimmed magnitude range for hotspot finding
        #x_full, y_full = proj.sphereToImage(data[basis_1], data[basis_2]) # If we want to use full magnitude range for significance evaluation
    
        # If fracdet map is available, use that information to either compute local density,
        # or in regions of spotty coverage, use the typical density of the region
        if self.fracdet is not None:
            # The following is copied from how it's used in compute_char_density
            fracdet_zero = np.tile(0., len(self.fracdet))
            cut = (self.fracdet != hp.UNSEEN)
            fracdet_zero[cut] = self.fracdet[cut]
    
            nside_fracdet = hp.npix2nside(len(self.fracdet))
            
            subpix_region_array = []
            for pix in np.unique(hp.ang2pix(self.nside,
                                            self.data[self.survey.catalog['basis_1']][iso_sel],
                                            self.data[self.survey.catalog['basis_2']][iso_sel],
                                            lonlat=True)):
                subpix_region_array.append(subpixel(self.pix_center, self.nside, nside_fracdet))
            subpix_region_array = np.concatenate(subpix_region_array)
    
            # Compute mean fracdet in the region so that this is available as a correction factor
            cut = (self.fracdet[subpix_region_array] != hp.UNSEEN)
            mean_fracdet = np.mean(self.fracdet[subpix_region_array[cut]])
    
            subpix_region_array = subpix_region_array[self.fracdet[subpix_region_array] > 0.99]
            subpix = hp.ang2pix(nside_fracdet, 
                                self.data[self.survey.catalog['basis_1']][cut_magnitude_threshold][iso_sel], 
                                self.data[self.survey.catalog['basis_2']][cut_magnitude_threshold][iso_sel],
                                lonlat=True)
    
            # This is where the local computation begins
            ra_peak, dec_peak = self.proj.imageToSphere(x_peak, y_peak)
            subpix_all = hp.query_disc(nside_fracet, hp.ang2vec(ra_peak, dec_peak, lonlat=True), np.radians(0.5))
            subpix_inner = hp.query_disc(nside_fracet, hp.ang2vec(ra_peak, dec_peak, lonlat=True), np.radians(0.3))
            subpix_annulus = subpix_all[~np.in1d(subpix_all, subpix_inner)]
            mean_fracdet = np.mean(fracdet_zero[subpix_annulus])
            print('mean_fracdet {}'.format(mean_fracdet))
            if mean_fracdet < 0.5:
                characteristic_density_local = characteristic_density
                print('characteristic_density_local baseline {}'.format(characteristic_density_local))
            else:
                # Check pixels in annulus with complete coverage
                subpix_annulus_region = np.intersect1d(subpix_region_array, subpix_annulus)
                print('{} percent pixels with complete coverage'.format(float(len(subpix_annulus_region)) / len(subpix_annulus)))
                if (float(len(subpix_annulus_region)) / len(subpix_annulus)) < 0.25:
                    characteristic_density_local = characteristic_density
                    print('characteristic_density_local spotty {}'.format(characteristic_density_local))
                else:
                    characteristic_density_local = float(np.sum(np.in1d(subpix, subpix_annulus_region))) \
                                                   / (hp.nside2pixarea(nside_fracdet, degrees=True) * len(subpix_annulus_region)) # deg^-2
                    print('characteristic_density_local cleaned up {}'.format(characteristic_density_local))
        else:
            # Compute the local characteristic density
            area_field = np.pi * (0.5**2 - 0.3**2)
            n_field = np.sum((angsep_peak > 0.3) & (angsep_peak < 0.5))
            characteristic_density_local = n_field / area_field
    
            # If not good azimuthal coverage, revert
            cut_annulus = (angsep_peak > 0.3) & (angsep_peak < 0.5) 
            #phi = np.degrees(np.arctan2(y_full[cut_annulus] - y_peak, x_full[cut_annulus] - x_peak)) # Use full magnitude range, NOT TESTED!!!
            phi = np.degrees(np.arctan2(y[cut_annulus] - y_peak, x[cut_annulus] - x_peak)) # Impose magnitude threshold
            h = np.histogram(phi, bins=np.linspace(-180., 180., 13))[0]
            if np.sum(h > 0) < 10 or np.sum(h > 0.5 * np.median(h)) < 10:
                #angsep_peak = np.sqrt((x - x_peak)**2 + (y - y_peak)**2)
                characteristic_density_local = characteristic_density
    
        print('Characteristic density local = {:0.1f} deg^-2 = {:0.3f} arcmin^-2'.format(characteristic_density_local, characteristic_density_local / 60.**2))
    
        return(characteristic_density_local)

    def find_peaks(self, iso_sel):
        """
        Convolve field to find characteristic density and peaks within the selected pixel
        """

        #characteristic_density = self.characteristic_density(iso_sel)
        characteristic_density = self.density
    
        x, y = self.proj.sphereToImage(self.data[self.survey.catalog['basis_1']][iso_sel], self.data[self.survey.catalog['basis_2']][iso_sel]) # Trimmed magnitude range for hotspot finding
        #x_full, y_full = proj.sphereToImage(data[basis_1], data[basis_2]) # If we want to use full magnitude range for significance evaluation
        delta_x = 0.01
        area = delta_x**2
        smoothing = 2. / 60. # Was 3 arcmin
        bins = np.arange(-8., 8. + 1.e-10, delta_x)
        #bins = np.arange(-4., 4. + 1.e-10, delta_x) # SM: not sure what to prefer here...
        centers = 0.5 * (bins[0: -1] + bins[1:])
        yy, xx = np.meshgrid(centers, centers)
    
        h = np.histogram2d(x, y, bins=[bins, bins])[0]
        
        h_g = scipy.ndimage.filters.gaussian_filter(h, smoothing / delta_x)
    
        # SM: If we can speed up this block that would be great
        factor_array = np.arange(1., 5., 0.05)
        rara, decdec = self.proj.imageToSphere(xx.flatten(), yy.flatten())
        cutcut = (hp.ang2pix(self.nside, rara, decdec, lonlat=True) == self.pix_center).reshape(xx.shape)
        threshold_density = 5 * characteristic_density * area
        for factor in factor_array:
            # This is reducing the contrast against the background through the arbitrary measurement 'factor'
            # until there are fewer than 10 disconnected peaks
            h_region, n_region = scipy.ndimage.measurements.label((h_g * cutcut) > (area * characteristic_density * factor))
            #print 'factor', factor, n_region, n_region < 10
            if n_region < 10:
                threshold_density = area * characteristic_density * factor
                break
    
        h_region, n_region = scipy.ndimage.measurements.label((h_g * cutcut) > threshold_density)
        #h_region = np.ma.array(h_region, mask=(h_region < 1))
    
        x_peak_array = []
        y_peak_array = []
        angsep_peak_array = []
    
        for index in range(1, n_region + 1): # loop over peaks
            #index_peak = np.argmax(h_g * (h_region == index))
            index_peak = np.ravel_multi_index(scipy.ndimage.maximum_position(input=h_g, labels=h_region, index=index), h_g.shape)
            x_peak, y_peak = xx.flatten()[index_peak], yy.flatten()[index_peak]
            #print index, np.max(h_g * (h_region == index))

            # SM: Could these numbers be useful?
            #index_max = scipy.ndimage.maximum(input=h_g, labels=h_region, index=index)
            #index_stddev = scipy.ndimage.standard_deviation(input=h_g, labels=h_region, index=index)
            #print('max: {}'.format(index_max))
            #print('stddev: {}'.format(index_stddev))
            
            #angsep_peak = np.sqrt((x_full - x_peak)**2 + (y_full - y_peak)**2) # Use full magnitude range, NOT TESTED!!!
            angsep_peak = np.sqrt((x-x_peak)**2 + (y-y_peak)**2)
    
            x_peak_array.append(x_peak)
            y_peak_array.append(y_peak)
            angsep_peak_array.append(angsep_peak)
        
        return x_peak_array, y_peak_array, angsep_peak_array
    
    def fit_aperture(self, iso_sel, x_peak, y_peak, angsep_peak):
        """
        Fit aperture by varing radius and computing the significance
        """

        characteristic_density_local = self.characteristic_density_local(iso_sel, x_peak, y_peak, angsep_peak)
    
        ra_peak_array = []
        dec_peak_array = []
        r_peak_array = []
        sig_peak_array = []
        n_obs_peak_array = []
        n_obs_half_peak_array = []
        n_model_peak_array = []
    
        size_array = np.arange(0.01, 0.3, 0.01)
        sig_array = np.zeros(len(size_array))
        
        size_array_zero = np.concatenate([[0.], size_array])
        area_array = np.pi * (size_array_zero[1:]**2 - size_array_zero[0:-1]**2)
    
        #n_obs_array = np.zeros(len(size_array))
        #n_model_array = np.zeros(len(size_array))
        #for ii in range(0, len(size_array)):
        #    n_obs = np.sum(angsep_peak < size_array[ii])
        #    n_model = characteristic_density_local * (np.pi * size_array[ii]**2)
        #    sig_array[ii] = np.clip(scipy.stats.norm.isf(scipy.stats.poisson.sf(n_obs, n_model)), 0., 37.5) # Clip at 37.5
        #    n_obs_array[ii] = n_obs
        #    n_model_array[ii] = n_model
        n_obs_array = np.array([np.sum(angsep_peak < size) for size in size_array])
        n_model_array = np.array([characteristic_density_local * (np.pi * size**2) for size in size_array])
        sig_array = np.array([np.clip(scipy.stats.norm.isf(scipy.stats.poisson.sf(n_obs, n_model)), 0., 37.5) for (n_obs,n_model) in zip(n_obs_array,n_model_array)])
    
        ra_peak, dec_peak = self.proj.imageToSphere(x_peak, y_peak)
    
        index_peak = np.argmax(sig_array)
        r_peak = size_array[index_peak]
        #if np.max(sig_array) >= 37.5:
        #    r_peak = 0.5
        n_obs_peak = n_obs_array[index_peak]
        n_model_peak = n_model_array[index_peak]
        n_obs_half_peak = np.sum(angsep_peak < (0.5 * r_peak))
    
        # Compile results
        print('Candidate: x_peak: {:12.3f}, y_peak: {:12.3f}, r_peak: {:12.3f}, sig: {:12.3f}, ra_peak: {:12.3f}, dec_peak: {:12.3f}'.format(x_peak, y_peak, r_peak, np.max(sig_array), ra_peak, dec_peak))
        ra_peak_array.append(ra_peak)
        dec_peak_array.append(dec_peak)
        r_peak_array.append(r_peak)
        #sig_peak_array.append(np.max(sig_array))
        sig_peak_array.append(sig_array[index_peak])
        n_obs_peak_array.append(n_obs_peak)
        n_obs_half_peak_array.append(n_obs_half_peak)
        n_model_peak_array.append(n_model_peak)
    
        return ra_peak_array, dec_peak_array, r_peak_array, sig_peak_array, n_obs_peak_array, n_obs_half_peak_array, n_model_peak_array
