#!/usr/bin/env python
"""
Generic python script.
"""
__author__ = "Sidney Mau"

import sys
import os
import glob
import yaml
import numpy as np
import healpy as hp
import scipy.interpolate

import fitsio as fits

import simple_adl.survey
import simple_adl.isochrone
from simple_adl.coordinate_tools import distanceModulusToDistance, angsep
from IPython.core.debugger import set_trace

#-------------------------------------------------------------------------------

def cut_isochrone_path(g, r, g_err, r_err, isochrone, mag_max, radius=0.1, return_all=False):
    """
    Cut to identify objects within isochrone cookie-cutter.
    """
    if np.all(isochrone.stage == 'Main'):
        # Dotter case
        index_transition = len(isochrone.stage)
    else:
        # Other cases
        index_transition = np.nonzero(isochrone.stage >= isochrone.hb_stage)[0][0] + 1    

    mag_1_rgb = isochrone.mag_1[0: index_transition] + isochrone.distance_modulus
    mag_2_rgb = isochrone.mag_2[0: index_transition] + isochrone.distance_modulus
    
    mag_1_rgb = mag_1_rgb[::-1]
    mag_2_rgb = mag_2_rgb[::-1]

    # Cut one way...
    f_isochrone = scipy.interpolate.interp1d(mag_2_rgb, mag_1_rgb - mag_2_rgb, bounds_error=False, fill_value = 999.)
    color_diff = np.fabs((g - r) - f_isochrone(r))
    cut_2 = (color_diff < np.sqrt(0.1**2 + r_err**2 + g_err**2))

     # ...and now the other
    f_isochrone = scipy.interpolate.interp1d(mag_1_rgb, mag_1_rgb - mag_2_rgb, bounds_error=False, fill_value = 999.)
    color_diff = np.fabs((g - r) - f_isochrone(g))
    cut_1 = (color_diff < np.sqrt(0.1**2 + r_err**2 + g_err**2))

    cut = np.logical_or(cut_1, cut_2)

    ## Cut for horizontal branch
    #mag_1_hb = isochrone.mag_1[isochrone.stage == isochrone.hb_stage][1:] + isochrone.distance_modulus
    #mag_2_hb = isochrone.mag_2[isochrone.stage == isochrone.hb_stage][1:] + isochrone.distance_modulus

    #f_isochrone = scipy.interpolate.interp1d(mag_2_hb, mag_1_hb - mag_2_hb, bounds_error=False, fill_value = 999.)
    #color_diff = np.fabs((g - r) - f_isochrone(r))
    #cut_4 = (color_diff < np.sqrt(0.1**2 + r_err**2 + g_err**2))

    #f_isochrone = scipy.interpolate.interp1d(mag_1_hb, mag_1_hb - mag_2_hb, bounds_error=False, fill_value = 999.)
    #color_diff = np.fabs((g - r) - f_isochrone(g))
    #cut_3 = (color_diff < np.sqrt(0.1**2 + r_err**2 + g_err**2))
    #
    #cut_hb = np.logical_or(cut_3, cut_4)

    #cut = np.logical_or(cut, cut_hb)

    #mag_bins = np.arange(17., 24.1, 0.1)
    mag_bins = np.arange(17., mag_max+0.1, 0.1)
    mag_centers = 0.5 * (mag_bins[1:] + mag_bins[0:-1])
    magerr = np.tile(0., len(mag_centers))
    for ii in range(0, len(mag_bins) - 1):
        cut_mag_bin = (g > mag_bins[ii]) & (g < mag_bins[ii + 1])
        magerr[ii] = np.median(np.sqrt(0.1**2 + r_err[cut_mag_bin]**2 + g_err[cut_mag_bin]**2))

    if return_all:
        return cut, mag_centers[f_isochrone(mag_centers) < 100], (f_isochrone(mag_centers) + magerr)[f_isochrone(mag_centers) < 100], (f_isochrone(mag_centers) - magerr)[f_isochrone(mag_centers) < 100]
    else:
        return cut

def write_output(results_dir, nside, pix_nside_select, best_ra_peak, best_dec_peak, best_r_peak, best_distance_modulus, 
                n_obs_peak, n_obs_half_peak, n_model_peak, 
                best_sig_peak, mc_source_id, mode, outfile):
    #writer = open(outfile, 'a') # append if exists
    #for ii in range(0, len(sig_peak_array)):
    #    # SIG, RA, DEC, MODULUS, r, n_obs, n_model, mc_source_id
    #    writer.write('{:10.3f}, {:10.3f}, {:10.3f}, {:10.3f}, {:10.3f}, {:10.3f}, {:10.3f}, {:10.3f}, {:10.3f}\n'.format(sig_peak_array[ii], 
    #                                                                                                                     ra_peak_array[ii], 
    #                                                                                                                     dec_peak_array[ii], 
    #                                                                                                                     distance_modulus_array[ii], 
    #                                                                                                                     r_peak_array[ii],
    #                                                                                                                     n_obs_peak_array[ii],
    #                                                                                                                     n_obs_half_peak_array[ii],
    #                                                                                                                     n_model_peak_array[ii],
    #                                                                                                                     mc_source_id_array[ii]))
    #data = [tuple(row) for row in np.stack([sig_peak_array, ra_peak_array, dec_peak_array, distance_modulus_array, r_peak_array, n_obs_peak_array, n_obs_half_peak_array, n_model_peak_array, mc_source_id_array], axis=-1)]
    #arr = np.array(data, dtype=[('SIG', float), ('RA', float), ('DEC', float), ('MODULUS', float), ('R', float), ('N_OBS', float), ('N_OBS_HALF', float ('N_MODEL', float), ('MC_SOURCE_ID', int))])
    #np.save(outfile, arr)
    #f = open(os.path.join(results_dir,outfile), 'ab')
    #np.savetxt(f, arr, delimiter=',')
    #f.close()
    
    #saving only the best
    data = np.array([best_sig_peak, best_ra_peak, best_dec_peak, best_distance_modulus, best_r_peak, n_obs_peak, n_obs_half_peak, n_model_peak, mc_source_id]) 
    f = open(os.path.join(results_dir,outfile), 'ab')
    if os.stat(os.path.join(results_dir,outfile)).st_size == 0:
        np.savetxt(f, [data], fmt="%.2f", delimiter=',', header="SIG,RA,DEC,MODULUS,R,N_OBS,N_OBS_HALF,N_MODEL,MC_SOURCE_ID", comments="")
    else:
         np.savetxt(f, [data], fmt="%.2f", delimiter=',')
    f.close()

def search_by_distance(survey, region, distance_modulus, iso_sel, extension=None):
    """
    Idea: 
    Send a data extension that goes to faint magnitudes, e.g., g < 24.
    Use the whole region to identify hotspots using a slightly brighter 
    magnitude threshold, e.g., g < 23, so not susceptible to variations 
    in depth. Then compute the local field density using a small annulus 
    around each individual hotspot, e.g., radius 0.3 to 0.5 deg.
    """

    #print('Distance = {:0.1f} kpc (m-M = {:0.1f})'.format(distanceModulusToDistance(distance_modulus), distance_modulus))

    #print('{} objects left after isochrone cut...'.format(len(data)))

    if (len(region.data[iso_sel]) == 0):
        return [], [], [], [], [], [], [], []

    ra_peak_array = []
    dec_peak_array = []
    r_peak_array = []
    sig_peak_array = []
    distance_modulus_array = []
    n_obs_peak_array = []
    n_obs_half_peak_array = []
    n_model_peak_array = []

    region.density = region.characteristic_density(iso_sel)
    x_peak_array, y_peak_array, angsep_peak_array = region.find_peaks(iso_sel)

    for x_peak, y_peak, angsep_peak in zip(x_peak_array, y_peak_array, angsep_peak_array):
        # Aperture fitting
        print('Fitting aperture to hotspot...')
        ra_peaks, dec_peaks, r_peaks, sig_peaks, n_obs_peaks, n_obs_half_peaks, n_model_peaks = region.fit_aperture(iso_sel, x_peak, y_peak, angsep_peak, extension)
        
        ra_peak_array.append(ra_peaks)
        dec_peak_array.append(dec_peaks)
        r_peak_array.append(r_peaks)
        sig_peak_array.append(sig_peaks)
        distance_modulus_array.append(distance_modulus*np.ones(len(ra_peaks)))
        n_obs_peak_array.append(n_obs_peaks)
        n_obs_half_peak_array.append(n_obs_half_peaks)
        n_model_peak_array.append(n_model_peaks)
        
    try:
        ra_peak_array = np.concatenate(ra_peak_array)
        dec_peak_array = np.concatenate(dec_peak_array)
        r_peak_array = np.concatenate(r_peak_array)
        sig_peak_array = np.concatenate(sig_peak_array)
        distance_modulus_array = np.concatenate(distance_modulus_array)
        n_obs_peak_array = np.concatenate(n_obs_peak_array)
        n_obs_half_peak_array = np.concatenate(n_obs_half_peak_array)
        n_model_peak_array = np.concatenate(n_model_peak_array)
    except ValueError:
        print('No arrays to concatenate')

    return ra_peak_array, dec_peak_array, r_peak_array, sig_peak_array, distance_modulus_array, n_obs_peak_array, n_obs_half_peak_array, n_model_peak_array

#-------------------------------------------------------------------------------

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--config',type=str,required=False,default='config.yaml',
                        help='Config file [.yaml]')
    parser.add_argument('--outfile',type=str,required=False,default='out.csv',
                        help='Output file [.csv]')
    parser.add_argument('--ra',type=float,required=False,
                        help='Right Ascension of target position [deg]')
    parser.add_argument('--dec',type=float,required=False,
                        help='Declination of target position [deg]')
    parser.add_argument('--nside',type=int,required=False,
                        help='Nside of HEALPix pixelization [int]')
    parser.add_argument('--ipix',type=int,required=False,
                        help='Pixel index [int]')
    parser.add_argument('--mod',type=float,required=False,
                        help='Fixed distance modulus search value')
    parser.add_argument('--mcid',type=int,required=False,
                        help='MC SOURCE ID of satellite to inject')
    args = parser.parse_args()

    if args.ra and args.dec:
        if args.nside or args.ipix:
            parser.error('Please specify either (ra, dec) or (nside, ipix).')
        ra, dec = args.ra, args.dec
    elif args.nside and args.ipix:
        if args.ra and args.dec:
            parser.error('Please specify either (ra, dec) or (nside, ipix).')
        ra, dec = hp.pix2ang(nside=args.nside, ipix=args.ipix, lonlat=True)
    else:
        parser.error('Please specify either (ra, dec) or (nside, ipix).')

    with open(args.config, 'r') as ymlfile:
        cfg = yaml.load(ymlfile, Loader=yaml.SafeLoader)
        survey = simple_adl.survey.Survey(cfg)

    #---------------------------------------------------------------------------

    region = simple_adl.survey.Region(survey, ra, dec)
    print('Search coordinates: (RA, Dec) = ({:0.2f}, {:0.2f})'.format(region.ra, region.dec))
    print('Search healpixel: {} (nside = {})'.format(region.pix_center, region.nside))

    #---------------------------------------------------------------------------

    if args.mcid:
        data = region.inject_satellite_sim(args.mcid)
    else:
        data = region.load_data(stars=True, galaxies=False)
    
    region.data = data
    
    # region.load_data(stars=True, galaxies=False)
    print('Found {} objects'.format(len(region.data)))
    if (len(region.data) == 0):
        print('Ending search.')
        exit()
    
    #---------------------------------------------------------------------------

    if args.mod:
        # Search at fixed distance modulus
        distance_modulus_search = args.mod
        iso_search = simple_adl.isochrone.Isochrone(survey=survey.isochrone['survey'],
                                                    band_1=survey.band_1.lower(),
                                                    band_2=survey.band_2.lower(),
                                                    age=12.0, #survey.isochrone['age'],
                                                    metallicity=0.00010, #survey.isochrone['metallicity'],
                                                    distance_modulus=distance_modulus_search)
        
        iso_selection = cut_isochrone_path(region.data[survey.mag_dered_1], 
                                           region.data[survey.mag_dered_2],
                                           region.data[survey.mag_err_1],
                                           region.data[survey.mag_err_2],
                                           iso_search,
                                           survey.catalog['mag_max'],
                                           radius=0.1)
                         
        results = search_by_distance(survey, region, distance_modulus_search, iso_selection)
    else:
        # Scan in distance moduli    
        distance_modulus_search_array = np.arange(16., survey.catalog['mag_max'], 0.5)

        #results = [np.empty((1,9)) for distance_modulus in distance_modulus_search_array]

        iso_search_array = [simple_adl.isochrone.Isochrone(survey=survey.isochrone['survey'],
                                                           band_1=survey.band_1.lower(),
                                                           band_2=survey.band_2.lower(),
                                                           age=12.0, #survey.isochrone['age'],
                                                           metallicity=0.00010, #survey.isochrone['metallicity'],
                                                           distance_modulus=distance_modulus)
                            for distance_modulus in distance_modulus_search_array]

        iso_selection_array = [cut_isochrone_path(region.data[survey.mag_dered_1], 
                                                  region.data[survey.mag_dered_2],
                                                  region.data[survey.mag_err_1],
                                                  region.data[survey.mag_err_2],
                                                  iso,
                                                  survey.catalog['mag_max'],
                                                  radius=0.1)
                               for iso in iso_search_array]

        #data_array = [region.data[iso_sel] for iso_sel in iso_selection_array]

        #ra_peak_array, dec_peak_array, r_peak_array, sig_peak_array, n_obs_array, n_obs_half_array, n_model_array = [search_by_distance(survey, region, data) for data in data_array]

        results = [search_by_distance(survey, region, distance_modulus, iso_sel) for (distance_modulus,iso_sel) in zip(distance_modulus_search_array,iso_selection_array)]
                         
    ra_peak_array, dec_peak_array, r_peak_array, sig_peak_array, distance_modulus_array, n_obs_peak_array, n_obs_half_peak_array, n_model_peak_array = np.asarray(results)

    # ra_peak_array = np.concatenate(ra_peak_array)
    # dec_peak_array = np.concatenate(dec_peak_array)
    # r_peak_array = np.concatenate(r_peak_array)
    # sig_peak_array = np.concatenate(sig_peak_array)
    # distance_modulus_array = np.concatenate(distance_modulus_array)
    # n_obs_peak_array = np.concatenate(n_obs_peak_array)
    # n_obs_half_peak_array = np.concatenate(n_obs_half_peak_array)
    # n_model_peak_array = np.concatenate(n_model_peak_array)

    if args.mcid:
        mc_source_id_array = np.full_like(distance_modulus_array, args.mcid)
    else:
        mc_source_id_array = np.zeros(len(distance_modulus_array))
    
    # Sort peaks according to significance
    index_sort = np.argsort(sig_peak_array)[::-1]
    ra_peak_array = ra_peak_array[index_sort]
    dec_peak_array = dec_peak_array[index_sort]
    r_peak_array = r_peak_array[index_sort]
    sig_peak_array = sig_peak_array[index_sort]
    distance_modulus_array = distance_modulus_array[index_sort]
    n_obs_peak_array = n_obs_peak_array[index_sort]
    n_obs_half_peak_array = n_obs_half_peak_array[index_sort]
    n_model_peak_array = n_model_peak_array[index_sort]
    mc_source_id_array = mc_source_id_array[index_sort]
    
    # Collect overlapping peaks
    for ii in range(0, len(sig_peak_array)):
        if sig_peak_array[ii] < 0:
            continue
        sep = angsep(ra_peak_array[ii], dec_peak_array[ii], ra_peak_array, dec_peak_array)
        sig_peak_array[(sep < r_peak_array[ii]) & (np.arange(len(sig_peak_array)) > ii)] = -1.
        #sig_peak_array[(sep < 0.5) & (np.arange(len(sig_peak_array)) > ii)] = -1. # 0.5 deg radius
    
    # Prune the list of peaks
    ra_peak_array = ra_peak_array[sig_peak_array > 0.]
    dec_peak_array = dec_peak_array[sig_peak_array > 0.]
    r_peak_array = r_peak_array[sig_peak_array > 0.]
    distance_modulus_array = distance_modulus_array[sig_peak_array > 0.]
    n_obs_peak_array = n_obs_peak_array[sig_peak_array > 0.]
    n_obs_half_peak_array = n_obs_half_peak_array[sig_peak_array > 0.]
    n_model_peak_array = n_model_peak_array[sig_peak_array > 0.]
    mc_source_id_array = mc_source_id_array[sig_peak_array > 0.]
    sig_peak_array = sig_peak_array[sig_peak_array > 0.] # Update the sig_peak_array last!
    
    for ii in range(0, len(sig_peak_array)):
        print('{:0.2f} sigma; (RA, Dec, d) = ({:0.2f}, {:0.2f}); r = {:0.2f} deg; d = {:0.1f}, mu = {:0.2f} mag), mc_source_id: {:0.2f}'.format(sig_peak_array[ii], 
                         ra_peak_array[ii], 
                         dec_peak_array[ii], 
                         r_peak_array[ii],
                         distanceModulusToDistance(distance_modulus_array[ii]),
                         distance_modulus_array[ii],
                         mc_source_id_array[ii]))
    
    # Write output
    if (len(sig_peak_array) > 0):
        write_output(survey.output['results_dir'], survey.catalog['nside'], region.pix_center, ra_peak_array, dec_peak_array,
                     r_peak_array, distance_modulus_array, 
                     n_obs_peak_array, n_obs_half_peak_array, n_model_peak_array, 
                     sig_peak_array, mc_source_id_array, 0, args.outfile)
    else:
        print('No significant hotspots found.')
        nan_array = [np.nan]
        write_output(survey.output['results_dir'], survey.catalog['nside'], region.pix_center,
                     nan_array, nan_array, nan_array, nan_array, 
                     nan_array, nan_array, nan_array, nan_array,
                     [mc_source_id], 0, args.outfile)
