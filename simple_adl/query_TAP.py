#!/usr/bin/env python
"""
Generic python script.
"""
__author__ = "Sidney Mau"

import numpy as np

# Astropy
from astropy import units as u
from astropy.coordinates import SkyCoord

#-------------------------------------------------------------------------------

# Redenning coefficients
R_g = 3.185
R_r = 2.140
R_i = 1.571

def query(service, ra, dec, radius=1.0, gmax=23.5, stars=True, galaxies=False):
    """Return data queried from Rubin TAP
    Parameters
    ----------
    service : TAP service [str]
    ra      : Right Ascension [deg]
    dec     : Declination [deg]
    radius  : radius around (ra, dec) [deg]

    Returns
    -------
    data : numpy recarray of data
    """

    # Define our reference position on the sky and cone radius in arcseconds
    # to use in all following examples
    coord = SkyCoord(ra=ra*u.degree, dec=dec*u.degree, frame='icrs')
    radius = radius * u.deg
    
    # Quality selection and star--galaxy separation adapted from
    # https://github.com/LSSTDESC/DC2-analysis/blob/master/tutorials/object_pandas_stellar_locus.ipynb

    snr_threshold = 25
    mag_err_threshold = 1/snr_threshold
    
    # bright_snr_threshold = 100
    # mag_err_threshold = 1/bright_snr_threshold
    
    safe_max_extended = 1.0
    
    if stars and not galaxies:
        query = f"""
            SELECT
                ra, dec,
                mag_g, mag_r,
                magerr_g, magerr_r,
                mag_g - {R_g} AS mag_corrected_g,
                mag_r - {R_r} AS mag_corrected_r,
                extendedness
            FROM dp01_dc2_catalogs.object
            WHERE CONTAINS(POINT('ICRS', ra, dec), CIRCLE('ICRS', {str(coord.ra.value)}, {str(coord.dec.value)}, {str(radius.value)})) = 1
            AND extendedness < {str(safe_max_extended)}
        """
        # query = f"""
        #     SELECT
        #         ra, dec, mag_g, mag_r, mag_i,
        #         magerr_g, magerr_r, magerr_i,
        #         tract, patch, good, clean,
        #         mag_i_cModel, psFlux_i
        #     FROM dp01_dc2_catalogs.object
        #     WHERE CONTAINS(POINT('ICRS', ra, dec), CIRCLE('ICRS', {str(coord.ra.value)}, {str(coord.dec.value)}, {str(radius.value)})) = 1
        #     AND (  ((mag_i - mag_i_cModel) < 0.03)
        #         OR ( ((mag_i - mag_i_cModel) < 0.1)
        #            AND (mag_i < 22)) )
        # """
    elif not stars and galaxies:
        query = f"""
            SELECT
                ra, dec, mag_g, mag_r, mag_i,
                magerr_g, magerr_r, magerr_i,
                tract, patch, extendedness, good, clean
            FROM dp01_dc2_catalogs.object
            WHERE CONTAINS(POINT('ICRS', ra, dec), CIRCLE('ICRS', {str(coord.ra.value)}, {str(coord.dec.value)}, {str(radius.value)})) = 1
            AND extendedness >= {str(safe_max_extended)}
        """
    else:
        print("stars _and_ galaxies are not supported at the moment; please choose one")

    # For more detailed analysis of results, converting
    # to a pandas dataframe is often very useful
    # results = service.search(query).to_table().to_pandas()
    job = service.submit_job(query)
    job.run()
    job.wait(phases=['COMPLETED', 'ERROR'])
    async_results = job.fetch_result()
    results = async_results.to_table().to_pandas()
    job.delete()

    good_snr = (results['magerr_g'] < mag_err_threshold) & (results['magerr_r'] < mag_err_threshold)
    good_results = results[good_snr]
    
    return good_results

    # query_stars = results[good_snr & star_selection]
    # query_galaxies = results[good_snr & galaxy_selection]

    # if stars and not galaxies:
    #     return query_stars
    # elif not stars and galaxies:
    #     return query_galaxies
    # elif stars and galaxies:
    #     return query_stars, query_galaxies
    # else:
    #     return None

    # return results  # deprecated return call; see if/else statements above


#-------------------------------------------------------------------------------

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    #parser.add_argument('--profile',type=str,required=False,
    #                    help='Profile for data lab query [str]')
    #parser.add_argument('--ra',type=float,required=True,
    #                    help='Right Ascension of target position [deg]')
    #parser.add_argument('--dec',type=float,required=True,
    #                    help='Declination of target position [deg]')
    #parser.add_argument('--radius',type=float,default=1.0,
    #                    help='Radius around target position [deg]')
    #parser.add_argument('--gmax',type=float,default=23.5,
    #                    help='Maximum g-band magnitude [mag]')
    #args = parser.parse_args()
    # data = query(args.profile, args.ra, args.dec, args.radius, args.gmax)

    # Get an instance of the TAP service
    from lsst.rsp import get_tap_service
    service = get_tap_service()
    assert service is not None
    # assert service.baseurl == "https://data.lsst.cloud/api/tap"

    # Define our reference position on the sky and cone radius in arcseconds
    # to use in all following examples
    coord = SkyCoord(ra=62.0*u.degree, dec=-37.0*u.degree, frame='icrs')
    radius = 0.1 * u.deg

    query = f"""
        SELECT
            ra, dec, mag_g, mag_r, mag_i,
            magerr_g, magerr_r, magerr_i,
            tract, patch, extendedness, good, clean
        FROM dp01_dc2_catalogs.object
        WHERE CONTAINS(POINT('ICRS', ra, dec),CIRCLE('ICRS', {str(coord.ra.value)}, {str(coord.dec.value)}, {str(radius.value)})) = 1
    """

    # For more detailed analysis of results, converting
    # to a pandas dataframe is often very useful
    results = service.search(query).to_table().to_pandas()
    # job = service.submit_job(query)
    # job.wait(phases=['COMPLETED', 'ERROR'])
    # results = job.fetch_result()

    # Quality selection and star--galaxy separation adapted from
    # https://github.com/LSSTDESC/DC2-analysis/blob/master/tutorials/object_pandas_stellar_locus.ipynb

    snr_threshold = 25
    mag_err_threshold = 1/snr_threshold
    good_snr = (results['magerr_g'] < mag_err_threshold) & (results['magerr_r'] < mag_err_threshold)

    bright_snr_threshold = 100
    mag_err_threshold = 1/bright_snr_threshold
    bright_snr = (results['magerr_g'] < mag_err_threshold) & (results['magerr_r'] < mag_err_threshold)

    safe_max_extended = 1.0
    star_selection = (results['extendedness'] < safe_max_extended)
    galaxy_selection = (results['extendedness'] >= safe_max_extended)

    stars = results[good_snr & star_selection]
    bright_stars = results[bright_snr & star_selection]
    galaxies = results[good_snr & galaxy_selection]

    print("%d stars (SNR > %.0f)" % (len(stars), snr_threshold))
    print("%d bright stars (SNR > %.0f)" % (len(bright_stars), bright_snr_threshold))
    print("%d galaxies (SNR > %.0f)" % (len(galaxies), snr_threshold))
 
    import pdb;pdb.set_trace()
