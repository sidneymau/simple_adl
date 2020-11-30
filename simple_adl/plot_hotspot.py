#!/usr/bin/env python
"""
Generic python script.
"""
__author__ = "Sid Mau"

# Python libraries
import yaml
import numpy as np
import healpy as hp
import scipy.ndimage

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

import simple_adl.survey
import simple_adl.isochrone
import simple_adl.coordinate_tools
from simple_adl.search import cut_isochrone_path

#-------------------------------------------------------------------------------

props = dict(facecolor='white', edgecolor='black', linewidth=1) #dict(facecolor='white', edgecolor='none', alpha=0.7)

#-------------------------------------------------------------------------------

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--config',type=str,required=False,default='config.yaml',
                        help='Config file [.yaml].')
    parser.add_argument('--outfile',type=str,required=False,default='out.png',
                        help='Output file [.png/.pdf]')
    parser.add_argument('--ra',type=float,required=True,
                        help='Right Ascension of target position [deg].')
    parser.add_argument('--dec',type=float,required=True,
                        help='Declination of target position [deg].')
    parser.add_argument('--mod',type=float,required=False,
                        help='Distance modulus [mag].')
    args = parser.parse_args()

    with open(args.config, 'r') as ymlfile:
        cfg = yaml.load(ymlfile, Loader=yaml.SafeLoader)
        survey = simple_adl.survey.Survey(cfg)

    #---------------------------------------------------------------------------

    region = simple_adl.survey.Region(survey, args.ra, args.dec)
    print('Plot coordinates: (RA, Dec) = ({:0.2f}, {:0.2f})'.format(region.ra, region.dec))

    region.load_data()
    print('Found {} objects'.format(len(region.data)))
    if (len(region.data) == 0):
        print('Ending search.')
        exit()

    data = region.data

    #---------------------------------------------------------------------------

    iso = simple_adl.isochrone.Isochrone(survey=survey.isochrone['survey'],
                              band_1=survey.band_1.lower(),
                              band_2=survey.band_2.lower(),
                              age=12.0, #survey.isochrone['age'],
                              metallicity=0.00010) #survey.isochrone['metallicity']
    if args.mod:
        iso.distance_modulus = args.mod
    #iso_sep = iso.separation(data[mag_1], data[mag_2])
    iso_filter = cut_isochrone_path(data[survey.mag_dered_1], data[survey.mag_dered_2], data[survey.mag_err_1], data[survey.mag_err_2], iso, mag_max=survey.catalog['mag_max'], radius=0.1, return_all=False)
    
    # projection of image
    proj = region.proj
    x, y = proj.sphereToImage(data[survey.catalog['basis_1']], data[survey.catalog['basis_2']])

    # hess
    mag = data[survey.mag_dered_1]
    color = data[survey.mag_dered_1] - data[survey.mag_dered_2]
    
    # Filters
    extension = 0.05
    r0 = 3.0 * extension # 3.0
    r1 = 5.0 * extension # 5.0
    r2 = np.sqrt(r0**2 + r1**2)
    angsep = simple_adl.coordinate_tools.angsep(args.ra, args.dec, data[survey.catalog['basis_1']], data[survey.catalog['basis_2']])
    inner = (angsep < r0)
    outer = ((angsep > r1) & (angsep < r2))
    background = (angsep > r2)

    #---------------------------------------------------------------------------

    fig, axs = plt.subplots(1, 2, figsize=(11, 4))
    fig.subplots_adjust(wspace=0.5)

    #---------------------------------------------------------------------------
    
    # Stellar histogram
    ax = axs[0]
    plt.sca(ax)
    
    bound = 0.5
    steps = 100
    bins = np.linspace(-bound, bound, steps)
    signal = np.histogram2d(x[iso_filter], y[iso_filter], bins=[bins, bins])[0]
    sigma = 0.01 * (0.25 * np.arctan(0.25 * r0 * 60. - 1.5) + 1.3)
    convolution = scipy.ndimage.filters.gaussian_filter(signal, sigma/(bound/steps)).T
    pc = ax.pcolormesh(bins, bins, convolution, cmap='Greys', rasterized=True)
    
    # search kernel
    #x, y = ra_proj, dec_proj
    #delta_x = 0.01
    #area = delta_x**2
    #smoothing = 2. / 60. # Was 3 arcmin
    #bins = np.arange(-8., 8. + 1.e-10, delta_x)
    #centers = 0.5 * (bins[0: -1] + bins[1:])
    #yy, xx = np.meshgrid(centers, centers)
    #h = np.histogram2d(x, y, bins=[bins, bins])[0]
    #h_g = scipy.ndimage.filters.gaussian_filter(h, smoothing / delta_x)
    #pc = ax.pcolormesh(bins, bins, h_g.T, cmap='Greys', rasterized=True)
    
    ax.text(0.05, 0.95, 'Stars', transform=ax.transAxes, verticalalignment='top', bbox=props)
    
    ax.set_xlim(0.5, -0.5)
    ax.set_ylim(-0.5, 0.5)
    ax.set_xlabel(r'$\Delta {\rm RA}$ (deg)')
    ax.set_ylabel(r'$\Delta {\rm Dec}$ (deg)')
    
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='5%', pad=0)
    cb = fig.colorbar(pc, cax=cax)
    #cb.set_label('Counts')

    #---------------------------------------------------------------------------
    
    ## Galactic histogram
    #ax = axs[1]
    #plt.sca(ax)
    #
    #bound = 0.5
    #steps = 100
    #bins = np.linspace(-bound, bound, steps)
    #signal = np.histogram2d(x[galaxies & iso_filter], y[galaxies & iso_filter], bins=[bins, bins])[0]
    #sigma = 0.01 * (0.25 * np.arctan(0.25 * r0 * 60. - 1.5) + 1.3)
    #convolution = scipy.ndimage.filters.gaussian_filter(signal, sigma/(bound/steps)).T
    #pc = ax.pcolormesh(bins, bins, convolution, cmap='Greys', rasterized=True)
    #
    ## search kernel
    ##x, y = ra_proj, dec_proj
    ##delta_x = 0.01
    ##area = delta_x**2
    ##smoothing = 2. / 60. # Was 3 arcmin
    ##bins = np.arange(-8., 8. + 1.e-10, delta_x)
    ##centers = 0.5 * (bins[0: -1] + bins[1:])
    ##yy, xx = np.meshgrid(centers, centers)
    ##h = np.histogram2d(x, y, bins=[bins, bins])[0]
    ##h_g = scipy.ndimage.filters.gaussian_filter(h, smoothing / delta_x)
    ##pc = ax.pcolormesh(bins, bins, h_g.T, cmap='Greys', rasterized=True)
    #
    #ax.text(0.05, 0.95, 'Galaxies', transform=ax.transAxes, verticalalignment='top', bbox=props)
    #
    #ax.set_xlim(0.5, -0.5)
    #ax.set_ylim(-0.5, 0.5)
    #ax.set_xlabel(r'$\Delta \alpha_{2000}$ (deg)')
    #ax.set_ylabel(r'$\Delta \delta_{2000}$ (deg)')
    #
    #divider = make_axes_locatable(ax)
    #cax = divider.append_axes('right', size='5%', pad=0)
    #cb = fig.colorbar(pc, cax=cax)
    ##cb.set_label('Counts')

    #---------------------------------------------------------------------------
    
    # Hess
    ax = axs[1]
    plt.sca(ax)
    
    xbins = np.arange(-0.3, 1.1, 0.1)
    ybins = np.arange(16., 24.0 + 0.5, 0.5)
    foreground = np.histogram2d(color[inner], mag[inner], bins=[xbins, ybins])
    background = np.histogram2d(color[outer], mag[outer], bins=[xbins, ybins])
    fg = foreground[0].T
    bg = background[0].T
    fg_abs = np.absolute(fg)
    bg_abs = np.absolute(bg)
    mask_abs = fg_abs + bg_abs
    mask_abs[mask_abs == 0.] = np.nan # mask common zeroes
    signal = fg - bg
    signal = np.ma.array(signal, mask=np.isnan(mask_abs)) # mask nan
    pc = ax.pcolormesh(xbins, ybins, signal, cmap='viridis', rasterized=True)
    if args.mod:
        plt.plot(iso.color, iso.mag+iso.distance_modulus, lw=1, c='k')
    
    ax.set_xlim(-0.3, 1.0)
    ax.set_ylim(24.0, 16.0)
    ax.set_xlabel(r'$g - r$ (mag)')
    ax.set_ylabel(r'$g$ (mag)')
    
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='5%', pad=0)
    cb = fig.colorbar(pc, cax=cax)
    #cb.set_label('Counts')

    #---------------------------------------------------------------------------
    
    #fig.savefig('./{}.pdf'.format(name), bbox_inches='tight')
    #file_name = 'candidate_{:0.2f}_{:0.2f}'.format(args.ra, args.dec)
    fig.savefig(args.outfile, bbox_inches='tight')
    plt.close(fig)
