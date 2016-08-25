#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
This is an updated version of sim.py created by Catherine Watkinson (Feb 2016).

This takes the sky model filename as an arguement allowing cosmological signal,
foregrounds and noise to be run in parallel.

User defined variables have also been centralised in USER_DEFINED_VARIABLES.py
"""
from __future__ import (absolute_import, print_function, division,
                        unicode_literals)
import os
import sys
import numpy as np
import pyfits
from PIL import Image
from subprocess import call, Popen, PIPE
import re
from math import sin, asin, cos, degrees, radians, log, pi
from user_defined_variables.USER_DEFINED import runtime_variables as rv


def setup_data_dirs(dir_list):
    """Utility function to create a set of directories."""
    for dir_name in dir_list:
        if not os.path.isdir(dir_name):
            os.makedirs(dir_name)


def channel_to_freq(idx):
    """Convert channel index to frequency, in hz"""
    freq0 = rv.user_obs['start_freq'] #cw 115.0e6
    inc = rv.user_obs['freq_inc'] #cw 0.5e6
    return freq0 + (idx * inc)


def fov_to_cellsize(fov_deg, size):
    """Convert image FoV and size in pixels to cellsize in arcseconds.

    Args:
        fov_deg (float): Image FoV, in degrees
        size (int): Image size in one dimension in pixels

    Returns:
        tuple (float, float): Image cellsize, in arcseconds,
            lm increment per pixel.
    """
    inc = (2 * sin(radians(fov_deg) * 0.5)) / size
    return degrees(asin(inc)) * 3600.0, inc


def image_lm_grid(size, fov_rad):
    """Returns the l, m coordinates for each pixel in the image."""
    _, inc = fov_to_cellsize(degrees(fov_rad), size)
    lm = np.arange(-size // 2, size // 2) * inc
    l, m = np.meshgrid(-lm, lm)
    return l, m


def lm_to_apparent_ra_dec(l, m, ra0_rad, dec0_rad):
    """Convert image direction cosines (l, m) to equatorial coordinates.

    Args:
        l (array, float): l direction cosines to convert.
        m (array, float): m direction cosines to convert.
        ra0_rad (float):
            Reference Right Ascension (image phase centre), in radians.
        dec0_rad (float):
            Reference Declination (image phase centre), in radians.

    Returns:
        tuple (array, array): Equatorial coordinates for each input l, m
        direction, in radians.
    """
    n = (1.0 - l**2 - m**2)**0.5
    sin_dec0 = sin(dec0_rad)
    cos_dec0 = cos(dec0_rad)
    dec_rad = np.arcsin(n * sin_dec0 + m * cos_dec0)
    ra_rad = ra0_rad + np.arctan2(l, n * cos_dec0 - m * sin_dec0)
    return ra_rad, dec_rad


def image_coords(size, fov_rad, ra0_rad, dec0_rad):
    """Obtain image equatorial coordinates for each pixel in an image."""
    l, m = image_lm_grid(size, fov_rad)
    return lm_to_apparent_ra_dec(l, m, ra0_rad, dec0_rad)


def cube_to_osm(cube_filename, fov_deg, channelID, ra0_deg, dec0_deg,
                osm_filename, upscale_factor, save_fits=False):
    """Convert image cube ot an oskar sky model file"""
    ra0_rad = radians(ra0_deg)
    dec0_rad = radians(dec0_deg)

    # Open the FITS image data cube
    # -------------------------------------------------------------------------
    hdulist = pyfits.open(cube_filename)
    cube = hdulist[0].data
    assert(channelID >= 0)
    assert(channelID < cube.shape[0])
    hdulist.close()

    # Extract image for given channel ID.
    # -------------------------------------------------------------------------
    im0 = np.array(cube[channelID, :, :], dtype=np.float32)

    # Rescale the raw image (im0) to obtain im1
    # -------------------------------------------------------------------------
    old_size = im0.shape[0]
    freqHz = channel_to_freq(channelID)
    size = int(np.ceil(upscale_factor * old_size))
    #im1 = imresize(im0, (size,size), interp='bicubic', mode='F')
    im_ = Image.frombuffer('F', im0.shape, im0.tostring(), 'raw', 'F', 0, 1)
    im_ = im_.resize((size,size), Image.NEAREST)
    im1 = np.array(im_.getdata(), dtype=np.float32).reshape((size, size))

    # Convert the rescaled image to an OSM file.
    # -------------------------------------------------------------------------
    # Obtain list of image coordinates for the rescaled map pixels
    ra, dec = image_coords(size, fov_deg * np.pi/180.0, ra0_rad, dec0_rad)
    ra  *= 180.0/np.pi
    dec *= 180.0/np.pi

    num_pixels = int(size*size)
    sky = np.zeros((num_pixels, 3))
    sky[:, 0] = ra.reshape((num_pixels,))
    sky[:, 1] = dec.reshape((num_pixels,))
    sky[:, 2] = im1.reshape((num_pixels,))

    # Remove sources with amplitude == 0.0
    sky = sky[sky[:, 2] != 0.0, :]

    # Convert to Jy/pixel
    old_cellsize, _ = fov_to_cellsize(fov_deg, old_size)
    new_cellsize, _ = fov_to_cellsize(fov_deg, size)
    # new pixel size is smaller than the old pixel size so the pixel area ratio
    # will be < 1 for all cases of up-scaling.
    kB = 1.3806488e-23
    c0 = 299792458.0
    pixel_area = radians(new_cellsize / 3600.0)**2
    # Convert from brightness temperature in K to Jy/Pixel
    # http://www.iram.fr/IRAMFR/IS/IS2002/html_1/node187.html
    sky[:, 2] *= 2.0 * kB * 1.0e26 * pixel_area * (freqHz/c0)**2

    if os.path.exists(osm_filename):
        os.remove(osm_filename)

    if np.__version__ == '1.4.1':
        np.savetxt(osm_filename, sky, fmt='%.10e, %.10e, %.10e')
    else:
        np.savetxt(osm_filename, sky, fmt='%.10e, %.10e, %.10e',
                   header = str(
                       'Channel = %i\n'
                       'Frequency = %e MHz\n'
                       'Number of sources = %i\n'
                       'Cellsize = %f arcsec (pixel separation at centre)\n'
                       'RA0 = %f\n'
                       'Dec0 = %f\n'
                       % (channelID, freqHz, len(sky), new_cellsize, ra0_deg,
                          dec0_deg)))

    # Save FITS maps of the selected channel slice
    if save_fits:
        if os.path.basename(cube_filename)[0] == 'f':
            rootname = os.path.basename(cube_filename)[:-5]
        elif os.path.basename(cube_filename)[0] == 'c':
            rootname = os.path.basename(cube_filename)[:-7]
        else:
            print(cube_filename)
            raise RuntimeError('Unknown cube filename')
        #outpath = 'os.path.dirname(cube_filename)'
        #rootname = os.path.basename(cube_filename) #cw added so below doesn't end up calling Models/Models/filename
        outpath = 'Models' #cw this allows the cube under analysis to be located anywhere

        img = '%s/IMG_%s_K_%03i.fits' % (outpath, rootname, channelID)
        #img = 'models/image_slice_%03i.fits' % (channelID)
        if os.path.exists(img): os.remove(img)
        im0 = np.reshape(im0, (1, im0.shape[0], im0.shape[1]))
        #im0[0,256,256] = 100
        hdu = pyfits.PrimaryHDU(im0)
        hdulist = pyfits.HDUList([hdu])
        hdr = hdulist[0].header
        hdr['BUNIT']  = ('K', 'Brightness unit')
        hdr['CTYPE1'] = 'RA---SIN'
        hdr['CRVAL1'] = ra0_deg
        hdr['CDELT1'] = -old_cellsize / 3600
        hdr['CRPIX1'] = im0.shape[1] // 2 + 1 # WARNING! Assumes even image dims
        hdr['CTYPE2'] = 'DEC--SIN'
        hdr['CRVAL2'] = dec0_deg
        hdr['CDELT2'] = old_cellsize / 3600.0
        hdr['CRPIX2'] = im0.shape[1] // 2 + 1 # WARNING! Assumes even image dims
        hdr['CTYPE3'] = 'FREQ'
        hdr['CRVAL3'] = freqHz
        hdr['CDELT3'] = 1
        hdr['CRPIX3'] = 1
        hdulist.writeto(img)
        hdulist.close()

        #rescaled_img = 'models/rescaled_image_slice_%03i.fits' % (channelID)
        rescaled_img = ('%s/IMG_%s_K_rescaled_%03i.fits' %
                        (outpath, rootname, channelID))
        if os.path.exists(rescaled_img): os.remove(rescaled_img)
        im1 = np.reshape(im1, (1, im1.shape[0], im1.shape[1]))
        hdu = pyfits.PrimaryHDU(im1)
        hdulist = pyfits.HDUList([hdu])
        hdr = hdulist[0].header
        hdr['BUNIT']  = ('K', 'Brightness unit')
        hdr['CTYPE1'] = 'RA---SIN'
        hdr['CRVAL1'] = ra0_deg
        hdr['CDELT1'] = -new_cellsize / 3600
        hdr['CRPIX1'] = im1.shape[1] // 2 + 1 # WARNING! Assumes even image dims
        hdr['CTYPE2'] = 'DEC--SIN'
        hdr['CRVAL2'] = dec0_deg
        hdr['CDELT2'] = new_cellsize / 3600
        hdr['CRPIX2'] = im1.shape[1] // 2 + 1 # WARNING! Assumes even image dims
        hdr['CTYPE3'] = 'FREQ'
        hdr['CRVAL3'] = freqHz
        hdr['CDELT3'] = 1
        hdr['CRPIX3'] = 1
        hdulist.writeto(rescaled_img)
        hdulist.close()

        # Convert to Jy/Beam
        uv_max = 3500 # m
        FWHM = 1.22 * (c0 / freqHz) / uv_max  # Radians
        beam_area = (pi * FWHM**2) / (4.0 * log(2))
        im2 = im1 * 2.0 * kB * 1.0e26 * beam_area * (freqHz / c0)**2
        clean_component_img = ('%s/IMG_%s_Jy_per_beam_%03i.fits' %
                               (outpath,rootname, channelID))
        if os.path.exists(clean_component_img): os.remove(clean_component_img)
        hdu = pyfits.PrimaryHDU(im2)
        hdulist = pyfits.HDUList([hdu])
        hdr = hdulist[0].header
        hdr['BUNIT']  = ('Jy/beam', 'Brightness unit')
        hdr['CTYPE1'] = 'RA---SIN'
        hdr['CRVAL1'] = ra0_deg
        hdr['CDELT1'] = -new_cellsize / 3600.0
        hdr['CRPIX1'] = im1.shape[1] / 2 + 1 # WARNING! Assumes even image dims
        hdr['CTYPE2'] = 'DEC--SIN'
        hdr['CRVAL2'] = dec0_deg
        hdr['CDELT2'] = new_cellsize / 3600.0
        hdr['CRPIX2'] = im1.shape[1] / 2 + 1 # WARNING! Assumes even image dims
        hdr['CTYPE3'] = 'FREQ'
        hdr['CRVAL3'] = freqHz
        hdr['CDELT3'] = 1
        hdr['CRPIX3'] = 1
        hdulist.writeto(clean_component_img)
        hdulist.close()

        # Convert to Jy/Pixel
        im3 = im1 * 2.0 * kB * 1.0e26 * pixel_area * (freqHz/c0)**2
        sky_model_img = ('%s/IMG_%s_Jy_per_pixel_%03i.fits' %
                         (outpath, rootname, channelID))
        if os.path.exists(sky_model_img): os.remove(sky_model_img)
        hdu = pyfits.PrimaryHDU(im3)
        hdulist = pyfits.HDUList([hdu])
        hdr = hdulist[0].header
        hdr['BUNIT']  = ('Jy/pixel', 'Brightness unit')
        hdr['CTYPE1'] = 'RA---SIN'
        hdr['CRVAL1'] = ra0_deg
        hdr['CDELT1'] = -new_cellsize / 3600.0
        hdr['CRPIX1'] = im1.shape[1] / 2 + 1 # WARNING! Assumes even image dims
        hdr['CTYPE2'] = 'DEC--SIN'
        hdr['CRVAL2'] = dec0_deg
        hdr['CDELT2'] = new_cellsize / 3600.0
        hdr['CRPIX2'] = im1.shape[1] / 2 + 1 # WARNING! Assumes even image dims
        hdr['CTYPE3'] = 'FREQ'
        hdr['CRVAL3'] = freqHz
        hdr['CDELT3'] = 1
        hdr['CRPIX3'] = 1
        hdulist.writeto(sky_model_img)
        hdulist.close()

    print('  Input file           = %s' % cube_filename)
    print('  Frequency            = %.4f MHz' % (freqHz / 1e6))
    print('  Cube image size      =', im0.shape[1])
    print('  No. sources          = %i [%i %s]' %
          (len(sky), size**2-len(sky), 'removed (==0.0)'))
    print('  Output sky model     =', osm_filename)
    print('  Writing FITS images  =', save_fits)
    if save_fits:
        print('  Beam area            =', beam_area)
        print('  Pixel area           =', pixel_area)

    return len(sky)


def set_setting(ini, key, value):
    """Set a setting into an OSKAR ini file."""
    call(["oskar_settings_set", "-q", ini, key, str(value)])


def run_interferometer(ini, verbose=True):
    """Run the OSKAR interferometry simulator."""
    if verbose:
        call(["oskar_sim_interferometer", ini])
    else:
        call(["oskar_sim_interferometer", "-q", ini])


def dict_to_settings(settings_dict, filename):
    """Convert a python dictionary to an OSKAR ini settings file."""
    for group in sorted(settings_dict.keys()):
        for key in sorted(settings_dict[group].keys()):
            key_ = group+key
            value_ = settings_dict[group][key]
            set_setting(filename, key_, value_)


def require_oskar_version(version):
    """Check for a specified version of OSKAR."""
    try:
        call(['oskar_sim_interferometer', '--version'], stdout=PIPE,
             stderr=PIPE)
    except OSError:
        raise Exception('OSKAR not found. Check your PATH settings.')
    proc = Popen('oskar_sim_interferometer --version', stdout=PIPE, shell=True)
    (out,err) = proc.communicate()
    out = out.strip('\n\r')
    ver = re.split('\.|-', out)
    ver[0] = int(ver[0])
    ver[1] = int(ver[1])
    ver[2] = int(ver[2])
    sVersion = '%i.%i.%i' % (version[0], version[1], version[2])
    if len(version) == 4: sVersion='%s-%s' % (sVersion, version[3])
    sVer = '%i.%i.%i' % (ver[0], ver[1], ver[2])
    if len(ver) == 4: sVer='%s-%s' % (sVer, ver[3])
    failMsg = ("ERROR: This script requires OSKAR %s [found version %s]." %
               (sVersion, out))
    if len(ver)!= len(version):
        print(failMsg)
    for i in range(0, len(ver)):
        if ver[i] != version[i]:
            print(failMsg)
    return ver


def evaluate_noise_rms_Jy(freqHz, bw, obs_length):
    """
    Evaluate baseline noise RMS, in Jy.

    see: http://www.skatelescope.org/uploaded/59513_113_Memo_Nijboer.pdf

    Args:
        freqHz (float): Frequency, in Hz
        bw (float): Bandwidth, in Hz
        obs_length (float): Observation length, in seconds

    Returns:
        Noise RMS in Jy
    """


    c0 = 299792458. #cw  m/s
    kB = 1.3806488e-23 #cw m^2 kg s^{-2} K^{-1}
    lambda_ = c0 / freqHz

    # Values from Peter Dewdney's spreadsheet.
    A_sparse = lambda_**2 / 3.0 #cw made 3 -> 3.0 out of paranoia
    A_physical = rv.user_interferometer['A_phys'] #cw 1.5625          # Element physical size within HBA tile (5mx5m tile / 16 antenna)
    eta_s = rv.user_interferometer['eta_s'] #cw 1.0                  # System efficiency
    #cw replaced station size with Nantenna:
    N_antenna = rv.user_interferometer['N_antenna'] #cw 24.0*16.0         # Total number of antennas (core station, single patch = 24 tiles of 16 antenna)
    T_recv = rv.user_interferometer['T_recv'] #cw 160.0                 # http://arxiv.org/pdf/1104.1577.pdf (140-180K)
    #T_recv = 180.0                # http://arxiv.org/pdf/1104.1577.pdf (140-180K)

    #print A_physical
    #print eta_s
    #print station_size
    #print T_recv

    #cw I need to adjust this for SKA as we know the frequency at which transition from dense to sparse occurs, make use of 'spa_den_trans':False parameter
    # Get effective area per dipole antenna, A_eff
    A_eff = np.minimum(A_sparse, A_physical)

    # Get system temperature.
    T_sky = 60.0 * lambda_**2.55 # "Standard" empirical T_sky
    T_sys = T_sky + T_recv

    # Get RMS noise per baseline for single polarisation.
    SEFD_station = ((2.0 * kB * T_sys) / (eta_s * A_eff * N_antenna)) #cw replaced station size with Nantenna
    SEFD_station *= 1e26 # Convert to Jy
    sigma_pq = SEFD_station / (2.0 * bw * obs_length)**0.5

    return sigma_pq


def create_settings(freqHz, sky_model_name, ra0_deg, dec0_deg, ms_name,
                    start_time, obs_length, num_time_steps, uvmax, add_noise,
                    noise_rms_Jy, noise_seed):
    """Create dictionary of simulation settings.
    cw This has been adjusted to pull all variables from the centralised
       User defined variable file in the ./User_defined_variables/ folder.
    """
    s = dict()
    s['simulator/'] = {
        'max_sources_per_chunk':rv.user_adv['max_chunks'], #cw 65536,
        'double_precision':rv.user_adv['dbl_precision'], #cw 'true',
        'keep_log_file':rv.user_adv['keep_log'] #cw 'true'
    }
    s['telescope/'] = {
        #'input_directory':'models/HBA_CS_layout_only.tm',
        'input_directory': rv.user_tel['input_dir'],#cw 'Telescope/HBA_CS.tm',
        'longitude_deg':rv.user_tel['long_deg'], #cw 6.868638889,
        'latitude_deg':rv.user_tel['lat_deg'], #cw 52.726305556,
        'pol_mode':rv.user_tel['pol_mod'], #cw 'Scalar',
        'normalise_beams_at_phase_centre':rv.user_tel['norm_beams'], #cw 'true',
        'allow_station_beam_duplication':rv.user_tel['beam_duplication'], #'true',
        'aperture_array/array_pattern/enable':rv.user_tel['enable_array_pattern'], #'true',
        'aperture_array/element_pattern/functional_type':rv.user_tel['element_pattern_type'], #'Dipole',
        'aperture_array/element_pattern/dipole_length':rv.user_tel['dipole_len'], #0.5,
        'aperture_array/element_pattern/dipole_length_units':rv.user_tel['dipole_len_units'], #'Wavelengths',
        #'gaussian_beam/fwhm_deg':rv.user_tel['gauss_fhwm'], #8,
        #'gaussian_beam/ref_freq_hz':rv.user_tel['gauss_ref_freq'], #125000000
        'station_type':rv.user_tel['stat_type'], #'Aperture array', 'Isotropic','Gaussian beam'
    }
    #print s['telescope/']['input_directory']
    #print s['telescope/']['longitude_deg']
    #print s['telescope/']['latitude_deg']
    #print s['telescope/']['pol_mode']
    #print s['telescope/']['normalise_beams_at_phase_centre']
    #print s['telescope/']['allow_station_beam_duplication']
    #print s['telescope/']['aperture_array/array_pattern/enable']
    #print s['telescope/']['aperture_array/element_pattern/functional_type']
    #print s['telescope/']['aperture_array/element_pattern/dipole_length']
    #print s['telescope/']['aperture_array/element_pattern/dipole_length_units']
    #print s['telescope/']['station_type']

    s['observation/'] = {
        'phase_centre_ra_deg':ra0_deg,
        'phase_centre_dec_deg':dec0_deg,
        'start_frequency_hz':freqHz,
        'num_channels':1,
        'start_time_utc':start_time,
        'length':obs_length,
        'num_time_steps':num_time_steps
    }
    if not add_noise:
        s['sky/'] = {
            'oskar_sky_model/file':sky_model_name,
            'advanced/apply_horizon_clip':rv.user_adv['horizon_clip'] #cw 'false'
        }
    s['interferometer/'] = {
        'channel_bandwidth_hz':rv.user_interferometer['channel_bw'], #cw 183e3,
        'time_average_sec':rv.user_interferometer['time_av'], #cw 10.0,
        'ms_filename':ms_name,
        'uv_filter_min': rv.user_interferometer['uv_min'], #cw 'min', 10.0,
        'uv_filter_max': rv.user_interferometer['uv_max'], #cw 'max', 800.0,
        'uv_filter_units': rv.user_interferometer['uv_units'] #cw 'Wavelengths',
    }
    #print s['interferometer/']['channel_bandwidth_hz']
    #print s['interferometer/']['time_average_sec']
    #print s['interferometer/']['uv_filter_min']
    #print s['interferometer/']['uv_filter_max']
    #print s['interferometer/']['uv_filter_units']
    if add_noise == True:
        noise = {
            #'oskar_vis_filename':ms_name[:-3]+'.vis',
            'noise/enable':add_noise,
            'noise/seed':noise_seed,
            'noise/values':'RMS flux density',
            'noise/values/rms':'Range',
            'noise/values/rms/start':noise_rms_Jy,
            'noise/values/rms/end':noise_rms_Jy,
            'noise/freq':'Range',
            'noise/freq/number':1,
            'noise/freq/start':freqHz,
            'noise/freq/inc':0
        }
        s['interferometer/'].update(noise)
    return s


if __name__ == '__main__':
    '''cw if len(sys.argv) != 5:
             print 'Usage: ./sim.py [channel] [generate data] [run sim] [mode 1=fg,2=cs_sf_tapered,3=noise]'

    start_channel = int(sys.argv[1])
    end_channel = start_channel
    generate_data = int(sys.argv[2])
    run_sim = int(sys.argv[3])
    mode = int(sys.argv[4])
    '''

    #cw start added block:
    '''
    Have enabled the filename to be passed to the program as it's called
    (allows the code to be run in parallel).
    Have also made the code run so that by default it runs in full mode.
    To stop the generation of data and running of simulation, then the user
    simply changes 1 1 to 0 0 in the $options variable of the submit file.
    '''
    if len(sys.argv) < 5:
        if len(sys.argv) != 7:
            print('Usage (cw updated): python ./runOSKAR.py '
                  '[sky model] [start channel] [end channel] '
                  '[mode 1=fg,2=cs_sf_tapered,3=noise] '
                  '[generate data (optional)] [run sim (optional)]')
            sys.exit(1)

    sky_root_name = str(sys.argv[1])
    start_channel = int(sys.argv[2])
    end_channel = int(sys.argv[3])
    mode = int(sys.argv[4])

    if len(sys.argv) < 7:
        generate_data = 1
        run_sim = 1
    else:
        generate_data = int(sys.argv[5])
        run_sim = int(sys.argv[6])
    print(sky_root_name, start_channel, end_channel, mode, generate_data, run_sim)

    # Sky model data.
    #cw these variables are now defined in User_defined_variables/USER_DEFINED.py
    fov_deg = rv.user_sky['fov_deg'] #10.0
    num_pixels_side = rv.user_sky['num_pixels_side'] #512
    ra0_deg = rv.user_sky['ra0_deg'] #0.0
    dec0_deg = rv.user_sky['dec0_deg'] #89.0
    upscale_factor = rv.user_sky['upscale_factor'] #2.0
    noise_obs_length = rv.user_sky['noise_obs_length'] #600.0 * 3600.0
    noise_bw = rv.user_sky['noise_bw'] #183.0e3

    print(fov_deg)
    print(num_pixels_side)
    print(ra0_deg)
    print(dec0_deg)
    print(upscale_factor)
    print(noise_obs_length)
    print(noise_bw) #'''
    #cw end added block

    if mode == 1:
        #cw (managed by [sky model argv]) sky_root_name = 'rough_LOSsigma_0_010_fg_zeromean'
        add_noise = False
    elif mode == 2:
        #cw (managed by [sky model argv])sky_root_name = 'cs_zeromean'
        add_noise = False
    elif mode == 3:
        #cw (managed by [sky model argv])sky_root_name = 'noise'
        add_noise = True
    else:
        print('ERROR: Invalid mode option.')
        os.exit(1)

    cube_filename = '%s.fits' % (sky_root_name)
    print(cube_filename)

    #cw Added creation of Models folder as it's use was hardcoded, but it's existence not
    setup_data_dirs(['Vis','Ini','Models'])
    rootname  = os.path.basename(sky_root_name) #cw addition to allow flexibility of sky model location
    if generate_data and mode != 3:
        for channelID in range(start_channel, end_channel + 1):
            # Generate sky model from slice of FITS cube.
            osm_filename  = 'Models/%s_%03i.osm' % (rootname, channelID) #cw replaced sky_root_name with rootname
            cube_to_osm(cube_filename, fov_deg, channelID, ra0_deg, dec0_deg,
                        osm_filename, upscale_factor, True)

    if run_sim:
        require_oskar_version([2, 6, 0,'trunk'])
        for channelID in range(start_channel, end_channel + 1):

            # Set up parameters.
            freqHz = channel_to_freq(channelID)
            ms_name = 'Vis/%s_%03i.ms' % (rootname, channelID) #cw replaced sky_root_name with rootname
            ini = 'Ini/%s_%03i.ini' % (rootname, channelID) #cw replaced sky_root_name with rootname
            osm_filename  = 'Models/%s_%03i.osm' % (rootname, channelID) #cw replaced sky_root_name with rootname
            start_time = rv.user_obs['start_time'] #cw '01-01-2016 00:00:00.000'
            obs_length = rv.user_obs['obs_length'] #cw 12.0 * 3600.0
            obs_interval = rv.user_obs['obs_interval'] #cw 10.0
            num_time_steps = int(np.ceil(obs_length / obs_interval))
            uvmax = (upscale_factor * num_pixels_side) / \
                    (2.0 * fov_deg * np.pi/180.0)
            noise_rms_Jy = evaluate_noise_rms_Jy(freqHz, noise_bw, noise_obs_length)

            # cw print start_time, obs_length, obs_interval

            # Create settings file.
            noise_seed = channelID + 100
            '''
            #cw
            The following creates a dictionary of the user defined settings
            and then generates an ini file that OSKAR may use.
            '''
            s = create_settings(freqHz, osm_filename, ra0_deg, dec0_deg,
                                ms_name, start_time, obs_length, num_time_steps,
                                uvmax, add_noise, noise_rms_Jy, noise_seed)
            dict_to_settings(s, ini)

            print('Running simulation for channel %i, freq = %.4f MHz' %
                  (channelID, freqHz / 1.0e6))
            run_interferometer(ini)
