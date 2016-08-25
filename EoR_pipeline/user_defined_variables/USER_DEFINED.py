#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
----------------------------------------------------------------------------
Here you will find various python dictionaries of user changeable variables.
You should not need to change anything beyond these variables in order to change
the details of the sky model, telescope and code operation.

However, these only include the most commonly used options, there are more,
these can be added here and then to create_settings (which generates the ini
file for OSKAR) in runOSKAR.py.

See http://www.oerc.ox.ac.uk/~ska/oskar2/OSKAR-Settings.pdf for more information
and an extensive list of options.

Written by Catherine Watkinson (February 2016)
----------------------------------------------------------------------------
"""


class runtime_variables(object):
    """Defines various parameters which control the simulation and imaging"""

    '''
    ------------------------------------------------------------------------
                                Sky model
    ------------------------------------------------------------------------
    'fov_deg' - FoV (degrees) of the passed lightcone;
    'num_pixels_side' - Pixels on lightcone side for each frequency-slice
        (must be integer);
    central RA ('ra0_deg') and DEC ('dec0_deg'), in degrees of the sky region
        modelled;
    'upscale_factor' - Factor by which to scale up the resolution;
    'noise_obs_length' - Time noise integrated over (seconds);
    'noise_bw' - bandwidth (Hz) of observation.
    '''
    user_sky = {
        'fov_deg': 10.0,
        'num_pixels_side': 512,
        'ra0_deg': 0.0,
        'dec0_deg': 89.0,
        'upscale_factor': 2,
        'noise_obs_length': (600.0 * 3600.0),
        'noise_bw': 183.0e3
    }

    '''
    --------------------------------------------------------------------------
                                Observation
    --------------------------------------------------------------------------
    'start_freq':   Central frequency (Hz) of lowest slice in lightcone;
    'freq_inc':     Frequency width of each slice of lightcone (i.e. pixel
                    depth in frequency)
    'start_time':   Start time/date of observation (UTC); length of
                    observation (s or 'h:m:s.z' where z is milliseconds);
    'obs_interval': Observation interval (s) providing the number of time
                    steps in the output data.
    '''
    user_obs = {
        'start_freq': 115.0e6,
        'freq_inc': 0.5e6,
        'start_time': '01-01-2016 00:00:00.000',
        'obs_length': (12.0 * 3600.0),
        'obs_interval': 10.0
    }

    '''
    -------------------------------------------------------------------------
                                Telescope
    -------------------------------------------------------------------------
    'input_dir' - path to telescope model;
    'long_deg' - longitude of telescope centre (degrees);
    'lat_deg' - latitude of telescope centre;
    'pol_mod' - polarisation mode, i.e. whether just simulate scalar
         polarisations or full;
    'norm_beams' - normalise beam centre to 1 if true;
    'beam_duplication' - if true and if all station identical then beam
        response is copied from first;
    'enable_array_pattern' - if true then contribution of array pattern to
        station beam response is evaluated;
    'element_pattern_type' - Type of functional pattern to apply to elements
        (if not specified numerically);
    'dipole_len' - length of dipole (if using);
    'dipole_len_units' - units of quoted dipole length;
    'stat_type' - station type, i.e. aperture array, isotropic or
        Gaussian beam;
    'gauss_fhwm' - FWHM of Gaussian beam (degrees);
    'gauss_ref_freq' - reference frequency for specified FWHM (Hz)
    '''
    user_tel = {
        'input_dir': 'Telescope/HBA_CS.tm',
        'long_deg': 6.868638889,
        'lat_deg': 52.726305556,
        'pol_mod': 'Scalar',
        'norm_beams': 'true',
        'beam_duplication': 'true',
        'enable_array_pattern': 'true',
        'element_pattern_type': 'Dipole',
        'dipole_len': 0.5,
        'dipole_len_units': 'Wavelengths',
        'stat_type': 'Aperture array',
        'gauss_fhwm': 8,
        'gauss_ref_freq': 125000000
    }

    '''
    ------------------------------------------------------------------------
                                Interferometer
    ------------------------------------------------------------------------
    'channel_bw' - Bandwidth of channel (Hz) used to simulate Bandwidth
        smoothing;
    'time_av' - correlator time-average duration (s) to simulate
        time-averaging smearing;
    'uv_min' - minimum baseline length to be evaluated;
    'uv_max' - maximum baseline length to be evaluated;
    'uv_units' - units of uv min and max
    'A_phys' - Physical size of element within a tile
    'spa_den_trans' - frequency (in Hz) of the transition from sparse to dense
        regime of the aperture array. Used for for calculating effective area.
        If you don't have this info, then set to False and an approximation
        will be used.
    'eta_s' - System efficiency
    'N_antenna' - Total number of antennas
    'T_recv' - Receiver temperature

    # NOTE THAT WE MIGHT WANT TO INCLUDE A_SPARSE IN HERE
    (see: http://www.skatelescope.org/uploaded/59513_113_Memo_Nijboer.pdf)
    '''
    user_interferometer = {
        'channel_bw': 183e3,
        'time_av': 10.0,
        'uv_min': 'min',
        'uv_max': 'max',
        'uv_units': 'Wavelengths',
        'A_phys': 1.5625,
        'spa_den_trans': False,
        'eta_s': 1.0,
        'N_antenna': (24.0 * 16.0),
        'T_recv': 160.0
    }

    '''
    ---------------------------------------------------------------------------
                                Imaging
    For more information on these options see:
    https://casa.nrao.edu/docs/TaskRef/clean-task.html
    ---------------------------------------------------------------------------
    'fov_deg' - FoV (degrees) of the CASA image produced from the OSKAR
        visibility set;
    'num_pixels_side' - x, y pixels on of the CASA image produced from the
        OSKAR visibility set (must be integer);
    'niter' - number of iteration over which to clean (i.e. deconvolve the
        dirty image), if set to zero no cleaning is done. Default is 500;
    'gridmode' - Apply corrections for non-coplanar effects during imaging
        using the W-Projection algorithm;
    'wprojplanes' - number of pre-computed w-planes used for the
         W-Projection algorithm
    'uvrange' - range of baselines to include when generating the image
        (https://casa.nrao.edu/Release3.4.0/docs/userman/UserMansu114.html).
        Default is '', i.e. all;
    'weighting' - decides how much weight is given to u-v  grid points to
        allow for different sampling densities. Default is 'natural'.
    '''
    user_image = {
        'fov_deg': 10.0,
        'num_pixels_side': 512,
        'niter': 500,
        'gridmode': 'widefield',
        'wprojplanes': 256,
        'uvrange': '0.01~0.8klambda',
        'weighting': 'uniform'
    }

    '''
    ---------------------------------------------------------------------------
                              ADVANCED USERS ONLY
    Most users should be able to get by without altering the below variables.
    However, more advanced users may well wish to have more control over some
    of them.
    ---------------------------------------------------------------------------

    In order of appearance: if true, clip sourses from below the horizon at
    each time step;
    '''
    user_adv = {
        'horizon_clip':'false',
        'max_chunks':65536,
        'dbl_precision':'true',
        'keep_log':'true'
    }
