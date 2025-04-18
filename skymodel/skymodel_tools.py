import pdb
import time
#import galsim
import numpy as np
from astropy import wcs
from astropy.io import fits


def get_image_size(fov, pixel_scale):  # ensure that image size is a power of 2
    image_size = np.ceil((fov * 60) / (pixel_scale))
    power2 = 0
    if power2:
        exponent = np.floor(np.log2(image_size))
        # choose exponent to make new image smaller or equal to FoV specified in ini file (to avoid accidentally making giant files)
        image_size = np.int(2 ** exponent)
    fov = (image_size * pixel_scale) / 60.0
    return fov, int(image_size)


def get_spectral_size(bw, dnu):  # linear in freq
    crpix3 = np.ceil((bw / dnu) / 2.0)  # round up
    n_chan = int(crpix3 * 2)
    return crpix3, n_chan


def get_spectral_sampling(config, rest_freq, ref_freq, dnu):
    dZ = ((config.getfloat("cosmology", "c") * rest_freq) / (ref_freq ** 2)) * dnu
    return dZ  #  optical velocity interval


def setup_wcs(config, ndim, cosmology=False, nu_axis=False):

    pixel_scale = config.getfloat("skymodel", "pixel_scale")# * galsim.arcsec
    fov = config.getfloat("field", "field_of_view")

    fov, image_size = get_image_size((fov), pixel_scale)# / galsim.arcsec)

    print(pixel_scale,image_size)
    #exit()

    
    dnu = config.getfloat("observation", "channel_width")

    # rest frequency defined only for spectral mode. Eliminate the need for this keyword in continuum

    try:
        rest_freq = config.getfloat("observation", "rest_freq")
    except:
        rest_freq=1.420405752e+9 
    
    base_freq = config.getfloat("observation", "lowest_frequency")
    top_freq = config.getfloat("observation", "highest_frequency")
    bw = top_freq - base_freq
    crpix3, n_chan = get_spectral_size(bw, dnu)

    ref_freq = base_freq + ((top_freq - base_freq) / 2)

    # ra_field = config.get('field', 'field_ra') # now using deg in ini file
    # ra_field_gs = galsim.Angle.from_hms(ra_field)
    # dec_field = config.get('field', 'field_dec')
    # dec_field_gs = galsim.Angle.from_dms(dec_field)

    ra_field_gs = config.getfloat("field", "field_ra")
    dec_field_gs = config.getfloat("field", "field_dec")

    dZ = get_spectral_sampling(config, rest_freq, ref_freq, dnu)

    ref_vel = config.getfloat("cosmology", "c") * (rest_freq / ref_freq - 1)

    msname = config.get("pipeline", "project_name") + ".ms"
    msname = config.get("pipeline", "data_path") + msname
    imagename = msname + ".image"

    w = wcs.WCS(naxis=ndim)

    
    
    if ndim == 4:
    # data cube with frequnecy and stokes dimensions    
        w.wcs.crpix = [
            np.ceil((image_size) / 2.0),  # crpix always integer
            np.ceil((image_size) / 2.0),
            1,
            1,
        ]

        w.wcs.cdelt = [
            -pixel_scale / 3600., #galsim.degrees,
            pixel_scale / 3600.,#galsim.degrees,
            dnu*1.e6,
            1,
        ]

        w.wcs.crval = [ra_field_gs, dec_field_gs, base_freq*1.e6, 1]
        w.wcs.ctype = ["RA---SIN", "DEC--SIN", "FREQ", "STOKES"]
        w.wcs.cunit = ["deg", "deg", "Hz", ""]

        #Initialized for total intensity. IQUV = 1234 should be used for CRVAL4.
        
    elif ndim == 3:

        # begin with a freq-based spectral axis
        w.wcs.crpix = [
            np.ceil((image_size) / 2),
            np.ceil((image_size) / 2),
            crpix3,
        ]
        w.wcs.cdelt = [-pixel_scale / 3600. , pixel_scale / 3600, dZ]
        w.wcs.crval = [ra_field_gs, dec_field_gs, ref_vel]
        w.wcs.ctype = [
            "RA---SIN",
            "DEC--SIN",
            "VOPT-F2W",
        ]  # see https://www.astro.rug.nl/software/kapteyn/spectralbackground.html (FELO_HEL behaves as VOPT_F2W)
        w.wcs.cunit = ["deg", "deg", "m/s"]
        w.wcs.restfrq = rest_freq

    elif ndim == 2:
        w.wcs.crpix = [np.ceil((image_size) / 2), np.ceil((image_size) / 2)]
        w.wcs.cdelt = [-pixel_scale / 3600, pixel_scale / 3600.]
        w.wcs.crval = [ra_field_gs, dec_field_gs]
        w.wcs.ctype = ["RA---SIN", "DEC--SIN"]
        w.wcs.cunit = ["deg", "deg"]

    return w
