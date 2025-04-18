"""
Script to convert a T-RECS catalogue into a continuum sky model FITS file.

Usage:
python skymodel.py example.ini

Author: Anna Bonaldi, Ian Harrison


"""


import logging
import multiprocessing
import os
import sys
import time
import numpy as np
import pickle
from astropy import units as uns
from astropy.coordinates import SkyCoord
from astropy.cosmology import LambdaCDM
from astropy.io import fits as astfits
from astropy.table import Table
import fitsio
from fitsio import FITS, FITSHDR
from numpy.core.defchararray import add as stradd
from numpy.core.defchararray import multiply as strmultiply
from multiprocessing import Manager
from multiprocessing import Process, current_process #AA
import skymodel.skymodel_tools as tools
from skymodel.continuum_morphology import make_img
from skymodel.skymodel_tools import setup_wcs



arcsectorad = (1.0 * uns.arcsec).to(uns.rad).value
degtoarcsec = (1.0 * uns.deg).to(uns.arcsec).value
degtorad= (1.0 * uns.deg).to(uns.rad).value

# Add outputs procuced by seleval parallel processes to a single final output
# .fits
# _P.fits
# _Q.fits
# _U.fits

def coadd(filename,n_cores,tag):

    print('coadding ',filename+tag )

    os.system("rm {0}".format(filename+tag+".fits"))
    #initialise with the first file
    data,h = fitsio.read(filename+"_1"+tag+".fits ",header=True)

    for i in range(1,n_cores):
        process_tag="_"+str(i+1)
        #print(process_tag)
        #read map 2
        data2,h2 = fitsio.read(filename+process_tag+tag+".fits ",header=True)

        data=data+data2
        #print(np.min(data),np.max(data))

    fitsf = FITS(filename+tag+".fits", "rw")
    fitsf.write(data, header=h)

    fitsf.close()
    
    
    # delete old files
    for i in range(0,n_cores):
        process_tag="_"+str(i+1)
        os.system("rm {0}".format(filename+process_tag+tag+".fits"))


# Combine outputs produced by seleval parallel processes to a single final output
# _maxflux.fits
# _z.fits 
def combine(filename,n_cores):

    print('combining ',filename+"_maxflux",filename+"_z" )

    os.system("rm {0}".format(filename+"_maxflux.fits"))
    os.system("rm {0}".format(filename+"_z.fits"))

    #initialise with the first file
    data_f,h_f = fitsio.read(filename+"_1_maxflux.fits ",header=True)
    data_z,h_z = fitsio.read(filename+"_1_z.fits ",header=True)

    for i in range(1,n_cores):
        process_tag="_"+str(i+1)
        #print(process_tag)
        #read map 2
        data_f2,h = fitsio.read(filename+process_tag+"_maxflux.fits ",header=True)
        data_z2,h = fitsio.read(filename+process_tag+"_z.fits ",header=True)

        
        data_z[data_f2 > data_f] = data_z2[data_f2 > data_f]
        data_f[data_f2 > data_f] = data_f2[data_f2 > data_f]
        
        
    

    fitsf_f = FITS(filename+"_maxflux.fits", "rw")
    fitsf_f.write(data_f, header=h_f)


    
    fitsf_z = FITS(filename+"_z.fits", "rw")
    fitsf_z.write(data_z, header=h_z)

    fitsf_f.close()
    fitsf_z.close()
    
    # delete old files
    for i in range(0,n_cores):
        process_tag="_"+str(i+1)
        os.system("rm {0}".format(filename+process_tag+"_maxflux.fits"))
        os.system("rm {0}".format(filename+process_tag+"_z.fits"))




    

def log_result(result):
    
    global cat
    (i, atlas_source, flux, unresolved,skipped_sources) = result
    cat["id"][i] = i
    cat["Atlas_source"][i] = atlas_source
    cat["New_flux"][i] = flux
    cat["Unresolved"][i] = unresolved
    skipped_sources=skipped_sources

def add_source_continuum(
    i,
    cat_gal,
    nobj,
    w_twod,
    config,
    pixel_scale_str,
    psf_maj_arcsec,
    arr_dims,
    all_gals_fname,
    base_freq,
    freqs,
    polarization,
    skipped_sources    
):
    
    
    
    logging.info(
        "..........Adding source {0} of {1} to skymodel..........".format(i + 1, nobj)
    )

    total_sources_added=0 #initialise this variable, if everything goes well make_img returns 1
    
    # how to ask process name:
    #process = multiprocessing.current_process()
    #pid=process.pid
    
    #x, y = w_twod.wcs_world2pix(
    #    cat_gal["RA"],
    #    cat_gal["DEC"],
    #    1,
    #)

    
    #x = float(x)
    #y = float(y)

    x=cat_gal["xs"]
    y=cat_gal["ys"]
    #print('vecchio',x,y)
    #print('nuovo',cat_gal["xs"],cat_gal["ys"])
    
    logging.info("RA, Dec: %f %f ", cat_gal["RA"], cat_gal["DEC"])
    logging.info("PA, flux: %f %f ", cat_gal["PA"], cat_gal["Total_flux"])
    logging.info("class:  %f", cat_gal["RadioClass"])
    logging.info("x, y,: %f %f ", x, y)
    
    logging.info("Continuum size from cat:  %f", cat_gal["Maj"])

    # get the postage for the source
    # it can be AGN from library, Gaussian lobe and Gaussian core, Sersic of simple Gaussian
    sub_img, atlas_source, unresolved, flux = make_img(
        config,
        cat_gal["Total_flux"],
        cat_gal["spectral_index"],
        base_freq,
        freqs,
        cat_gal["Maj"],
        cat_gal["Min"],
        cat_gal["PA"],
        float(pixel_scale_str),
        psf_maj_arcsec,
        cat_gal["RadioClass"],
        cat_gal["corefrac"],
        cat_gal["ranid"],
    )

    
    sub_img_shape = sub_img.shape
    sub_img_size = sub_img.shape[1]
    logging.info("postage stamp size %f", sub_img_size)

    # works out the bounds for the postage stamp in the FoV image
    l_bounds = np.array([0, y, x]) - np.array(
        [0, (sub_img.shape[1] / 2), (sub_img.shape[2] / 2)]
    )
    u_bounds = np.array([len(freqs), y, x]) + np.array(
        [
            0,
            sub_img_shape[1] - (sub_img.shape[1] / 2),
            sub_img_shape[2] - (sub_img.shape[2] / 2),
        ]
    )

    logging.info("Lower bounds, upper bounds: %s, %s", l_bounds, u_bounds)

    l_bounds = np.floor(l_bounds).astype(np.int)
    u_bounds = np.floor(u_bounds).astype(np.int)

    logging.info(
        "Lower bounds, upper bounds, int: %s, %s",
        l_bounds,
        u_bounds,
    )
    logging.info("Subcube shape: %s", sub_img.shape)

    # add it to the large cube
    img3 = sub_img
    blc0 = l_bounds[0]
    blc1 = l_bounds[1]
    blc2 = l_bounds[2]
    trc0 = u_bounds[0] - 1
    trc1 = u_bounds[1] - 1
    trc2 = u_bounds[2] - 1

    # the top bounds are all -1 the true values, since the large cube is added using Fortran indexing
    trcs = np.array([trc0, trc1, trc2])
    blcs = np.array([blc0, blc1, blc2])

    # pixels from the image to top coordinates of field
    top_excess = arr_dims - (
        trcs + 1
    )  
    bottom_excess = blcs
    excess = np.hstack((bottom_excess, top_excess))
    # initialise indicator to say whether the galaxy overlaps an edge
    overlap = False  
    # the galaxy is clipped if it overlaps an edge
    for coord in excess:
        if coord < 0:
            overlap = True
            logging.info("Subcube is overlapping the edge: cropping to fit")            
            break

    if overlap:
        start_list = np.copy(bottom_excess)
        end_list = np.copy(top_excess)
        np.putmask(start_list, bottom_excess < 0, (-bottom_excess))
        np.putmask(start_list, bottom_excess >= 0, 0)
        start0, start1, start2 = start_list
        np.putmask(end_list, top_excess >= 0, img3.shape)
        end0, end1, end2 = end_list
        img3 = img3[start0:end0, start1:end1, start2:end2]
        np.putmask(blcs, bottom_excess < 0, 0)
        np.putmask(trcs, top_excess < 0, arr_dims - 1)
        blc0, blc1, blc2 = blcs
        trc0, trc1, trc2 = trcs

    logging.info(
        "BLC, TRC: %f %f %f, %f %f %f ",
        blc0,
        blc1,
        blc2,
        trc0,
        trc1,
        trc2,
    )

    # to prevent crashes in case the object is cropped so much one of the dimensions is zero
    if ((trc1+1 <=blc1) or (trc2+1<=blc2)):
        logging.info("Skipping this source as outside of field")
        skipped_sources=skipped_sources+1
        return (i, atlas_source, flux, unresolved,skipped_sources)

    
    # write the info for this object to files
    fitsf = FITS(all_gals_fname+".fits", "rw")
    fitsf_f = FITS(all_gals_fname + "_maxflux.fits", "rw")
    fitsf_z = FITS(all_gals_fname + "_z.fits", "rw")
    #testpola

              
        # the redshift map contains the redshift of the brightest source on the LoS. 
        # this is judged by looking at the dummy map _maxflux and comparing if with the postage stamp. 
        # the z map is updated only where the postage stamp is brighter than what recorder in _maxflux

        # the flux is I for Polarization == False and P for POlarization == True
        # this is guaranteed by polafract which is initialised as 1. for polarization == false
        
        # read the recorded values for flux and redshift at the postage location
    img3_pola=img3*cat_gal['polafrac']
    
    flux_old = fitsf_f[0][0:1, blc1 : trc1 + 1, blc2 : trc2 + 1]
    z_old = fitsf_z[0][0:1, blc1 : trc1 + 1, blc2 : trc2 + 1]

        # initialise the new arrays
    flux_new = z_old * 0.0
    flux_new[0] = img3[0]*cat_gal['polafrac']  # at the lowest frequency
        # for total intensity only, polafrac=1 and so the flux is total intensity flux
        
    zvalue = cat_gal["z"]
    img_z = z_old
    img_f = flux_old

        # update only where postage brighter than record
    img_z[flux_new > flux_old] = zvalue
    img_f[flux_new > flux_old] = flux_new[flux_new > flux_old]

    fitsf_f[0].write(img_f, 0, blc1, blc2, 0, trc1, trc2)
    fitsf_z[0].write(img_z, 0, blc1, blc2, 0, trc1, trc2)

        # adding the source to the total map
        # if running in parallel, this step can go out of synch
    region = fitsf[0][blc0 : trc0 + 1, blc1 : trc1 + 1, blc2 : trc2 + 1]
    img3 += region
    fitsf[0].write(img3, blc0, blc1, blc2, trc0, trc1, trc2)

    fitsf.close()
    fitsf_f.close()
    fitsf_z.close()

    if (polarization == True):
        
        fitsf_p = FITS(all_gals_fname + "_P.fits", "rw")
        fitsf_q = FITS(all_gals_fname + "_Q.fits", "rw")
        fitsf_u = FITS(all_gals_fname + "_U.fits", "rw")
                    
        img3_q=img3_pola*np.cos(cat_gal['EVPA']/ 180. * np.pi)
        img3_u=img3_pola*np.sin(cat_gal['EVPA']/ 180. * np.pi)
            
            
        region = fitsf_p[0][blc0 : trc0 + 1, blc1 : trc1 + 1, blc2 : trc2 + 1]
            
        img3_pola += region
        
        fitsf_p[0].write(img3_pola, blc0, blc1, blc2, trc0, trc1, trc2)

        region = fitsf_q[0][blc0 : trc0 + 1, blc1 : trc1 + 1, blc2 : trc2 + 1]
                            
        img3_q += region
        fitsf_q[0].write(img3_q, blc0, blc1, blc2, trc0, trc1, trc2)



        region = fitsf_u[0][blc0 : trc0 + 1, blc1 : trc1 + 1, blc2 : trc2 + 1]
        img3_u += region
        fitsf_u[0].write(img3_u, blc0, blc1, blc2, trc0, trc1, trc2)

        fitsf_p.close()
        fitsf_q.close()
        fitsf_u.close()

                        

        
    logging.info("")

    return (i, atlas_source, flux, unresolved)

### run the continuum sky model
def runSkyModel(config,process,total_cores):
    """Simulate a sky model from a T-RECS catalogue.

    Parameters
    ----------
    config : configparser
        ConfigParser configuration containing necessary sections.

    """

    # set up parallel processing name tags 

    process_tag=""
    #process_add=0
    if total_cores >1:
        process_tag="_"+str(process)

        time.sleep(60*process) #stagger processes to avoid too much reading and writing on disk at the same time  
    tstart = time.time()    
    # Set up logging
    logfilename = "logs/%s.log" % config.get("field", "fits_prefix")+process_tag
    os.system("rm %s" % logfilename)
    log = logfilename
    logging.basicConfig(
        filename=log,
        level=logging.DEBUG,
        format="%(asctime)s %(message)s",
        datefmt="%d/%m/%Y %H:%M:%S",
    )

    
    logging.info("Beginning simulation")


    # read keywords from the config file
    n_cores = int(config.getfloat("pipeline", "n_cores"))
    logging.info("Running with %d cores", n_cores)
    doplot = config.getboolean("pipeline", "doplot")

    # set up cosmology
    H = config.getfloat("cosmology", "H")
    M = config.getfloat("cosmology", "M")
    L = config.getfloat("cosmology", "L")
    cosmo = LambdaCDM(H0=H, Om0=M, Ode0=L)

    mother_seed = int(config.get("pipeline", "mother_seed"))

    logging.info("mother_seed %s",mother_seed)
    data_path_large_files = (
        config.get("pipeline", "data_path_large_files")
        + config.get("pipeline", "project_name")
        + "/"
    )
    data_path = (
        config.get("pipeline", "base_dir")
        + config.get("pipeline", "data_path")
        + config.get("pipeline", "project_name")
        + "/"
    )

    # create output path if not existing
    if not os.path.exists(data_path):
        os.system('mkdir -p '+ data_path)

   
    # set image properties
    psf_maj_arcsec = config.getfloat("observation", "simple_psf_maj")
    psf_maj = psf_maj_arcsec * arcsectorad
    #print(psf_maj_arcsec,psf_maj,psf_maj_arcsec*arcsectorad)
    psf_min = config.getfloat("observation", "simple_psf_min") * arcsectorad
    psf_pa = config.getfloat("observation", "simple_psf_pa") * degtorad

    pixel_scale = config.getfloat("skymodel", "pixel_scale")
    pixel_scale_str = str(pixel_scale).split()[0]
    fov = config.getfloat("field", "field_of_view")
    logging.info("FoV from ini file, arcmin: %f", fov)

    # set sky coordinates
    ra_field_gs = config.getfloat("field", "field_ra")

    # convert to range +/- 180 to enable cutoffs later
    if ra_field_gs > 180.0:
        ra_field_gs -= 360.0
    dec_field_gs = config.getfloat("field", "field_dec")
    global cat
    
    # set spectral properties for continuum
    logfreq = False
    if config.getboolean("observation", "dologfreq") == True:
        logfreq = True


    base_freq = config.getfloat("observation", "lowest_frequency")  # *1.e6 #Hz
    base_freqname = config.get("observation", "lowest_frequency")
    top_freq = config.getfloat("observation", "highest_frequency")  # *1.e6 #Hz
    top_freqname = config.get("observation", "highest_frequency")
    fov, image_size = tools.get_image_size(fov, pixel_scale)

    logging.info("Image_size, pixels: %f", image_size)
    logging.info("Final base_freq, Hz: %f", base_freq)
    logging.info("Final top_freq, Hz: %f", top_freq)

    dnu = config.getfloat("observation", "channel_width")

    # this is the definittion of linearly-spaced frequencies
    if (logfreq == False):
        nfreqs = int((top_freq - base_freq) / dnu) + 1
        freqs = np.zeros(nfreqs).astype(np.float32)
        freq = base_freq
        for ff in range(nfreqs):
            freqs[ff] = freq
            freq = freq + dnu


    # this is the definition of logarithmically-spaced frequencies
    if (logfreq == True):
        nfreqs=0
        freq=base_freq

        while freq <= top_freq:
            freq=freq*dnu
            nfreqs=nfreqs+1

        freqs = np.zeros(nfreqs).astype(np.float32)
        freq = base_freq
        for ff in range(nfreqs):
            freqs[ff] = freq
            freq = freq * dnu

        #print(freqs)
        #print(nfreqs)

                    
    n_chan = nfreqs
    arr_dims = np.array([n_chan, image_size, image_size]).astype(np.int)

    logging.info(
        "Final array dimensions: %d %d %d" % (n_chan, image_size, image_size)
    )
    logging.info(
        "Final array size, elements: %.3e" % (n_chan * image_size * image_size)
    )
    logging.info(
        "Final array size, bytes: %.3e" % (n_chan * image_size * image_size * 4)
    )

    # calling the appropriate wcs for continuum.
    # based on header_4D plus extra keywords
    
    w_twod = setup_wcs(config, ndim=2, cosmology=cosmo)
    w_fourd = setup_wcs(config, ndim=4)
    
    header_fourd = w_fourd.to_header()
    header_fourd["BUNIT"] = "JY/PIXEL"


    bmaj=psf_maj / degtorad  #they are in radians, get to degree  
    bmin = psf_min / degtorad #they are in radians, get to degree  
    bpa = psf_pa / degtorad #they are in radians, get to degree  

    header_fourd["BMAJ"] = bmaj
    header_fourd["BMIN"] = bmin
    header_fourd["BPA"] = bpa

    # initialse empty cubes
    all_gals_fname = (
        data_path_large_files + config.get("field", "fits_prefix")+ process_tag
    )

    # CHOICE1: never overwrite files
    #if os.path.exists(all_gals_fname+".fits"):
    #    print ('**** message from pipeline: '+all_gals_fname+'.fits already exists') 
    #    print ('**** message from pipeline: not running skymodel_continuum this time')
    #    return

    # CHOICE2: overwrite files
    os.system("rm {0}".format(all_gals_fname+ ".fits"))
    os.system("rm {0}".format(all_gals_fname + "_z.fits"))
    os.system("rm {0}".format(all_gals_fname + "_maxflux.fits"))
    os.system("rm {0}".format(all_gals_fname + "_P.fits"))
    os.system("rm {0}".format(all_gals_fname + "_Q.fits"))#a
    os.system("rm {0}".format(all_gals_fname + "_U.fits"))#a

    
    logging.info("Creating empty image file, {0} ...".format(all_gals_fname))

    # exploting the expand_if_needed functionality of the image write function in fitsio

    logging.info(
        "Padded array dimensions: %d %d %d" % (n_chan, image_size, image_size)
    )
    logging.info(
        "Padded array size, elements: %.2e" % (n_chan * image_size * image_size)
    )
    logging.info(
        "Padded array size, bytes: %.2e" % (n_chan * image_size * image_size * 4)
    )

    test_array_size = 0
    if test_array_size:
        data = np.zeros((n_chan, image_size, image_size)).astype(np.float32)
        print("nbytes, 32 bits precision: ", data.nbytes)

    header_dict = {}
    logging.info("Making header:")

    # convert astropy header to dict for fitsio
    for key in header_fourd:
        header_dict[key] = header_fourd[key]
        logging.info("%s = %s" % (key, header_dict[key]))

    # start pixel coords need to be supplied in reverse order when expanding - fiddly but works
    # fitsio code changed to reverse axes ordering when writing (going from C to fortran order)

    blc0 = image_size - 1  # -1 since the large cube is added using Fortran indexing
    blc1 = image_size - 1
    blc2 = n_chan - 1

    # initialise two empy objects: one nfreq X npix X npix cube and one npix X npix map
    img2 = np.zeros((blc2 + 1, blc0 + 1, blc1 + 1)).astype(
        np.float32
    )  # empty array for the continuum cube
    img2_1D = np.zeros((1, blc0 + 1, blc1 + 1)).astype(
        np.float32
    )  # empy array for 1D maps



    cr_x, cr_y = w_twod.wcs_world2pix(
        ra_field_gs,
        dec_field_gs,
        1,
    )

    logging.info("Check world2pix crpix1, crpix2: %f %f", cr_x, cr_y)
    
    
  
    # here we initialise 3 files:
    # 1) continuum cube with all sources for the chosen continuum channels
    # 2) redshift map containing the redshift of the brightest source on the LoS - needed for HI absorption
    # 3) map of the brightest fluxes along the LoS - This is needed for obtaining 3 and could just be a trowaway dummy map. For the moment I am keeping for debugging purposes

    # initialise files 
    fitsf = FITS(all_gals_fname+".fits", "rw")
    fitsf_f = FITS(all_gals_fname + "_maxflux.fits", "rw")
    fitsf_z = FITS(all_gals_fname + "_z.fits", "rw")


    #modify the headers to reflect logaritmically-spaced frequency
    if (logfreq == True):
        header_dict['CTYPE3'] = 'FREQ-LOG'
        header_dict['CDELT3'] = np.log(dnu)


    # write those 3 empty files    
    fitsf.write(img2, header=header_dict)
    fitsf_f.write(img2_1D, header=header_dict)
    fitsf_z.write(img2_1D, header=header_dict)

    fitsf.close()
    fitsf_f.close()
    fitsf_z.close()


    # in the POLARIZATION case:
    # we also have 3 more files;
    # 4) polarised intensity
    # 5) Q Stokes
    # 6) U Stokes
    # the _maxflux and _z maps as in 2) and 3) refer to polarised flux instead of total intensity flux. This product is used for rotation measure propagation. 
    
    polarization = False
    if config.getboolean("skymodel", "dopolarization") == True:
        polarization = True

    if (polarization == True):
        fitsf_p = FITS(all_gals_fname + "_P.fits", "rw")
        fitsf_q = FITS(all_gals_fname + "_Q.fits", "rw")
        fitsf_u = FITS(all_gals_fname + "_U.fits", "rw")

        fitsf_p.write(img2, header=header_dict)
        #modify CRVAL4 for U
        header_dict['CRVAL4']=2
        
        fitsf_q.write(img2, header=header_dict)

        #modify CRVAL4 for U
        header_dict['CRVAL4']=3
        fitsf_u.write(img2, header=header_dict)
    

        fitsf_p.close()
        fitsf_q.close()
        fitsf_u.close()
        
    
    # READING T-RECS CATALOGUE AND IMPORTING RELEVANT QUANTITIES
    # this flag indicates whether the catalogue is a continuum (False) or
    # continuum X HI (True) catalogue
    HI_cross = False  #initialised as continuum only, to be revised when inspecting the catalogue

    # Load the catalogue and import column names
    cat_file_name = config.get("field", "catalogue")
    logging.info("Loading catalogue from {0} ...".format(cat_file_name))
    cat = Table()
    cat_read = Table.read(cat_file_name)  # remove ascii
    keywords = cat_read.colnames
    logging.info("Done")
    # use the HI mass keyword to identify an HIXcontinuum catalogue
    if "MHI" in keywords:
        HI_cross = True
    
    cat["id"] = np.zeros(len(cat_read))
    source_prefix = "TRECS-"
    cat["Source_id"] =cat_read["ID_cont"]  # source unique identifier
    cat["RA"] = cat_read["longitude"]  # deg
    cat["ra_offset"] = cat["RA"] - ra_field_gs  # deg
    cat["ra_offset"].unit = "deg"
    cat["DEC"] = cat_read["latitude"]  # deg
    cat["dec_offset"] = cat_read["latitude"] - dec_field_gs  # deg
    cat["dec_offset"].unit = "deg"
    z = cat_read["redshift"]
  
    # Each source frequency behaviour is approximated as a power law within the cube. A spectral index is computed between the lowest and highest specified frequencies.
    # This approximation is OK for channels within the same band.
    # For frequencies belonging do different bands, perform multiple runs of the code.

    cat["Total_flux"] = cat_read["I" + base_freqname] * 1.0e-3  # Jy
    cat["Total_flux"].unit = "Jy"
    cat["flux2"] = cat_read["I" + top_freqname] * 1.0e-3  # Jy
    cat["flux2"].unit = "Jy"
    cat["spectral_index"] = np.log10(
        cat["Total_flux"] / cat["flux2"]
    ) / np.log10(base_freq / top_freq)


    #Initialise polarization fraction as 1 so that P=I if polarization ==False
    cat["polafrac"] = cat["Total_flux"]*0.+1.
    
    if (polarization == True):
         cat["polafrac"] = cat_read["P" + base_freqname] * 1.0e-3 /cat["Total_flux"]
         # generate polarization angle
         # convention it is 0 at North and anticlockwise
         np.random.seed(mother_seed + 1093548)
         evpa = np.random.uniform(low=0., high=359.99, size=len(cat)
         )  
         cat["EVPA"] = evpa

    # source size and morphology description
    maj = cat_read["size"]  # arcsec
    cat["Maj"] = maj
    cat["Maj"].unit = "arcsec"
    q = cat_read["axis ratio"]
    
    cat["Min"] = maj * q
    cat["Min"].unit = "arcsec"
    cat["Rs"] = cat_read["Rs"]
    rdcl = cat_read["RadioClass"]  # to be used in source selection
    cat["RadioClass"] = rdcl

    #source position angle
    pa = cat_read["PA"]  # position angle

    if HI_cross == True:
        # PA in continuum is the HI PA rotated by 90 degs
        pa_copy = pa + 90.0
        pa_copy[pa_copy > 359.0] = pa_copy[pa_copy > 359.0] - 360.0
        pa_copy[pa == -100] = -100.0
        pa[rdcl > 3] = pa_copy[rdcl > 3]  # AGN PA 90degs from HI PA.
        pa_1 = cat_read["PA_1"]
        pa[pa == -100] = pa_1[pa == -100]
        pa_1 = 0
        pa_copy = 0
        
    # PA not defined for AGN, here generate random
    np.random.seed(mother_seed + 1)
    pa_2 = np.random.uniform(low=0., high=359.99, size=len(cat))
    pa[pa == -100] = pa_2[pa == -100]
    pa_2 = 0
    
    cat["PA"] = pa
    cat["PA"].unit = "deg"

    # A cross-catalogue has the continuum sources as additional columns.
    # couterparts have values on the primary columns which should be used. Continuum-only sources need top be loaded from the additional columns

    if HI_cross == True:
        z_1 = cat_read["redshift_1"]
        z[z == -100] = z_1[z == -100]
        q1 = cat_read["axis ratio_1"]
        q[q == -100] = q1[q == -100]

        
    cat["z"] = z


    # SELECT THE RELEVANT PORTION OF THE T-RECS CATALOGUE
    # read the relevant quantities to implement flux cuts
    if config.getboolean("continuum", "highfluxcut") == True:
        highflux = config.getfloat("continuum", "highfluxcut_value")
        flux_sel_freq = config.get("continuum", "fluxcut_frequency")
        cat["flux_selection"] = cat_read["I" + flux_sel_freq] * 1.0e-3  # Jy

    if config.getboolean("continuum", "lowfluxcut") == True:
        lowflux = config.getfloat("continuum", "lowfluxcut_value")
        flux_sel_freq = config.get("continuum", "fluxcut_frequency")
        cat["flux_selection"] = cat_read["I" + flux_sel_freq] * 1.0e-3  # Jy

    # read the relevant quantities to implement redshift range selection
    if config.getboolean("continuum", "zrange") == True:
        zmin = config.getfloat("continuum", "zmin")
        zmax = config.getfloat("continuum", "zmax")

 
    ###APPLY CATALOGUE SELECTION
    

  
    # select continuum sources
    cat = cat[(cat["RadioClass"] != -100)]
    
    # exclude too big; memory problem and not realistic
    cat = cat[(cat["Maj"] < 2000.0)]

    
    if config.getboolean("continuum", "highfluxcut") == True:
        # Exclude sources brigther than a value
        print("applying high flux cut")
        len_old = len(cat)
        cat = cat[(cat["flux_selection"] < highflux)]
        print("number of sources excluded")
        print(len_old - len(cat))

    if config.getboolean("continuum", "lowfluxcut") == True:
        # Exclude sources fainter than a value
        print("applying low flux cut")
        len_old = len(cat)
        cat = cat[(cat["flux_selection"] > lowflux)]
        print("number of sources excluded")
        print(len_old - len(cat))

    if config.getboolean("continuum", "zrange") == True:
        # Exclude sources outside a redshift range
        print("applying redshift range cut")
        len_old = len(cat)
        cat = cat[(cat["z"] >= zmin)]
        cat = cat[(cat["z"] < zmax)]
        print("number of sources excluded")
        print(len_old - len(cat))
        print(min(cat["z"]),max(cat["z"]))
    
        
    # define additional source attributes not contained in the TRECS cat
    cat["Atlas_source"] = np.zeros(len(cat)).astype(np.str)
    cat["Unresolved"] = np.zeros(len(cat)).astype(np.str)
    cat["New_flux"] = np.zeros(len(cat)).astype(np.str)
    np.random.seed(mother_seed + 100)

    # initialise core fraction. Steep-spectrum AGN don't use it as it is determined by the postage stamp.
    corefrac = np.random.normal(
        loc=0.75, scale=0.1, size=len(cat)
    )  
    cat["corefrac"] = corefrac
    np.random.seed(mother_seed + 1000)

    # this random number is used later to associate sources to postage stamps
    ranid = np.random.uniform(low=0, high=1, size=len(cat))
 
    ranid[cat["Rs"] <= 0.5] = ranid[cat["Rs"] <= 0.5] - 10
    ranid[cat["Rs"] > 0.5] = ranid[cat["Rs"] > 0.5] + 10
    cat["ranid"] = ranid

    #apply catalogue selection based on pixel projection as otherwise 
    #I get empty borders
    
    nobj=len(cat)
    selected=np.zeros(nobj)
    xs=np.zeros(nobj)
    ys=np.zeros(nobj)
    #print('ecco',nobj,cr_x,cr_y)
    logging.info('Selecting for pixel projections')
    for i, cat_gal in enumerate(cat):
        x, y = w_twod.wcs_world2pix(cat_gal["RA"],cat_gal["DEC"],1,)    
        xs[i] = float(x)
        ys[i] = float(y)
        # select ony spurces inside the FoV. keeping a 10 pixel padding for large sources 
        if (x >= -10) and (x <=2*cr_x+10) and (y >= -10) and (y <=2*cr_y+10):
            selected[i]=1
    cat=cat[selected==1]
    xs=xs[selected==1]
    ys=ys[selected==1]
    logging.info('Done')
    cat["xs"] = xs
    cat["ys"] = ys
            
    nobj = len(cat)
    logging.info('Catalogue total size %s',len(cat))

    
    # save the truth catalogue after selections on a pickle file
    # if run in parallel save just for the first process
    if (total_cores ==1):
        
        logging.info('Writing the truth catalogue')
        file = open(all_gals_fname+"_catalogue_pickle", 'wb')
        pickle.dump(cat, file)
        logging.info('Done')
            
    # if running in parallel: split the catalogue in similar-sized portions
    # the selection is done randomly so that all processes get sources in no particular redshift order

    if (total_cores >1):
        #if (process ==1):
            # if multiple processes, write the picke file only once
            
            #logging.info('Writing the truth catalogue')
            #file = open(all_gals_fname+"_catalogue.pickle", 'wb')
            #pickle.dump(cat, file)
            #logging.info('Done')
            

        print("number of objects to share between processes",nobj)
        np.random.seed(mother_seed + 3479)
        ranid2 = np.random.uniform(low=0., high=1., size=nobj)
        cat["ranid2"] = ranid2
        
        print('process ',process,'of ',total_cores)
        delta=1./total_cores
        cores_cut = (cat["ranid2"] >= delta*(process-1)) * (cat["ranid2"] < delta*process)


        #print('cores cut',cores_cut)
        print(delta*(process-1),delta*process)

        cat=cat[cores_cut]
        nobj=len(cat)
        print('Catalogue size for this process: ',nobj)
        logging.info('Catalogue size for this process: %f',nobj)

    # EXECUTE THE MAIN SIMULATION
    skipped_sources=0
    for i, cat_gal in enumerate(cat):
        add_source_continuum(i,cat_gal,nobj,w_twod,config,pixel_scale_str,psf_maj_arcsec,arr_dims,all_gals_fname,base_freq,freqs,polarization,skipped_sources)

    #todo: skipped_sources does not work, always returns zero.
    logging.info("Main loop done")
    logging.info("Skipped sources: %s",skipped_sources)
    

    print ('summed cube:', np.sum(astfits.getdata(all_gals_fname+".fits")))
    tend = time.time()
    logging.info("...done in {0} seconds.".format(tend - tstart))
    print("skymodel_continuum finished in {0} seconds.".format(tend - tstart))

# run the coadd and combine functions to get a single output set from
# outputs of multiple processes
def runCoadd(config,process,total_cores):
    """Simulate a sky model from a T-RECS catalogue.

    Parameters
    ----------
    config : configparser
        ConfigParser configuration containing necessary sections.

    """

    #edit confing object for the process number

               
    n_cores = int(config.getfloat("pipeline", "n_cores"))
   
    data_path_large_files = (
        config.get("pipeline", "data_path_large_files")
        + config.get("pipeline", "project_name")
        + "/"
    )
    data_path = (
        config.get("pipeline", "base_dir")
        + config.get("pipeline", "data_path")
        + config.get("pipeline", "project_name")
        + "/"
    )


    
    
    all_gals_fname = (
        data_path_large_files + config.get("field", "fits_prefix")
    )
    

    coadd(all_gals_fname,n_cores,"")

    if config.getboolean("skymodel", "dopolarization") == True:
        coadd(all_gals_fname,n_cores,"_P")
        coadd(all_gals_fname,n_cores,"_Q")
        coadd(all_gals_fname,n_cores,"_U")
    
    combine(all_gals_fname,n_cores)





if __name__ == "__main__":

    config = ConfigParser.ConfigParser()
    config.read(sys.argv[1])


