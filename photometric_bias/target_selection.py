import astropy.io.fits as fits
import numpy as np
import matplotlib.pyplot as plt
import healpy as hp
import math
from astropy.table import Table
import tarfile
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from astropy.io import fits
from astropy.table import Table
#import corner
import matplotlib.colors
import math
import astropy
import healpy as hp
import re
import argparse
import urllib.request, urllib.error, urllib.parse
import astropy.io.ascii as ascii
import sys,os

def target_selection(filename,tractlist):
    ra_all = np.array([])
    dec_all = np.array([])
    tract_all = np.array([])
    id_all = np.array([])

    g_all = np.array([])
    r_all = np.array([])
    i_all = np.array([])
    z_all = np.array([])
    y_all = np.array([])
    
    for tract_num in tractlist:
        num = str(tract_num)
        file = ("database/s21-colorterm/tracts_all/%s_bias.fits"%num)
        if os.path.exists(file):
            hdu = fits.open(file)
            data = hdu[1].data
            ra = data["i_ra"]
            dec = data["i_dec"]
            tract = data["tract"]

            id = data["object_id"]

            g_mag = data["forced_g_cmodel_mag"]
            r_mag = data["forced_r_cmodel_mag"]
            i_mag = data["forced_i_cmodel_mag"]
            z_mag = data["forced_z_cmodel_mag"]
            y_mag = data["forced_y_cmodel_mag"]

            g_err = data["forced_g_cmodel_magerr"]

            i_cmodel = data["meas_i_cmodel_mag"]
            i_psf = data["meas_i_psfflux_mag"]


####################################################extinction            
            g_a = data["a_g"]
            r_a = data["a_r"]
            i_a = data["a_i"]
            z_a = data["a_z"]
            y_a = data["a_y"]
            target_ra = (ra%360)/180*math.pi
            target_dec = math.pi/2 - dec/180*math.pi
            healpix = hp.ang2pix(256,target_dec, target_ra)
            
            hdu1 = fits.open('csfd_ebv.fits')
            data1 = hdu1[1].data
            im=data1["T"]
            m = np.ndarray.flatten(im)
            im = hp.ud_grade(m, 256)

            #g_a = 3.214*im[healpix]
            #r_a = 2.165*im[healpix]
            #i_a = 1.592*im[healpix]
            #z_a = 1.211*im[healpix]
            #y_a = 1.064*im[healpix]

            g = g_mag - g_a
            r = r_mag - r_a
            i = i_mag - i_a
            z = z_mag - z_a
            y = y_mag - y_a
        
            ####target selection###########
            #quality = (g_err<g_mag*0.05-1.1)&(g>22.5)&(g<24.5)
            quality = (g_err<g_mag*0.05-1.1)&(g_mag>22.5)
            #redshift = ((g-r)<0.15)|((i-y)>2.0*(g-r)-0.15)
            redshift = ((g_mag-r_mag)<0.15)|((i_mag-y_mag)>2.0*(g_mag-r_mag)-0.15)
            galaxy = (i_cmodel - i_psf < -0.15)
            selection = quality & redshift & galaxy
            ###################################
            _ra = ra[selection]
            _dec = dec[selection]
            _tract = tract[selection]
            _id = id[selection]
            
            _g = g_mag[selection]
            _r = r_mag[selection]
            _i = i_mag[selection]
            _z = z_mag[selection]
            _y = y_mag[selection]
            
            ra_all = np.concatenate([ra_all,_ra])
            dec_all = np.concatenate([dec_all,_dec])
            tract_all = np.concatenate([tract_all,_tract])
            id_all = np.concatenate([id_all, _id])
            
            g_all = np.concatenate([g_all, _g])
            r_all = np.concatenate([r_all, _r])
            i_all = np.concatenate([i_all, _i])
            z_all = np.concatenate([z_all, _z])
            y_all = np.concatenate([y_all, _y])

    t = Table([ra_all, dec_all, id_all, g_all, r_all, i_all, z_all, y_all, tract_all], names=('ra', 'dec', 'id', 'forced_g_cmodel_mag', 'forced_r_cmodel_mag', 'forced_i_cmodel_mag', 'forced_z_cmodel_mag', 'forced_y_cmodel_mag', 'tract'))
    t.write(filename, format='fits')