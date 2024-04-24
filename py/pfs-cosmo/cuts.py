'''
'''



def apply_cuts(objects, survey='cmx'): 
    ''' apply target selections on objects 


    parameters
    ----------



    returns
    -------

    '''
    
    _prepare_hsc_imaging(objects) 




def _prepare_hsc_imaging(objects): 
    '''
    '''

    du = fits.open(file)
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



