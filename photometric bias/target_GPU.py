import os
import math
import numpy as np
from astropy.io import fits
from astropy.table import Table
import healpy as hp
from multiprocessing import Pool

filename = "output_test.fits"
sql = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]
def process_file(sql_num):
    file = f"../../../mnt/data_cat5/yuka/database/s21_colorterm/sql_all/{num}_bias.fits"
    print(file)
    if os.path.exists(file):
        with fits.open(file) as hdu:
            data = hdu[1].data
            ra = data["i_ra"]
            dec = data["i_dec"]
            tract = data["tract"]
            id = data["object_id"]
            g_mag = data["meas_g_cmodel_mag"]
            r_mag = data["meas_r_cmodel_mag"]
            i_mag = data["meas_i_cmodel_mag"]
            z_mag = data["meas_z_cmodel_mag"]
            y_mag = data["meas_y_cmodel_mag"]
            g_err = data["meas_g_cmodel_magerr"]
            i_cmodel = data["meas_i_cmodel_mag"]
            i_psf = data["meas_i_psfflux_mag"]

            g_a = data["a_g"]
            r_a = data["a_r"]
            i_a = data["a_i"]
            z_a = data["a_z"]
            y_a = data["a_y"]

            g = g_mag - g_a
            print(g)
            r = r_mag - r_a
            i = i_mag - i_a
            z = z_mag - z_a
            y = y_mag - y_a

            quality = (i_mag > 22.5) & (i_mag < 24.0)
            redshift = ((g_mag - r_mag) < 0.15) | ((i_mag - z_mag) > 2.0 * (g_mag - r_mag) - 0.15)
            galaxy = (i_cmodel - i_psf < -0.08)
            location = (210.0 < ra) & (ra < 215.0) & (dec < 1.0) & (dec > -1.0)
            selection = quality & redshift & galaxy & location

            return ra[selection], dec[selection], tract[selection], id[selection], \
                   g_mag[selection], r_mag[selection], i_mag[selection], z_mag[selection], y_mag[selection]
    return None

def main():
    with Pool(processes=20) as pool:  # Adjust number of processes based on your CPU
        results = pool.map(process_file, sql)
        
    #results = []
    #with ProcessPoolExecutor() as executor:
        #futures = [executor.submit(process_file, sql_num) for sql_num in sql]
        #for future in futures:
            #result = future.result()
            #if result:
                #results.append(result)

    # Flatten results
    ra_all, dec_all, tract_all, id_all, g_all, r_all, i_all, z_all, y_all = map(np.concatenate, zip(*results))

    # Write results to a file
    t = Table([ra_all, dec_all, id_all, g_all, r_all, i_all, z_all, y_all, tract_all],
              names=('ra', 'dec', 'id', 'meas_g_cmodel_mag', 'meas_r_cmodel_mag', 'meas_i_cmodel_mag',
                     'meas_z_cmodel_mag', 'meas_y_cmodel_mag', 'tract'))
    t.write(filename, format='fits', overwrite=True)
    return

if __name__ == '__main__':
    main()

# Example usage
# target_selection("output.fits", tract_list)
