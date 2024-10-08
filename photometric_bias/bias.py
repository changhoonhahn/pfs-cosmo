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
import hscReleaseQuery as HQR
import target_GPU as selection
import target_selection as tract_selection
import sys,os
import mask_dot as Mask
from matplotlib.backends.backend_pdf import PdfPages
from multiprocessing import Pool

version =   20190924.1
args    =   None
doDownload  = False
doUnzip =   False
starUnzip = False
diffver =   '-colorterm'
nside = 256
area = hp.nside2pixarea(nside,degrees=True)
SQL = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--file", default ="output_test.fits", help = "photometry-data filename")
    #parser.add_argument("--photometry", nargs = "*", choices = ["gseeing", "rseeing", "iseeing", "zseeing", "yseeing", "gdepth", "rdepth", "idepth", "zdepth", "ydepth", "star", "star1", "star2", "star3", "extinction"], help = "properties")
    parser.add_argument('--user', '-u', required=True, help='specify your STARS account')
    parser.add_argument('--password-env', default='HSC_SSP_CAS_PASSWORD', help='environment variable for STARS')
    parser.add_argument('sql-file', type=argparse.FileType('r'), help='SQL file')
    
    parser.add_argument('--api-url',
                        default='https://hscdata.mtk.nao.ac.jp/datasearch/api/catalog_jobs/',
                        help='for developers')
    parser.add_argument('--format', '-f', dest='out_format', default='fits',
                        choices=['csv', 'csv.gz', 'sqlite3', 'fits'],
                        help='specify output format')
    parser.add_argument('--release_year', '-r', default = 's21',
                        choices='s15 s16 s17 s18 s19 s20 s21'.split(),
                        help='the release year')
    parser.add_argument('--nomail', '-M', action='store_true',
                        help='suppress email notice')

    parser.add_argument('--preview', '-p', action='store_true',
                        help='quick mode (short timeout)')
    parser.add_argument('--skip-syntax-check', '-S', action='store_true',
                        help='skip syntax check')
    parser.add_argument('--delete-job', '-D', action='store_true',
                        help='delete job after your downloading')


    #If only using parts of the whole HSC data
    parser.add_argument("--tracts", nargs = "*", help = "properties", default = '')

    pdf = PdfPages("bias_i_meas_new.pdf")

    global args,release_version,prefix,prefix2,area, tractlist,sql

    args = parser.parse_args()
    sql         =   args.__dict__['sql-file'].read()
    release_year   =   "s21"
    release_version =   'dr4-citus'
    prefix  =   'database/%s%s/sql_all' %(release_year,diffver)
    prefix2 =   'database/%s%s/tracts_all' %(release_year,diffver)
    if not os.path.exists(prefix):
        os.system('mkdir -p %s' %prefix)
    if not os.path.exists(prefix2):
        os.system('mkdir -p %s' %prefix2)

    filename = args.file
    #photometry = args.photometry
    photometry = ["gseeing", "rseeing", "iseeing", "zseeing", "yseeing", "gdepth", "rdepth", "idepth", "zdepth", "ydepth", "star", "star1", "star2", "star3", "extinction"]
    tractlist = args.tracts
    ngroups = 1

    property = {}
    global rangex
    rangex = {}
    for photo in photometry:
        #dict[photo] = data[photo]
        property[photo] = []
        
    rangex["depth"] = np.linspace(24,28,51) #define xaxis for histogram
    rangex["size"] = np.linspace(0.5,2.5,51)
    rangex["extinction"] = np.linspace(0.0,5.0,51)
    rangex["star"] = np.linspace(0.0,15000,21)

        ###################################################################    
    #if no tracts, all tract in HSC database will be downloaded
    if (tractlist==''):
        tractname= 'fieldTractInfoS21.csv'
        tracts      =   ascii.read(tractname)['tract']
        tractlist = list(tracts)
        ngroups = 20

    tracts2     =   HQR.chunkNList(list(tractlist),ngroups)

    ##############################################################
    #download photometry-datafile from HSC database
    if doDownload:
        global credential
        credential  =   {'account_name': args.user, 'password': HQR.getPassword(args)}
        for ig,tractL in enumerate(tracts2):
            #restart
            #if(ig<194): continue
            print('Group: %s' %ig)
            HQR.downloadTracts(args, prefix, release_version, credential, sql,ig,tractL)

    if doUnzip:
        for ig,tractL in enumerate(tracts2):
            print('unzipping group: %s' %ig)
            HQR.separateTracts(args, prefix,prefix2,ig,tractL)

    #################### get mask and all the healpixels that will be included in the analysis

    eff_area = np.array([])
    density = np.array([])
    healpix = np.array([], int)
    healpix_flags = np.array([], int)
    mask_flags = np.array([])
    flag = np.array([])

    if not os.path.exists("mask"):
        os.system("mkdir -p mask")
    for i in tractlist:
        outfname    =   os.path.join('mask/tract_%s.fits' %(i))
        if not (os.path.exists('mask/tract_%s.fits' %(i))| os.path.exists('mask/tract_%s_left.fits' %(i))):
            ra,dec = Mask.area(i)
            if (max(ra)-min(ra)>300):
                file = "tract_%s_left.fits"%(i)
                command = "HSC-SSP_brightStarMask_Arcturus/venice-4.0.3/bin/venice -r -xmin %s -xmax %s -ymin %s -ymax %s -coord spher -o mask/tract_%s_left.fits"%(0, min(ra), min(dec), max(dec),i)
                os.system(command)

                command1 = "HSC-SSP_brightStarMask_Arcturus/venice-4.0.3/bin/venice -m HSC-SSP_brightStarMask_Arcturus/reg/masks_all.reg -cat mask/tract_%s_left.fits -xcol ra -ycol dec -f all -flagName isOutsideMask -o mask/tract_%s_flagged_left.fits"%(i,i)
                os.system(command1)
                
                command = "HSC-SSP_brightStarMask_Arcturus/venice-4.0.3/bin/venice -r -xmin %s -xmax %s -ymin %s -ymax %s -coord spher -o mask/tract_%s_right.fits"%(max(ra), 360, min(dec), max(dec),i)
                os.system(command)                
                command1 = "HSC-SSP_brightStarMask_Arcturus/venice-4.0.3/bin/venice -m HSC-SSP_brightStarMask_Arcturus/reg/masks_all.reg -cat mask/tract_%s_right.fits -xcol ra -ycol dec -f all -flagName isOutsideMask -o mask/tract_%s_flagged_right.fits"%(i,i)
                os.system(command1)

                hdu = fits.open("mask/tract_%s_flagged_left.fits" %i)
                data = hdu[1].data
                ra_left = data["ra"]
                dec_left = data["dec"]
                flag_left = data["IsOutsideMask"]
                hdu.close()
                
                hdu = fits.open("mask/tract_%s_flagged_right.fits" %i)
                data = hdu[1].data
                ra_right = data["ra"]
                dec_right = data["dec"]
                flag_right = data["IsOutsideMask"]
                hdu.close()
                
                
                ra = np.concatenate([ra_left,ra_right])
                dec = np.concatenate([dec_left,dec_right])
                flag = np.concatenate([flag_left, flag_right])

                t = Table([ra, dec, flag], names=('ra', 'dec', "IsOutsideMask"))
                t.write('mask/tract_%s_flagged.fits' %(i), format='fits')
            
            
            else:
                command = "HSC-SSP_brightStarMask_Arcturus/venice-4.0.3/bin/venice -r -xmin %s -xmax %s -ymin %s -ymax %s -coord spher -o mask/tract_%s.fits"%(min(ra), max(ra), min(dec), max(dec),i)
                os.system(command)
                command1 = "HSC-SSP_brightStarMask_Arcturus/venice-4.0.3/bin/venice -m HSC-SSP_brightStarMask_Arcturus/reg/masks_all.reg -cat mask/tract_%s.fits -xcol ra -ycol dec -f all -flagName isOutsideMask -o mask/tract_%s_flagged.fits"%(i,i)
                os.system(command1)
        else:
            print("already have mask")


    for i in tractlist:
        hdu = fits.open("mask/tract_%s_flagged.fits" %(i))
        data = hdu[1].data
        mask_ra = data["ra"]/180*math.pi
        mask_dec = math.pi/2 - data["dec"]/180*math.pi
        mask_flag = data["IsOutsideMask"]
        healpix_mask = hp.ang2pix(nside,mask_dec,mask_ra)
        hdu.close()

        unique_values, counts = np.unique(healpix_mask, return_counts=True)
        healpix_above_threshold = unique_values[counts > 16000]
        healpix = np.concatenate((healpix, healpix_above_threshold))    ####exclude the healpixels that are on the edge of the tract

        mask_condition = np.isin(healpix_mask, healpix_above_threshold)   ####exclude flags that are on the edge healpix
        mask_flags = np.concatenate((mask_flags, mask_flag[mask_condition]))
        healpix_flags = np.concatenate((healpix_flags, healpix_mask[mask_condition])) #healpix of the selected flags
	
    print("Get the healpixel for the observation area")    
    healpix = np.unique(healpix)    ############################### The entire healpix

    
    npix = hp.nside2npix(nside)
    # select flags outside stellar mask
    pixels_flagged = healpix_flags[mask_flags == 1.0]
    pixel_count = np.bincount(pixels_flagged, minlength=npix) # number of flags with flag=1 (for the entire healpixel)

    all_count = np.bincount(healpix_flags, minlength=npix) # number of flags for the entire healpixel

    eff_area = pixel_count[all_count>0]/all_count[all_count>0]*area  ############################## effective area of the entire healpixel

    #flag the healpixels that are on the edge of the observation area
    ra1 = np.array([])
    dec1 = np.array([])
    mask_edge = np.array([True]*len(healpix))
    for healpy in healpix:
        #location = hp.pix2ang(256, int(healpy), lonlat=True)
        location = hp.pix2ang(256, int(healpy))
        ra1 = np.append(ra1, location[1])
        dec1 = np.append(dec1, location[0])

    #ra1 = ra1 - (ra1>300)*360
    #autumn = ((ra1>30)&(ra1<38) &(dec1<2.6)&(dec1 > -6.0))|((ra1<27.5)&(ra1>-27)&(dec1>-1.1)&(dec1<2.62))
    #spring = (ra1<217.7)&(ra1>130)&(dec1<3.0)&(dec1>-1.67)
    #north = (dec1 > 42.7)&(dec1<44) &(ra1 >205) &(ra1<245)
    ra1 = ra1 - 2*math.pi*(ra1>5)
    autumn = ((ra1>0.52)&(ra1<0.68) &(dec1<1.69)&(dec1 > 1.525))|((ra1<0.48)&(ra1>-0.5)&(dec1>1.525)&(dec1<1.59))
    spring = (ra1<3.8)&(ra1>2.25)&(dec1<1.6)&(dec1>1.5)
    north = (dec1 > 0.8)&(dec1<0.825) &(ra1 >3.5) &(ra1<4.35)
    mask_edge = autumn | spring | north
    print(ra1)
    print(dec1)

    #######################################
    #get the location of the target galaxies
    if os.path.exists(filename):
        hdu = fits.open(filename)
        data = hdu[1].data
    else:
        if (args.tracts==''):
            command = f"python target_GPU.py --file {filename}"
            os.system(command)
            #selection.target_selection(filename,list(tractlist))
            hdu = fits.open(filename)
            data = hdu[1].data
        else:
            tract_selection.target_selection(filename,list(tractlist))
            hdu = fits.open(filename)
            data = hdu[1].data            

    target_tract = data["tract"]
    target_ra = (data["ra"]%360)/180*math.pi
    target_dec = math.pi/2 - data["dec"]/180*math.pi
    target_healpix = hp.ang2pix(nside,target_dec,target_ra)
    print("Get the location of the target galaxies")
    hdu.close()

    number_target = np.bincount(target_healpix, minlength=npix)
    density = number_target[healpix]/eff_area                      ################### density of healpixels



    ################################################################
    #get the location and the photometry of the entire objects
    if (args.tracts==''):
        ra_all, dec_all, i_cmodel, i_psf, g_seeing, r_seeing, i_seeing, z_seeing, y_seeing, g_depth, r_depth, i_depth, z_depth, y_depth = Get_Photometry(SQL)
    else:
        ra_all, dec_all, i_cmodel, i_psf, g_seeing, r_seeing, i_seeing, z_seeing, y_seeing, g_depth, r_depth, i_depth, z_depth, y_depth = Get_Photometry_tract(tractlist)

    ra_all = (ra_all%360)/180*math.pi
    dec_all = math.pi/2 - dec_all/180*math.pi
    healpix_all = hp.ang2pix(nside, dec_all, ra_all) ######### get the healpixel for entire object
    all_count = np.bincount(healpix_all, minlength=npix)    ###### number of object in each healpixel

    property["gseeing"] = healpix_photometry(healpix_all, g_seeing, healpix, all_count)
    property["rseeing"] = healpix_photometry(healpix_all, r_seeing, healpix, all_count)
    property["iseeing"] = healpix_photometry(healpix_all, i_seeing, healpix, all_count)
    property["zseeing"] = healpix_photometry(healpix_all, z_seeing, healpix, all_count)
    property["yseeing"] = healpix_photometry(healpix_all, y_seeing, healpix, all_count)

    property["gdepth"] = healpix_photometry(healpix_all, g_depth, healpix, all_count)
    property["rdepth"] = healpix_photometry(healpix_all, r_depth, healpix, all_count)
    property["idepth"] = healpix_photometry(healpix_all, i_depth, healpix, all_count)
    property["zdepth"] = healpix_photometry(healpix_all, z_depth, healpix, all_count)
    property["ydepth"] = healpix_photometry(healpix_all, y_depth, healpix, all_count)

    ########get star
    star_pixel = healpix_all[(i_cmodel - i_psf > -0.08)&(i_cmodel<22.0)]
    star_count = np.bincount(star_pixel, minlength=npix)
    property["star"] = star_count[healpix]/eff_area

    star1_pixel = healpix_all[(i_cmodel - i_psf > -0.05)&(i_cmodel<22.0)]
    star1_count = np.bincount(star1_pixel, minlength=npix)
    property["star1"] = star1_count[healpix]/eff_area

    star2_pixel = healpix_all[(i_cmodel - i_psf > -0.03)&(i_cmodel<22.0)]
    star2_count = np.bincount(star2_pixel, minlength=npix)
    property["star2"] = star2_count[healpix]/eff_area

    star3_pixel = healpix_all[(i_cmodel - i_psf > -0.01)&(i_cmodel<22.0)]
    star3_count = np.bincount(star3_pixel, minlength=npix)
    property["star3"] = star3_count[healpix]/eff_area
    print("Get stellar density")
    ################extinction
    hdu = fits.open('csfd_ebv.fits')
    data = hdu[1].data
    im=data["T"]
    m = np.ndarray.flatten(im)
    im = hp.ud_grade(m, nside)
    property["extinction"] = im[healpix]
    print("Get dust extinction")
    ##########################################################

    val = list(property.values())
    val.append(list(density))
    val.append(list(healpix))
    val.append(list(mask_edge))

    nam = list(property.keys())
    nam_all = tuple(nam + ["density" , "healpix", "edge_flag"])
    t = Table(val, names=nam_all)
    t.write("property_new.fits", format='fits', overwrite = True)

    fig = plt.figure(figsize = (18,6))

    i = 1
    print("Get plots")
    ax1 = fig.add_subplot(2,6,1)
    ax1.set_xlabel("stellar density [deg-2]")
    ax1.set_title("stellar density")
    for k,v in property.items():
        if ("star" in k):
            histogram_star(fig, density[mask_edge], np.array(v)[mask_edge], ax1)
        else:
            i+=1
            histogram(fig, density[mask_edge],k,np.array(v)[mask_edge],i)


    fig.subplots_adjust(hspace=0.5, wspace=0.5)
    plt.tight_layout()
    pdf.savefig(fig)
    pdf.close()

        
        
def histogram_star(fig, dens,value, ax):
    x = rangex["star"]
    number = []
    mean_density =[]
    SD_density = []
    x1 = []

    #propertyごとにhealpixをわけてmean density, SD, number of healpix　を計算
    for i in range(len(x)-1):
        healpix = (x[i]<=value)&(value<x[i+1])
        number.append(np.sum(healpix))
        mean_density.append(np.average(dens[healpix]))
        SD_density.append(np.std(dens[healpix]))

        x1.append((x[i]+x[i+1])/2)

    x1 = np.array(x1)
    number = np.array(number)
    mean_density = np.array(mean_density)[number>100]
    SD_density = np.array(SD_density)[number>100]
    density_ave = np.average(mean_density)
    print(mean_density/density_ave)
    
    ax.errorbar(x1[number>100],mean_density/density_ave - 1.0,yerr = SD_density/density_ave ,fmt = "o")
    ax.plot(x1[number>100], mean_density/density_ave - 1.0)
    ax1 = ax.twinx()
    ax1.hist(value, bins=20, range =(np.min(x1[number>100]),np.max(x1[number>100])), weights=[1/len(value)]*len(value))
    ax.set_ylim(-0.5,0.5)
    ax1.set_ylim(0.0,1.0)
    ax1.tick_params(labelright=False)
    return
    

def histogram(fig, dens,key,value,i):
    print(dens)
    if "seeing" in key:
        ax = fig.add_subplot(2,6,i)
        x = rangex["size"]
        if key == "gseeing":
            ax.set_xlabel("gband psf size")
            ax.set_title("gband psf size [arcsec]")
        elif key == "rseeing":
            ax.set_xlabel("rband psf size")
            ax.set_title("rband psf size [arcsec]")
        elif key == "iseeing":
            ax.set_xlabel("iband psf size")
            ax.set_title("iband psf size [arcsec]")
        elif key == "zseeing":
            ax.set_xlabel("zband psf size")
            ax.set_title("zband psf size [arcsec]")
        elif key == "yseeing":
            ax.set_xlabel("yband psf size")
            ax.set_title("yband psf size [arcsec]")
    elif "depth" in key:
        ax = fig.add_subplot(2,6,i)
        x = rangex["depth"]
        if key == "gdepth":
            ax.set_xlabel("gband depth")
            ax.set_title("gband depth [mag]")
        elif key == "rdepth":
            ax.set_xlabel("rband depth")
            ax.set_title("rband depth [mag]")
        elif key == "idepth":
            ax.set_xlabel("iband depth")
            ax.set_title("iband depth [mag]")
        elif key == "zdepth":
            ax.set_xlabel("zband depth")
            ax.set_title("zband depth [mag]")
        elif key == "ydepth":
            ax.set_xlabel("yband depth")
            ax.set_title("yband depth [mag]")
    elif "extinction" in key:
        ax = fig.add_subplot(2,6,i)
        x = rangex["extinction"]
        ax.set_xlabel("E(B-V)")
        ax.set_title("dust extinction")

    number = []
    mean_density =[]
    SD_density = []
    x1 = []

    #propertyごとにhealpixをわけてmean density, SD, number of healpix　を計算
    for i in range(len(x)-1):
        healpix = (x[i]<=value)&(value<x[i+1])
        number.append(np.sum(healpix))
        mean_density.append(np.average(dens[healpix]))
        SD_density.append(np.std(dens[healpix]))

        x1.append((x[i]+x[i+1])/2)

    x1 = np.array(x1)
    number = np.array(number)
    mean_density = np.array(mean_density)[number>100]
    SD_density = np.array(SD_density)[number>100]
    density_ave = np.average(mean_density)
    print(mean_density/density_ave)
    
    ax.errorbar(x1[number>100],mean_density/density_ave - 1.0,yerr = SD_density/density_ave ,fmt = "o")
    ax.plot(x1[number>100], mean_density/density_ave - 1.0)
    ax1 = ax.twinx()
    ax1.hist(value, bins=20, range =(np.min(x1[number>100]),np.max(x1[number>100])), weights=[1/len(value)]*len(value))
    ax.set_ylim(-0.5,0.5)
    ax1.set_ylim(0.0,1.0)
    ax1.tick_params(labelright=False)

    return

def download(tract):
    tractstr = str(tract)
    
    job         =   None
    sqlU        =   sql.replace('{$tract}',tractstr)
    
    #fits file name will be tract number.fits
    outfname    =   '%s_bias.fits'%(tractstr)
    outfname    =   os.path.join(prefix,outfname)
    if os.path.exists(outfname):
        print('already have output')
        return
    print('querying data')
    job         =   submitJob(credential, sqlU, args.out_format)
    blockUntilJobFinishes(credential, job['id'])
    print('downloading data')
    fileOut     =   open(outfname,'w')
    fileBuffer  =   fileOut.buffer
    download(credential, job['id'], fileBuffer)
    #if args.delete_job:
        #deleteJob(credential, job['id'])
    print('closing output file')
    fileBuffer.close()
    fileOut.close()
    del fileBuffer
    del fileOut
    return

def healpix_photometry(healpix_all, photometry, healpix, count):
    print("calculating photometry")
    # 各healpix_targetのカウント（点の数）を取得
    npix = hp.nside2npix(nside)
    _property = np.bincount(healpix_all, weights=photometry, minlength=npix)
    Property = _property[healpix]/count[healpix]
    return Property

def process_sql(sql_num):
    file = f"database/s21-colorterm/sql_all/{sql_num}_bias.fits"
    print(file)
    if os.path.exists(file):
        with fits.open(file) as hdu:
            data = hdu[1].data
            ra = data["i_ra"]
            dec = data["i_dec"]
            i_cmodel = data["meas_i_cmodel_mag"]
            i_psf = data["meas_i_psfflux_mag"]

            g_seeing = data["gseeing"]
            r_seeing = data["rseeing"]
            i_seeing = data["iseeing"]
            z_seeing = data["zseeing"]
            y_seeing = data["yseeing"]

            g_depth = data["g_depth"]
            r_depth = data["r_depth"]
            i_depth = data["i_depth"]
            z_depth = data["z_depth"]
            y_depth = data["y_depth"]

            return ra, dec, i_cmodel, i_psf, g_seeing, r_seeing, i_seeing, z_seeing, y_seeing, g_depth, r_depth, i_depth, z_depth, y_depth
    return None

def process_tract(tract):
    file = f"database/s21-colorterm/tracts_all/{tract}_bias.fits"
    print(file)
    if os.path.exists(file):
        with fits.open(file) as hdu:
            data = hdu[1].data
            ra = data["i_ra"]
            dec = data["i_dec"]
            i_cmodel = data["meas_i_cmodel_mag"]
            i_psf = data["meas_i_psfflux_mag"]

            g_seeing = data["gseeing"]
            r_seeing = data["rseeing"]
            i_seeing = data["iseeing"]
            z_seeing = data["zseeing"]
            y_seeing = data["yseeing"]

            g_depth = data["g_depth"]
            r_depth = data["r_depth"]
            i_depth = data["i_depth"]
            z_depth = data["z_depth"]
            y_depth = data["y_depth"]

            return ra, dec, i_cmodel, i_psf, g_seeing, r_seeing, i_seeing, z_seeing, y_seeing, g_depth, r_depth, i_depth, z_depth, y_depth
    return None

def Get_Photometry(sql):
    with Pool(processes=20) as pool:  # Adjust number of processes based on your CPU
        results = pool.map(process_sql, sql)
        
    #results = []
    #with ProcessPoolExecutor() as executor:
        #futures = [executor.submit(process_file, sql_num) for sql_num in sql]
        #for future in futures:
            #result = future.result()
            #if result:
                #results.append(result)

    # Flatten results
    ra, dec, i_cmodel, i_psf, g_seeing, r_seeing, i_seeing, z_seeing, y_seeing, g_depth, r_depth, i_depth, z_depth, y_depth = map(np.concatenate, zip(*results))

    # Write results to a file
    return ra, dec, i_cmodel, i_psf, g_seeing, r_seeing, i_seeing, z_seeing, y_seeing, g_depth, r_depth, i_depth, z_depth, y_depth

def Get_Photometry_tract(tract_list):
    with Pool(processes=20) as pool:  # Adjust number of processes based on your CPU
        results = pool.map(process_tract, tract_list)
        
    #results = []
    #with ProcessPoolExecutor() as executor:
        #futures = [executor.submit(process_file, sql_num) for sql_num in sql]
        #for future in futures:
            #result = future.result()
            #if result:
                #results.append(result)

    # Flatten results
    ra, dec, i_cmodel, i_psf, g_seeing, r_seeing, i_seeing, z_seeing, y_seeing, g_depth, r_depth, i_depth, z_depth, y_depth = map(np.concatenate, zip(*results))

    # Write results to a file
    return ra, dec, i_cmodel, i_psf, g_seeing, r_seeing, i_seeing, z_seeing, y_seeing, g_depth, r_depth, i_depth, z_depth, y_depth



if __name__ == '__main__':
    main()
