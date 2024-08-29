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
import target_selection as selection
import sys,os
import mask_dot as Mask
from matplotlib.backends.backend_pdf import PdfPages

version =   20190924.1
args    =   None
doDownload  =   True
doUnzip =   True
starUnzip = False
diffver =   '-colorterm'
nside = 256
area = hp.nside2pixarea(nside,degrees=True)

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--file", nargs ="*", help = "photometry-data filename")
    parser.add_argument("--photometry", nargs = "*", choices = ["gseeing", "rseeing", "iseeing", "zseeing", "yseeing", "gdepth", "rdepth", "idepth", "zdepth", "ydepth", "star", "extinction"], help = "properties")
    parser.add_argument('--user', '-u', required=True, help='specify your STARS account')
    parser.add_argument('--password-env', default='HSC_SSP_CAS_PASSWORD', help='environment variable for STARS')
    parser.add_argument('sql-file', type=argparse.FileType('r'), help='SQL file')
    parser.add_argument('--star_sql', type=argparse.FileType('r'), help='SQL file for stars')
    
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

    pdf = PdfPages("PFS_photometry_bias_20.pdf")

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

    filename = args.file[0]
    photometry = args.photometry
    tractlist = args.tracts
    ngroups = 1

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

##################################################################    
    #get the location of tha target galaxies
    if os.path.exists(filename):
        hdu = fits.open(filename)
        data = hdu[1].data
    else:
        selection.target_selection(filename,list(tractlist))
        hdu = fits.open(filename)
        data = hdu[1].data

    target_tract = data["tract"]
    target_ra = (data["ra"]%360)/180*math.pi
    target_dec = math.pi/2 - data["dec"]/180*math.pi

    hdu.close()
    
    ################################################################
    #################################################################
    #download mask from Arcturus
    #for i in tractlist:
        #outfname    =   os.path.join("mask",'tract_%s.fits' %(i))
        #if os.path.exists(outfname):
            #print('already have file for tract: %s' \
                    #%i)
            #continue
        #ra,dec = Mask.area(str(i))
        #command = "venice-4.0.3/bin/venice -r -xmin %s -xmax %s -ymin %s -ymax %s -coord spher -o tract_%s.fits"%(min(ra), max(ra), min(dec), max(dec),i)
        #os.system(command)

    #for i in tractlist:
        #outfname    =   os.path.join("mask",'tract_%s_flagged.fits' %(i))
        #if os.path.exists(outfname):
            #print('already have file for tract: %s' \
                    #%i)
            #continue
        #ra,dec = Mask.area(str(i))
        #command1 = "venice-4.0.3/bin/venice -m reg/tracts/BrightStarMask-%s-HSC-I.reg -cat tract_%s.fits -xcol ra -ycol dec -f all -flagName isOutsideMask -o tract_%s_flagged.fits"%(i,i,i)
        #os.system(command1)

    #os.system("mkdir -p mask")
    #os.system("mv tract*.fits mask")

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

                command1 = "HSC-SSP_brightStarMask_Arcturus/venice-4.0.3/bin/venice -m reg/masks_all.reg -cat mask/tract_%s_left.fits -xcol ra -ycol dec -f all -flagName isOutsideMask -o mask/tract_%s_flagged_left.fits"%(i,i)
                os.system(command1)
                
                command = "HSC-SSP_brightStarMask_Arcturus/venice-4.0.3/bin/venice -r -xmin %s -xmax %s -ymin %s -ymax %s -coord spher -o mask/tract_%s_right.fits"%(max(ra), 360, min(dec), max(dec),i)
                os.system(command)
                
                command1 = "HSC-SSP_brightStarMask_Arcturus/venice-4.0.3/bin/venice -m reg/masks_all.reg -cat mask/tract_%s_right.fits -xcol ra -ycol dec -f all -flagName isOutsideMask -o mask/tract_%s_flagged_right.fits"%(i,i)
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
                command1 = "HSC-SSP_brightStarMask_Arcturus/venice-4.0.3/bin/venice -m reg/masks_all.reg -cat mask/tract_%s.fits -xcol ra -ycol dec -f all -flagName isOutsideMask -o mask/tract_%s_flagged.fits"%(i,i)
                os.system(command1)
        else:
            print("already have mask")

    #os.system("mkdir -p mask")
    #os.system("mv tract*.fits mask")
    ################################################################ 
    
    property = {}
    global rangex
    rangex = {}
    for photo in photometry:
        #dict[photo] = data[photo]
        property[photo] = []
        
    rangex["depth"] = np.linspace(24,28,101) #define xaxis for histogram
    rangex["size"] = np.linspace(0.5,2.5,101)
    rangex["extinction"] = np.linspace(0.0,5.0,101)
    rangex["star"] = np.linspace(0.0,150000,31)

    eff_area = []
    density = []
    healpix = np.array([])
    
    for tract in tractlist:
        if os.path.exists("database/s21-colorterm/tracts_all/%s_bias.fits" %tract):
        #if star density in property, get the location of stars in tract##############
            if ("star" in photometry):
                hdu = fits.open("database/s21-colorterm/tracts_all/%s_bias.fits" %tract)
                data = hdu[1].data

                cmodel = data["meas_i_cmodel_mag"]
                psf = data["meas_i_psfflux_mag"]

                ext =(cmodel - psf > -0.08)
                #location
                star_ra = (data["i_ra"][ext]%360)/180*math.pi
                star_dec = math.pi/2 - data["i_dec"][ext]/180*math.pi
                global star_pixel
                star_pixel = hp.ang2pix(nside, star_dec, star_ra)
                hdu.close()
            ###############################################################################
            #get galaxy property of all galaxy in tract
            hdu = fits.open("database/s21-colorterm/tracts_all/%s_bias.fits" %tract)
            data1 = hdu[1].data

            global all_pixel
            #location
            all_ra = (data1["i_ra"]%360)/180*math.pi
            all_dec = math.pi/2 - data1["i_dec"]/180*math.pi
            all_pixel = hp.ang2pix(nside, all_dec, all_ra)
        
            #photometry
            g_mag = data1["forced_g_cmodel_mag"]
            r_mag = data1["forced_r_cmodel_mag"]
            i_mag = data1["forced_i_cmodel_mag"]
            z_mag = data1["forced_z_cmodel_mag"]
            y_mag = data1["forced_y_cmodel_mag"]

            g_err = data1["forced_g_cmodel_magerr"]
            r_err = data1["forced_r_cmodel_magerr"]
            i_err = data1["forced_i_cmodel_magerr"]
            z_err = data1["forced_z_cmodel_magerr"]
            y_err = data1["forced_y_cmodel_magerr"]

            hdu.close()
            ####################################################################################
            #target galaxyと all galaxyをhealpixelに分ける
            healpix_all = hp.ang2pix(nside,all_dec,all_ra)
            healpix_target = hp.ang2pix(nside,target_dec[target_tract == float(tract)],target_ra[target_tract == float(tract)])
            healpix_range = np.unique(healpix_all)

            #mask dotsをhealpixelに分ける###########
            hdu = fits.open("mask/tract_%s_flagged.fits" %(tract))
            data = hdu[1].data
            mask_ra = data["ra"]/180*math.pi
            mask_dec = math.pi/2 - data["dec"]/180*math.pi
            mask_flag = data["IsOutsideMask"]
            healpix_mask = hp.ang2pix(nside,mask_dec,mask_ra)
            hdu.close()
            ######################################################################################

            #端のpixelを選ぶ
            ra_edge = (np.array(Mask.area(np.array(tract).astype(np.int64)))[0]%360)/180*math.pi
            dec_edge = math.pi/2-np.array(Mask.area(np.array(tract).astype(np.int64)))[1]/180*math.pi
            ra_ = np.linspace(np.min(ra_edge),np.max(ra_edge),100)
            dec_ = np.linspace(np.min(dec_edge),np.max(dec_edge),100)
            dd,rr = np.meshgrid(dec_,ra_)
            r = rr.flatten()
            d = dd.flatten()
            #print(r)

            r1 = ((sorted(ra_edge)[-2]-0.0001 < r) & (r < np.max(r) + 0.0001))
            r2 = ((sorted(ra_edge)[1]+0.0001 > r) & (r > np.min(r)-0.0001))
            d1 = ((sorted(dec_edge)[-2]-0.0001 < d) & (d < np.max(d) + 0.0001))
            d2 = ((sorted(dec_edge)[1] + 0.0001 > d) & (d > np.min(d) - 0.0001))

            #print(d[((sorted(ra_edge)[-2]-0.0002 < r) & (r < np.max(r) + 0.0002))|((sorted(dec_edge)[-2] < d) & (d < np.max(d)))
                       #|((sorted(ra_edge)[1]+0.0002 > r) & (r > np.min(r)-0.0002))|((sorted(dec_edge)[1] > d) & (d > np.min(d)))])
            healpix_edge = hp.ang2pix(nside, d[(r1|r2|d1|d2)], r[(r1|d1|r2|d2)])
            if (- np.min(ra_edge) + np.max(ra_edge)>3):
                ra_ = np.linspace(np.max(ra_edge) - 2*math.pi,np.min(ra_edge),100)
                dec_ = np.linspace(np.min(dec_edge),np.max(dec_edge),100)
                dd,rr = np.meshgrid(dec_,ra_)
                r = rr.flatten()
                d = dd.flatten()

                r1 = ((sorted(ra_edge)[1]-0.0001 < r) & (r < np.max(r) + 0.0001))
                r2 = ((sorted(ra_edge)[-2]-2*math.pi+0.0001 > r) & (r > np.min(r)-0.0001))
                d1 = ((sorted(dec_edge)[-2]-0.0001 < d) & (d < np.max(d) + 0.0001))
                d2 = ((sorted(dec_edge)[1] + 0.0001 > d) & (d > np.min(d) - 0.0001))

            #print(d[((sorted(ra_edge)[-2]-0.0002 < r) & (r < np.max(r) + 0.0002))|((sorted(dec_edge)[-2] < d) & (d < np.max(d)))
                       #|((sorted(ra_edge)[1]+0.0002 > r) & (r > np.min(r)-0.0002))|((sorted(dec_edge)[1] > d) & (d > np.min(d)))])
                r = r+(r<0)*2*math.pi
        #print(r[(r1|d1|r2|d2)])
                healpix_edge = hp.ang2pix(nside, d[(r1|r2|d1|d2)], r[(r1|d1|r2|d2)])
        
            for pixel in healpix_range:
                #effective areaを求める################################
                flag = mask_flag[healpix_mask == pixel]
                eff = np.sum(flag)/len(flag)*area
                eff_area.append(eff)
                if (pixel in healpix_edge):
                    eff = 0
                ######################################################
                #selected galaxy density
                if eff>0:
                    dens = np.sum(healpix_target == pixel)/eff
                    density.append(dens)
                    healpix = np.append(healpix, pixel)
                #####################################################
                #photometry of each healpixe
                    for photo in photometry:
                        #property[photo].append(healpix_photometry(eff, tract, pixel, photo))
                        property[photo].append(healpix_photometry(all_pixel, eff, tract, pixel, photo))


    val = list(property.values())
    val.append(density)
    val.append(list(healpix))
    print(val)
    nam = list(property.keys())
    nam_all = tuple(nam + ["density" , "healpix"])
    print(nam_all)
    t = Table(val, names=nam_all)
    t.write("property_20.fits", format='fits')
    
    #total_area = np.sum(eff_area)

    fig = plt.figure(figsize = (16,12))
    i = 1
    for k,v in property.items():
        histogram(fig, np.array(density),k,v,i)
        i += 1

    fig.subplots_adjust(hspace=0.5, wspace=0.5)
    plt.tight_layout()
    pdf.savefig(fig)
    pdf.close()
    
    return
        
        


def histogram(fig, dens,key,value,i):
    ax = fig.add_subplot(2,6,i)
    if "seeing" in key:
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
        x = rangex["extinction"]
        ax.set_xlabel("E(B-V)")
        ax.set_title("dust extinction")
    elif key == "star":
        x = rangex["star"]
        ax.set_xlabel("stellar density [deg-2]")
        ax.set_title("stellar density")
    
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
    mean_density = np.array(mean_density)[number>10]
    SD_density = np.array(SD_density)[number>10]
    density_ave = np.average(mean_density)
    print(mean_density/density_ave)
    
    ax.errorbar(x1[number>10],mean_density/density_ave,yerr = SD_density/density_ave ,fmt = "o")
    ax.plot(x1[number>10], mean_density/density_ave)
    ax.hist(value, bins=20, range =(np.min(x1[number>10]),np.max(x1[number>10])), weights=[1/len(value)]*len(value))
    ax.set_ylim(0.0,2.5)

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

def healpix_photometry(all_pixel, eff, tract, healpix, photo):
    if "seeing" in photo:
        hdu = fits.open("database/s21-colorterm/tracts_all/%s_bias.fits" %tract)
        data = hdu[1].data
        if photo=="gseeing":
            g_size = data["gseeing"]
            photo_galaxy = g_size[all_pixel == healpix]
        elif photo=="rseeing":
            r_size = data["rseeing"]
            photo_galaxy = r_size[all_pixel == healpix]
        elif photo=="iseeing":
            i_size = data["iseeing"]
            photo_galaxy = i_size[all_pixel == healpix]
        elif photo=="zseeing":
            z_size = data["zseeing"]
            photo_galaxy = z_size[all_pixel == healpix]
        elif photo=="yseeing":
            y_size = data["yseeing"]
            photo_galaxy = y_size[all_pixel == healpix]

        a = np.average(photo_galaxy)
        hdu.close()

    elif "depth" in photo:
        hdu = fits.open("database/s21-colorterm/tracts_all/%s_bias.fits" %tract)
        data = hdu[1].data
        if photo == "gdepth":
            g_depth = data["g_depth"]
            photo_galaxy = g_depth[all_pixel == healpix]
        elif photo == "rdepth":
            r_depth = data["r_depth"]
            photo_galaxy = r_depth[all_pixel == healpix]
        elif photo == "idepth":
            i_depth = data["i_depth"]
            photo_galaxy = i_depth[all_pixel == healpix]
        elif photo == "zdepth":
            z_depth = data["z_depth"]
            photo_galaxy = z_depth[all_pixel == healpix]
        elif photo == "ydepth":
            y_depth = data["y_depth"]
            photo_galaxy = y_depth[all_pixel == healpix]
            
        a = np.average(photo_galaxy)
        hdu.close()
        
    elif (photo == "star"):
        a = np.sum(star_pixel == healpix)/eff
        f = open('few_star.txt', 'a')
        if (a<10000):
            f.write(str(healpix) + "," + str(a) + "¥n")
            f.close()

    elif (photo == "extinction"):
        hdu = fits.open('csfd_ebv.fits')
        data = hdu[1].data
        im=data["T"]
        m = np.ndarray.flatten(im)
        im = hp.ud_grade(m, nside)
        a = im[healpix]
    return a



if __name__ == '__main__':
    main()
