import os
import sys
import csv
import json
import time
import getpass
import argparse
import urllib.request, urllib.error, urllib.parse
import astropy.io.ascii as ascii
import re
import glob
from astropy.table import Table
import astropy.io.fits as fits
import numpy as np

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--tracts",  nargs ="*", help="specify tracts", default = '')
    
    global args
    args = parser.parse_args()
    tracts = args.tracts
    print("done")

    if (tracts == ''):
        tractname= 'fieldTractInfoS21.csv'
        tracts      =   ascii.read(tractname)['tract']
        tracts = list(tracts)

    for i in tracts:
        prefix = "mask/tract_%s*.fits"%(i)
        if not os.path.exists(prefix):
            ra,dec = area(i)
            if (max(ra)-min(ra)>300):
                file = "tract_%s_left.fits"%(i)
                command = "venice-4.0.3/bin/venice -r -xmin %s -xmax %s -ymin %s -ymax %s -coord spher -o tract_%s_left.fits"%(0, min(ra), min(dec), max(dec),i)
                os.system(command)

                command1 = "venice-4.0.3/bin/venice -m reg/masks_all.reg -cat tract_%s_left.fits -xcol ra -ycol dec -f all -flagName isOutsideMask -o tract_%s_flagged_left.fits"%(i,i)
                os.system(command1)
                
                command = "venice-4.0.3/bin/venice -r -xmin %s -xmax %s -ymin %s -ymax %s -coord spher -o tract_%s_right.fits"%(max(ra), 360, min(dec), max(dec),i)
                os.system(command)
                
                command1 = "venice-4.0.3/bin/venice -m reg/masks_all.reg -cat tract_%s_right.fits -xcol ra -ycol dec -f all -flagName isOutsideMask -o tract_%s_flagged_right.fits"%(i,i)
                os.system(command1)

            else:
                command = "venice-4.0.3/bin/venice -r -xmin %s -xmax %s -ymin %s -ymax %s -coord spher -o tract_%s.fits"%(min(ra), max(ra), min(dec), max(dec),i)
                os.system(command)
                command1 = "venice-4.0.3/bin/venice -m reg/masks_all.reg -cat tract_%s.fits -xcol ra -ycol dec -f all -flagName isOutsideMask -o tract_%s_flagged_left.fits"%(i,i)
                os.system(command1)
        else:
            print("already have mask")

    #os.system("mkdir -p mask")
    os.system("mv tract*.fits mask")

    return

def area(tract):
    ra = []
    dec = []
    filelist = glob.glob("tracts_patches_W-*.txt")
    #filelist = ["tracts_patches_W-autumn.txt", "tracts_patches_W-spring.txt"]
    line = "Tract: %s  Corner"%(tract)
    for filename in filelist:
        print(filename)
        f = open(filename,"r")
        while True:
            data = f.readline()
            if (data == ""):
                break
            if (line in data):
                a=re.split("[(:,)]",data)
                if (float(a[6])>360):
                    ra.append(float(a[6])-360)
                else:
                    ra.append(float(a[6]))
                dec.append(float(a[7]))
        f.close()
    
    print(ra)
    print(dec)
    return ra,dec

if __name__ == '__main__':
    main()
