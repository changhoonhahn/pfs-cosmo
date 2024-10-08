import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
import os
import glob
import pandas as pd
import tarfile

class select_image:
    def __init__(self,fits,filename,boolean):
        self.fits=fits
        self.name=filename
        self.bool=boolean
    def select(self):
        hdu=fits.open(self.fits)
        data=hdu[1].data
        ra=data["ra"]
        dec=data["dec"]
        d = {'ra':ra[self.bool][:30:] , 
             'dec': dec[self.bool][:30:]}
        df = pd.DataFrame(data=d)
        np.savetxt(self.name,df,header="ra, dec")
        
class show_image:
    def __init__(self,tarname,directory):
        self.tar=tarname
        self.dir=directory
        
    def image(self):
        with tarfile.open(self.tar+'.tar','r') as t:
            t.extractall(path=self.dir)
        filedir = "./"+self.dir+"/**"
        files = glob.glob(filedir+'/*')
        filenames = []
        
        for file in files:
            filename = os.path.splitext(os.path.basename(file))[0]
            filenames.append(".\\"+self.dir+"\\"+self.tar+"\\"+filename)
        
        fig=plt.figure(figsize=(50,50))
        for i in range (30):
            filename=filenames[i]+'.fits'
            data1=fits.getdata(filename,1)
            ax=fig.add_subplot(5,6,i+1)
            ax.imshow(data1,cmap='gray')
