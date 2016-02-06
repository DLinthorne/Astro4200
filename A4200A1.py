
from __future__ import print_function
from photutils import daofind
from photutils import datasets
from astropy.stats import sigma_clipped_stats
import astropy as ap
from astropy.io import fits
import os
import matplotlib.pylab as plt
import numpy as np
from tqdm import tqdm
from numba import autojit, jit
from photutils import CircularAperture
from astropy.visualization import SqrtStretch
from astropy.visualization.mpl_normalize import ImageNormalize
from astropy import coordinates as coord
from astropy import units as u
from photutils import SkyCircularAperture as SCA
from astropy.coordinates import SkyCoord as SC
from photutils import aperture_photometry










''' DIRECTORY WHERE YOUR FITS FILES ARE! MAKE SURE ITS CORRECT :) '''
Directory = ('/home/nick/Desktop/School/Astro4200/MZCYG')

''' LIST OF ALL OF THE STARS IN THE DIRECTORY THAT ARE FITS FILES. '''
files = []
for file in os.listdir(Directory):
    if file.endswith('.fit'):
        files.append(file)
files.sort()




'''This  '''
Starnames = []
for i in range(len(files)):
    Starnames.append(files[i].split('.'))
    
    
    
'''Uses the command fits2bitmap in the system terminal to change all fits files into png images. 
   DO THIS ONLY ONE TIME, THEN COMMENT IT OUT!
   UNTESTED ON MACOS, if it works, OMG NICE. IF not, well theres probably a mac terminal command to do the
   same thing. In that case, replace the argument in the brackets with the mac command. 
'''
#for i in range(len(files)):
#    os.system('fits2bitmap {} -o {}'.format(Directory+'/'+files[i], Directory+'/'+Starnames[i][0]))




''' This line opens one of the fits files into memory. To observe its contituents, we need to use 
    the library functions to continue. '''
        
        
#Everyimage = []
#for i in range(0,2):
 
Star = fits.open('{}/{}'.format(Directory,files[55]))
stardata = Star[0].data
IMGheader = Star[0].header
pixelx, pixely = Star[0].header['CRPIX1'], Star[0].header['CRPIX2']
''' pixelx and pixely are the centers of the image, as array indices. THEY ARENT IN EVERY
    FITS FILE UNFORTUNATELY!!    I collected these values, since they are right next to the
    cepheid variable. 
    
    Issues with fits headers: Not consistent, some header extensions change names, or are not
    referenced in every single file. So we can't fully automate this whole process easily.'''


worldcoord = ap.wcs.WCS(header=Star[0].header)

#coordinates = worldcoord.wcs_pix2world(750,780,1)
px = np.arange(0,stardata.shape[0],1, dtype=int)
py = np.arange(0,stardata.shape[1],1, dtype=int)


pxtest = np.arange(0,500,1, dtype=int)
pytest = np.arange(0,500,1, dtype=int)

imagecoord = ap.wcs.WCS(header=IMGheader).wcs_pix2world(px,py,1)

'''Initializing the empty arrays for the coordframe function '''
RAarr = np.empty(shape=(1536,1536))
DECarr= np.empty(shape=(1536,1536))

'''This jit tag MIGHT speed up the function, though not by much because the function has to call
   an external library (ap.wcs.WCS) like millions of times to calculate array elements. '''
@jit(cache=True)
def coordframe(px,py,RAarr,DECarr):
    '''
        Takes the pixel indices of the image, and uses the header file to convert them into
        two arrays, one has the right ascension values for every pixel, and the other has 
        declination values for every pixel.
        
        BEWARE: As it is written now, these values are flipped upside down from the image.
        ie: The top left array element in the RA and DEC outputs are NOT the top left pixels in
        the image.
        
    '''
    for i in range(len(px)):
        for j in range(len(py)):
           RAarr[i,j] , DECarr[i,j] = worldcoord.wcs_pix2world(px[i],py[j],1)
    return RAarr, DECarr






'''Using the header info from the fits file, this line calculates the RA-DEC coordinates for each corner
   of the image, in DEGREES. '''
coords_for_corners = ap.wcs.WCS(header=Star[0].header).calc_footprint(header=Star[0].header)

'''This block takes the RA-DEC coordinates for the target object out of the header file,
   and converts them into floating-point numbers in DEGREES. 
   As well as the FWHM of the source in the image, for using when finding objects with
   aperture photometry.  '''

targetRA = repr(Star[0].header['RA'])
targetRA = targetRA.replace("'", "")
targetRA = np.array(targetRA.split(' '), dtype=float)
targetRAdeg = 15.0*(targetRA[0] + (targetRA[1]/60. )  +  (targetRA[2]/3600.)     )

targetDEC= repr(Star[0].header['DEC'])
targetDEC = targetDEC.replace('+', '')
targetDEC = targetDEC.replace("'", "")
targetDEC = np.array(targetDEC.split(' '), dtype=float)
targetDECdeg = (targetDEC[0] + (targetDEC[1]/60.) + (targetDEC[2]/3600.))

targetFWHM = float(repr(Star[0].header['FWHM']))






#print 'Pixel reference for Cepheid: x= {}, y= {}'.format(pixelx,pixely)


#subtract = stardata1+stardata
#plt.imshow(subtract,cmap='Spectral',origin='upper')


def add_images(stardata,added):
    stardata = stardata + added
    return stardata
'''
for i in range(40,45):
    intermed = fits.open('{}/{}'.format(Directory,files[i]))
    
    total = add_images(stardata,intermed[0].data)
'''

#Everyimage.append(stardata)
    
'''Doesn't work. ignore this function.'''
@jit(cache=True)    
def sourcefind(stardata):
    for i in range(0,1536):
        for j in range(0,1536):
            if float(stardata[j,i]) >= 8000.:
                print(j,i)


def cut_image(stardata,pixelx,pixely):
    '''
    Function to cut around my image, capturing at least the variable star in every image
    Taking into account that the image might have different resolutions.
    '''
    if stardata.shape[0] == 1024:
        xl , xr = int(pixelx-30), int(pixelx)
        yb, yt  = int(pixely-30), int(pixely)
        cutdata = stardata[xl:xr,yb:yt]
    elif stardata.shape[0] == 1536:
        xl , xr = int(pixelx-40), int(pixelx)
        yb, yt  = int(pixely), int(pixely+40)
        cutdata = stardata[yb:yt,xl:xr]
    
    return cutdata



cutdata = cut_image(stardata,pixelx,pixely)

    
#templist = []
#templist.append(sourcefind(stardata))


'''
pixel = np.arange(0,len(stardata), 1)


for i in tqdm(range(len(stardata)-2)):
        
    x = stardata[:,i]
    
    plt.plot(pixel, x, marker='.',c='b')

for i in tqdm(range(len(stardata) - 2)):
    xx = stardata[i,:]
    plt.plot(pixel,xx,marker='.', c='r')
'''


'''
    statistics on the image, determines background noise, standard deviation and mean photon counts
    to find all objects in the image that are above some threshold, and are NEAR some FWHM value.
    typing  "sources" without quotes in the spyde terminal will print all of the sources, x,y coords
    of their locations in the image, and a peak photon count, and corresponding flux and magnitude values 
    for the sources as well.
    
'''

mean, median, std = sigma_clipped_stats(cutdata, sigma =5.0, iters=10)

sources = daofind(cutdata - median, fwhm=targetFWHM, threshold=2.*std)

#sources = sources[(sources['peak'] > 1000.)]



#positions = SC(ra=targetRAdeg * u.deg, dec = targetDECdeg * u.deg, frame='fk5')

positions = (sources['xcentroid'], sources['ycentroid'])

apertures = CircularAperture(positions, r=4.0)

norm = ImageNormalize(stretch=SqrtStretch())


'''Uncomment these plot commands to view the image. plt.imshow plots the whole image, and apertures.plot
   plots circles around all of the things in the image that have been determined to be sources.'''
#plt.imshow(cutdata, cmap='Spectral', origin='lower', norm=norm)
#apertures.plot(color='blue',lw=1.5, alpha=0.5)


































    


