from __future__ import print_function
import astropy as ap
from astropy.io import fits
import os
import matplotlib.pyplot as plt
import numpy as np
from photutils import CircularAperture
from photutils import aperture_photometry, CircularAnnulus
from astropy.table import hstack


'''
    THIS CODE IS ABSOLUTELY NOT OPTIMIZED FOR ANYTHING, and currently messing with something will
    probably break it. Though a few seconds of looking at the error code will help ya fix, or change it 
    to a better method.
    
    If the object oriented approach doesn't simplify this enough, I suppose you could rip the code apart
    by deleting the 'class _____' crap, and remove 'self' or '@classmethod' from every function, then it 
    will work just like normal procedural code.



'''


class thissucks(object):
    
    def __init__(self, directory,sourceRA,sourceDEC):
        self.directory = directory
        self.files = self.files()
        self.Star = fits.open('{}/{}'.format(self.directory,self.files[0]))
        self.stardata = self.Star[0].data
        self.IMGheader = self.Star[0].header
        self.worldcoord = ap.wcs.WCS(header=self.IMGheader)
        self.sourceRA = sourceRA
        self.sourceDEC = sourceDEC


    def cut_vals(self,tbl):
        '''
            I don't bother using annoying astropy table objects, they are ridiculous and
            don't behave with python nicely. So this function cuts everything out that isnt 
            a string of digits in the table, to extract the background subtracted photon counts.
        '''
        nums = []
        for t in repr(tbl).split():
            try:
                nums.append(float(t))
            except:
                pass
        return nums


    def files(self):
        '''Load the files in your directory into an array... '''
        files = []
        for file in os.listdir(self.directory):
            if file.endswith('.fit'):
                files.append(file)
        files.sort()
        return files

    def calculate_fluxes(self,stardata,apertures,aper_annulus):
        '''
            Calculates the photon flux in an aperture, and an annulus
            around the aperture to subtract the background.
            
            As far as I can tell, the output value is still just a 'photon count', not technically a 
            photon flux. Possible modifications will start here, maybe uncommenting the apertures.area()
            which calculates the aperture phot_count divided by area of aperture.
            giving photons per area.
            
            I think we would further have to supplement that by dividing by the exposure time, and some 
            other wavelength value I cant think of, to get PHOTON FLUX ( photons per sec per cm^2 per wavelength)
        '''
        flux_table = aperture_photometry(stardata, apertures)
        bkg_table  = aperture_photometry(stardata, aper_annulus)
        
        phot_table = hstack([flux_table, bkg_table], table_names=['raw','bkg'])
        bkg_mean = phot_table['aperture_sum_bkg'] / aper_annulus.area()
        
        bkg_sum = bkg_mean * apertures.area()
        final_sum = phot_table['aperture_sum_raw'] - bkg_sum
        phot_table['residual_aperture_sum'] = final_sum
        #phot_table['res_aper_div_area'] = final_sum/apertures.area()
        #print(phot_table)
        return self.cut_vals(phot_table)
        

    def lightcurve(self,sourceRA,sourceDEC):
        '''
            Stores the background subtracted photon fluxes in an aperture around the sourceRA and sourceDEC.
            
            r_in :  inner radius of the annulus, in pixels
            r_out:  outer radius of the annulus, in pixels.
                        The r_in value should be strictly greater than the radius of the aperture.
                        I used an aperture radius of r=6 ( see below), this should be fine for all
                        of everyones stars, references or not.
           
            worldcoord : a world-coordinate-system object that uses the fits header stuff to transform
                         a location on the image in pixels into RA and DEC values in degrees, and vice versa.
                         This is how the sourceRA and sourceDEC 'know' exactly where there cepheids and references
                         are.
           
           
        '''
        
        starlist = []
        #testtable = np.empty(shape=(len(files),5))
        testrun = []
        for i in range(len(self.files)):
            Star = fits.open('{}/{}'.format(self.directory,self.files[i]))
            stardata = Star[0].data
            IMGheader = Star[0].header
            worldcoord = ap.wcs.WCS(IMGheader)
            exposure = Star[0].header['EXPTIME']
            starlist.append(Star[0].header['OBJECT'])
            zeromag = Star[0].header['ZMAG']
            juliandate = Star[0].header['JD']
            aper_annulus = CircularAnnulus((sourceRA, sourceDEC), r_in=7., r_out = 10.)
            apertures = CircularAperture((worldcoord.wcs_world2pix(sourceRA,sourceDEC,0)), r=6)
            
            testrun.append((self.calculate_fluxes(stardata,apertures, aper_annulus)[3], exposure,juliandate , zeromag))
        
        return testrun


    def mag_calc_noref(self,sourceRA,sourceDEC):
        '''
            Calculating the magnitude of the source star with NO reference star, using photon counts,
            exposure time, and zero point magnitude calibration: ALL header items.
            
            Output: 2 column array: 1st column is magnitudes, 2nd column is julian dates from fits files.
        '''
        testrun = self.lightcurve(sourceRA,sourceDEC)
        maglist = []
        for i in range(len(testrun)):
            magnitude = -2.5*np.log10((testrun[i][0]/testrun[i][1])) + testrun[i][3]
            maglist.append((magnitude, testrun[i][2]))
        
        maglist = np.array(maglist, dtype=float)
        return maglist
    
    def mag2(self,source,reference):
        '''
            Calculating the magnitude of the source star including a reference star.
            but no calibration constant is yet figured out.
        '''
        maglist = []    
        for i in range(len(source)):
            magnitude = -2.5*np.log10(source[i][0]/reference[i][0])
            maglist.append((magnitude, source[i][2]))
        return np.array(maglist,dtype=float)
    
    @classmethod
    def dms_to_deg(self,dmsval):
        ''' ENTER THE DEG/MIN/SEC value as a STRING.'''
        deg,mins,secs = dmsval.split(' ')    
        deg = float(deg)
        mins = float(mins)
        secs = float(secs)
        total = deg + (mins/60.) + (secs/3600.)
        return total

    @classmethod
    def hms_to_deg(self,hmsval):
        ''' ENTER THE DEG/MIN/SEC value as a STRING.'''
        hour,mins,secs = hmsval.split(' ')    
        hour = float(hour)
        mins = float(mins)
        secs = float(secs)
        total = 15.*(hour + (mins/60.) + (secs/3600.))
        return total


    def show_image(self,sourceRA,sourceDEC,refRA,refDEC):
        print('RED circle is the cepheid. WHITE circle is the reference object(s).')
        print('Add more reference stars by defining ref1aper = blahblahbelow, and ref1aper.plot(etc...)')
        #aper_annulus = CircularAnnulus((sourceRA, sourceDEC), r_in=6., r_out = 8.)

        apertures = CircularAperture((self.worldcoord.wcs_world2pix(sourceRA,sourceDEC,0)), r=6)
        ref1aper  = CircularAperture((self.worldcoord.wcs_world2pix(refRA,refDEC,0)),     r=6)
        #ref2aper  = CircularAperture((worldcoord.wcs_world2pix(ref2RA,ref2DEC,0)),     r=7)
        #ref3aper  = CircularAperture((worldcoord.wcs_world2pix(ref3RA,ref3DEC,0)),     r=3.5)
        #darkaper  = CircularAperture((worldcoord.wcs_world2pix(darkRA,darkDEC,0)),     r=3.5)
        
        
        fig = plt.figure()
        fig.add_subplot(111, projection = self.worldcoord)
        plt.imshow(self.stardata,origin='lower', cmap='Spectral')
        apertures.plot(color='red',lw=1.5, alpha=0.5)
        ref1aper.plot(color='white', lw=2.5, alpha=0.5)
        #apertures2.plot(color='white',lw=2.0,alpha=0.5)
        #largeaperture.plot(color='red',  lw=1.5, alpha=0.5)





#Since hms/dms_to_deg() functions are marked as class methods, we can actually call the class object without 
#initializing it, fancy object oriented toolbag tricks eh. I may be digressing too much already...
sourceRA = thissucks.hms_to_deg('21 21 54.73')
sourceDEC= thissucks.dms_to_deg('37 27 33.1')
ref1RA = 320.5825
ref1DEC= 37.55194445

#Directory where the pot of gold is:
directory = ('/home/nick/Desktop/School/Astro4200/MZCYG')

#initializing the object source, which has the methods needed to calculate lightcurves and whatever.
source = thissucks(directory,sourceRA,sourceDEC)



mystartest = source.mag_calc_noref(sourceRA,sourceDEC)

referencestar = source.mag_calc_noref(ref1RA,ref1DEC)



























































