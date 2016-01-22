
import astropy as ap
from astropy.io import fits
import os


''' This is the directory on MY computer, for your computer, change this to where your 
fits files are.  '''
Directory = ('/home/nick/Desktop/School/Astro4200/MZCYG')

''' LIST OF ALL OF THE STARS IN THE DIRECTORY THAT ARE FITS FILES. '''
files = []
for file in os.listdir(Directory):
    if file.endswith('.fit'):
        files.append(file)



''' This line opens one of the fits files into memory. To observe its contituents, we need to use 
    the library functions to continue. '''
Star = fits.open('{}/{}'.format(Directory,files[0]))



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





    


