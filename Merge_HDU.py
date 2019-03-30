from astropy.utils.data import get_pkg_data_filename
import astropy
from astropy.io import fits
from astropy.io.fits import getheader
import numpy as np

# ------------------------------------------------------------------------------------------------------------
# This program takes indiviual HDU in different fits files and puts them together in one individual fits file.
# ------------------------------------------------------------------------------------------------------------

# Get the file names:
fitsfile0 = get_pkg_data_filename('007_0.fits')
fitsfile1 = get_pkg_data_filename('007_1.fits')
fitsfile2 = get_pkg_data_filename('007_2.fits')

# Open the fits files, the first will be the primary hdu and the other ones will be the secondary HDUs
primary_hdu = fits.open(fitsfile0)
hdu1 = fits.open(fitsfile1)
hdu2 = fits.open(fitsfile2)

# Load the data
data0 = primary_hdu[0].data
data1 = hdu1[0].data
data2 = hdu2[0].data

#hduprim = fits.ImageHDU(data=data0)
hdu11 = fits.ImageHDU(data=data1)
hdu22 = fits.ImageHDU(data=data2)

# Get the header, in this case it comes from the primary fits file
hdr = fits.Header()
hdr = getheader('007_0.fits')
hduprim = fits.PrimaryHDU(data=data0,header=hdr)

# put all the HDUs together
hdul = fits.HDUList([hduprim,hdu11,hdu22])

# write it out
hdul.writeto('007_ready.fits')

