import scipy.io as sio
from scipy.io import readsav
import astropy.units as u
from astropy.coordinates import SkyCoord
import sunpy.map
from sunpy.coordinates import frames
import matplotlib.pyplot as plt
import numpy as np
import sys
import skimage

# Function for out own use in converted .sav data to SunPy maps
def sav_to_maps(filename):
    """ Read .sav file data into a sunpy map.

    Parameters:
    ----------

    filename (string): Path to input data file (.sav format)

    Returns:
    -------

    flux_<...> (SunPy map): flux in specified wavelength
    velocity (SunPy map): measured velocity (negative is up)
    velocity_clean (SunPy map): measured velocity (negative is up) with 
                                large-scale structure removed
    """
    fake_coord = SkyCoord(0*u.arcsec, 0*u.arcsec, obstime='2013-10-28 08:24',
                 observer='earth', frame=frames.Helioprojective)
    fake_header = sunpy.map.make_fitswcs_header(data=np.empty((512, 512)), coordinate=fake_coord,
                                       reference_pixel=[0, 0]*u.pixel,
                                       scale=[2, 2]*u.arcsec/u.pixel,
                                       telescope='Fake Telescope', instrument='Fake Instrument',
                                       wavelength=1000*u.angstrom)

    data = readsav(filename)
    flux_4170 = sunpy.map.Map(data['rosa_wl'], fake_header)
    flux_4310 = sunpy.map.Map(data['rosa_gband'], fake_header)
    flux_7200 = sunpy.map.Map(data['ibis_wl'], fake_header)
    velocity =  sunpy.map.Map(data['velocity'], fake_header) # velocities in km/sec in 7090 Å spectral line
    velocity_clean =  sunpy.map.Map(data['velocity_cln'], fake_header)#  velocities after removal of the large scale structures (which are related more to oscillations than the granulation/convective overshoot)
    
    return flux_4170, flux_4310, flux_4170, flux_7200, velocity, velocity_clean

def sav_to_numpy(filename, instrument, band):
    """ Read .sav file data into a numpy array.

        Parameters:
        ----------

        filename (string): Path to input data file (.sav format)
        instrument (string): The telescope/instrument (ie. IBIS, DKSIST).
        band (string): Band or data field to read.

        Returns:
        -------

        data (numpy array): Array of data values.
    """
    if instrument != 'IBIS':
        print('This functionality has so far only been implemented ' + 
              'for IBIS data.')
        sys.ext(1)
    
    bands = ['rosa_wl', 'rosa_gband', 'ibis_wl', 'velocity', 'velocity_cln']
    if band not in bands:
        print('Band must be one of: ')
        print('rosa_wl      --> 4170 A')
        print('rosa_gband   --> 4310 A')
        print('ibis_wl      --> 7200 A')
        print('velocity     --> velocity in km/sec at 7090 Å')
        print('velocity_cln --> velocity with large-scale structure removed')
        sys.exit(1)

    # TO DO: catch error if file is not in sav format
    file = sio.readsav(filename)
    data = file[band]
    
    return data

def get_threshold(data, method):
    """ Get the threshold value using given skimage segmentation type.

        Parameters:
        ----------

        data (numpy array): data to threshold
        method (string): skimage thresholding method - 'otsu', 'li', 'isodata', 'mean', 'minimum', 'yen', 'triangle'

        Returns:
        -------

        threshold (float): threshold
    """    

    methods = ['otsu', 'li', 'isodata', 'mean', 'minimum', 'yen', 'triangle']
    if method not in methods:
        print('Method must be one of: ' + str(methods))
        sys.exit(1)
    if method == 'otsu':
        threshold = skimage.filters.threshold_otsu(data)
    elif method == 'li':
        threshold = skimage.filters.threshold_li(data)
    elif method == 'yen':
        threshold = skimage.filters.threshold_yen(data)
    elif method == 'mean':
        threshold = skimage.filters.threshold_mean(data)
    elif method == 'minimum':
        threshold = skimage.filters.threshold_minimum(data)
    elif method == 'triangle':
        threshold = skimage.filters.threshold_triangle(data)
    elif method == 'isodata':
        threshold = skimage.filters.threshold_isodata(data)

    return threshold

def remove_middles(segmented_image):
    """ Remove the erronous idenfication of intergranule material in the middle of granules
        that pure thresdhold segmentation produces

        Parameters:
        ----------

        segmented_image (numpy array): the segmented image containing incorrect middles

        Returns:
        -------

        segmented_image_fixed (numpy array): the segmented image without incorrect middles
    """   

    segmented_image = segmented_image/255 #[80:90, 130:140]/255 
    segmented_image_fixed = np.copy(segmented_image)
    labeled_seg = skimage.measure.label(segmented_image+1, connectivity=2)
    values = np.unique(labeled_seg)
    for value in values:
        mask = np.zeros_like(segmented_image)
        mask[labeled_seg==value] = 1
        if np.sum(np.multiply(mask, segmented_image)) == 0:  # check that is a 0 (black) region
            if np.sum(mask) < 100: # check that region is small [bad criteria, change later on] 
                segmented_image_fixed[mask==1] = 1

    return segmented_image_fixed
        
  

        

    
