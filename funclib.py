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
import scipy


def sav_to_map(filename, field):
    """ Read .sav file data into a sunpy map.

    Parameters:
    ----------

    filename (string): Path to input data file (.sav format)
    field (string): Field of sav file to access

    Returns:
    -------

    data: SunPy map containing the data and arbitrary coordinate header
    """

    try:
        data = readsav(filename)
    except FileNotFoundError:
        raise FileNotFoundError('Cannot find '+filename)
    except Exception:
        raise Exception('Data does not appear to be in correct .sav format')

    if field not in data.keys():
        print('Field '+ field +' is not in file keys ', data.keys())
        sys.exit(1)

    fake_coord = SkyCoord(0*u.arcsec, 0*u.arcsec, obstime='2013-10-28 08:24',
                 observer='earth', frame=frames.Helioprojective)
    fake_header = sunpy.map.make_fitswcs_header(data=np.empty((512, 512)),
                                                coordinate=fake_coord,
                                                reference_pixel=[0, 0]*u.pixel,
                                                scale=[2, 2]*u.arcsec/u.pixel,
                                                telescope='Fake Telescope',
                                                instrument='Fake Instrument',
                                                wavelength=1000*u.angstrom)

    data_map = sunpy.map.Map(data[field], fake_header)
    
    return data_map

def sav_to_numpy(filename, instrument, field):
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

    try:
        data = readsav(filename)
    except FileNotFoundError:
        raise FileNotFoundError('Cannot find '+filename)
    except Exception:
        raise Exception('Data does not appear to be in correct .sav format')

    if instrument != 'IBIS':
        print('This functionality has so far only been implemented ' + 
              'for IBIS data.')
        sys.ext(1)
    
    if field not in data.keys():
        print('Field '+ field +' is not in file keys ', data.keys())
        sys.exit(1)

    # TO DO: catch error if file is not in sav format
    file = sio.readsav(filename)
    data = file[field]
    
    return data

def segment(data_map, skimage_method):
    """ Segment optical image of the solar photosphere into tri-value maps 
    with 0 = inter-granule, 0.5 = faculae, 1 = granule. 

    Parameters:
    ----------

    data_map (SunPy map): SunPy map containing data to segment
    skimage_method (string): skimage thresholding method - options are 'otsu', 
                             'li', 'isodata', 'mean', 'minimum', 'yen', 
                             'triangle'

    Returns:
    -------

    data_map (SunPy map): SunPy map containing segmentated image (with the
                          original header)
"""

    data = data_map.data
    header = data_map.meta

    # apply median filter
    median_filtered = scipy.ndimage.median_filter(data, size=3)

    # apply threshold
    threshold = get_threshold(median_filtered, skimage_method)

    # initial skimage segmentation
    segmented_image = np.uint8(median_filtered > threshold)  

    # fix the extra IGM bits in the middle of granules
    segmented_image_fixed = remove_middles(segmented_image)

    # mark faculae
    segmented_image_markfac = mark_faculae(segmented_image_fixed, data)

    # show pipeline process
    fig, ((ax0, ax1), (ax2, ax3)) = plt.subplots(2, 2, figsize=(30, 30))
    s1 = 20; s2 =26
    fig.suptitle('Intermediate processesing steps ', fontsize=s2)    
    im0 = ax0.imshow(data/np.max(data))
    ax0.set_title('scaled input image', fontsize=s1)
    plt.colorbar(im0, ax=ax0)
    im1 = ax1.imshow(segmented_image, cmap='gray')
    ax1.set_title('direct '+skimage_method+' skimage segmentation', fontsize=s1)
    plt.colorbar(im1, ax=ax1)
    im2 = ax2.imshow(segmented_image_fixed, cmap='gray')
    ax2.set_title('wrong middles removed', fontsize=s1)
    plt.colorbar(im2, ax=ax2)
    im3 = ax3.imshow(segmented_image_markfac, cmap='gray')
    ax3.set_title('faculae identified', fontsize=s1)
    plt.colorbar(im3, ax=ax3)
    plt.axis('off')
    plt.savefig('intermediate_outputs.png')

    # convert segmentated image back into SunPy map with original header
    segmented_map = sunpy.map.Map(segmented_image, header) 

    return segmented_map

def get_threshold(data, method):
    """ Get the threshold value using given skimage segmentation type.

        Parameters:
        ----------

        data (numpy array): data to threshold
        method (string): skimage thresholding method - options are 'otsu', 
                        'li', 'isodata', 'mean', 'minimum', 'yen', 'triangle'

        Returns:
        -------

        threshold (float): threshold
    """    

    methods = ['otsu', 'li', 'isodata', 'mean', 'minimum', 'yen', 'triangle']
    if method not in methods:
        print('Method must be one of: ' + str(methods))
        sys.exit(1)
    if method == 'otsu': # works ok, but classifies low-flux ganules as IG
        threshold = skimage.filters.threshold_otsu(data)
    elif method == 'li': # slightly better than otsu
        threshold = skimage.filters.threshold_li(data)
    elif method == 'yen': # poor - overidentifies IG
        threshold = skimage.filters.threshold_yen(data)
    elif method == 'mean': # similar to li
        threshold = skimage.filters.threshold_mean(data)
    elif method == 'minimum': # terrible - identifies almost all as granule
        threshold = skimage.filters.threshold_minimum(data)
    elif method == 'triangle': # does not work well - overidentifies IG worse than yen
        threshold = skimage.filters.threshold_triangle(data)
    elif method == 'isodata': # similar to otsu
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

    if len(np.unique(segmented_image)) > 2:
        print('segmented_image must have only values of 1 and 0')
        sys.exit(1)

    segmented_image = segmented_image #[80:90, 130:140]/255 
    segmented_image_fixed = np.copy(segmented_image)
    labeled_seg = skimage.measure.label(segmented_image+1, connectivity=2)
    values = np.unique(labeled_seg)
    for value in values:
        mask = np.zeros_like(segmented_image)
        mask[labeled_seg==value] = 1
        if np.sum(np.multiply(mask, segmented_image)) == 0:  # check that is a 0 (black) region
            if len(segmented_image_fixed[mask==1]) < 100: # check that region is small [bad criteria, change later on] 
                segmented_image_fixed[mask==1] = 1

    return segmented_image_fixed
        
  
def mark_faculae(segmented_image, data):
    """ Mark faculae seperatly from granules - give them a value of 0.5 not 1

        Parameters:
        ----------

        segmented_image (numpy array): the segmented image containing incorrect middles

        Returns:
        -------

        segmented_image_fixed (numpy array): the segmented image with faculae marked as 0.5
    """   

    if len(np.unique(segmented_image)) > 2:
        print('segmented_image must have only values of 1 and 0')
        sys.exit(1)

    segmented_image = segmented_image#[80:90, 130:140]
    segmented_image_fixed = np.copy(segmented_image.astype(float))
    labeled_seg = skimage.measure.label(segmented_image+1, connectivity=2)
    values = np.unique(labeled_seg)
    for value in values:
        mask = np.zeros_like(segmented_image)
        mask[labeled_seg==value] = 1
        if np.sum(np.multiply(mask, segmented_image)) > 0:  # check that is a 1 (white) region
            region_size = len(segmented_image_fixed[mask==1])
            tot_flux = np.sum(data[mask==1])
            if region_size < 20: # check that region is small [bad criteria, change later on] 
                if tot_flux/region_size > 4000 : # check that avg flux very high
                    segmented_image_fixed[mask==1] = 0.5
        
    return segmented_image_fixed

