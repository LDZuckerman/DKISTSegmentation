import glob
import os

import astropy.io.fits as fits
import astropy.units as u
import matplotlib.pyplot as plt
import matplotlib.lines as lines
import matplotlib.patheffects as mpe
import numpy as np
import scipy
import scipy.io as sio
import skimage
import sunpy
import sunpy.map
from astropy.coordinates import SkyCoord
from sunpy.coordinates import frames


def open_file(filename):
    """
    Convenience function for catching file read errors.
    ----------
    Parameters:
        filename: (string) Path to file.
    ----------
    Returns:
        Nothing
    """
    if not os.path.exists(filename):
        raise Exception('Input file ' + filename + ' could not be found.')

    if os.path.isdir(filename):
        raise Exception('Input file ' + filename + ' is a directory.')

    if not os.access(filename, os.R_OK):
        raise Exception('Input file ' + filename + ' is not readable.')

    if not (filename.endswith('.sav') | filename.endswith('.gz')):
        raise Exception('Input file must be a .sav file.')

    if os.path.getsize(filename) == 0:
        raise Exception('Input file must not be empty.')

    return


def save_to_fits(segmented_map, data_map, out_file, out_dir):
    """
    Save input sunpy map and output segmented sunpy map to fits file
    ----------
    Parameters:
        segmented_map (sunpy map): segmented data
        data_map (sunpy map): input data
        out_file (str): desired name of output fits file
        out_dir (str): filepath for the fits file
    ----------
    Returns:
        None: creates fits file with first extension being the
              segmented map and second extension being the input map
    """
    if not os.path.exists(out_dir):
        try:
            os.mkdir(out_dir)
        except Exception:
            raise OSError('Could not make directory ' + out_dir)
    try:
        filename = out_dir + out_file
    except Exception:
        raise TypeError('Appears that out_dir or out_file are not strings')

    try:
        segmented_map.save(filename, overwrite=True)
    except Exception:
        raise TypeError('segmented_map must be a sunpy map')

    fits.append(filename, data_map.data)


def sav_to_map(filename, field):
    """
    Read .sav file data into a sunpy map.
    ----------
    Parameters:
        filename (string): Path to input data file (.sav format)
        field (string): Field of sav file to access
    ----------
    Returns:
        data: SunPy map containing the data and arbitrary coordinate header
    """

    try:
        data = sio.readsav(filename)
    except FileNotFoundError:
        raise FileNotFoundError('Cannot find ' + filename)
    except Exception:
        raise Exception('Data ' + filename + ' does not appear to be in '
                        + 'correct .sav format')

    if field not in data.keys():
        raise Exception('Field ' + field +
                        ' is not in file keys ', data.keys())

    fake_coord = SkyCoord(0 * u.arcsec,
                          0 * u.arcsec,
                          obstime='2013-10-28 08:24',
                          observer='earth',
                          frame=frames.Helioprojective)
    fake_header = \
        sunpy.map.make_fitswcs_header(data=np.empty((512, 512)),
                                      coordinate=fake_coord,
                                      reference_pixel=[0, 0] * u.pixel,
                                      scale=[2, 2] * u.arcsec / u.pixel,
                                      telescope='Fake Telescope',
                                      instrument='Fake Instrument',
                                      wavelength=1000 * u.angstrom)

    data_map = sunpy.map.Map(data[field], fake_header)

    return data_map


def fits_to_map(filename):
    """
    Read .fits file data into a sunpy map.
    ----------
    Parameters:
        filename (string): Path to input data file (.fits format)
    ----------
    Returns:
        data: SunPy map containing the data and arbitrary coordinate header
    """

    try:
        hdu = fits.open(filename)
        data = hdu[0].data
    except FileNotFoundError:
        raise FileNotFoundError('Cannot find ' + filename)
    except Exception:
        raise Exception('Data does not appear to be in correct .fits format')

    fake_coord = SkyCoord(0 * u.arcsec,
                          0 * u.arcsec,
                          obstime='2013-10-28 08:24',
                          observer='earth',
                          frame=frames.Helioprojective)
    fake_header = \
        sunpy.map.make_fitswcs_header(data=np.empty((512, 512)),
                                      coordinate=fake_coord,
                                      reference_pixel=[0, 0] * u.pixel,
                                      scale=[2, 2] * u.arcsec / u.pixel,
                                      telescope='Fake Telescope',
                                      instrument='Fake Instrument',
                                      wavelength=1000 * u.angstrom)

    data_map = sunpy.map.Map(data, fake_header)

    return data_map


def sav_to_numpy(filename, instrument, field):
    """
    Read .sav file data into a numpy array.
    ----------
    Parameters:
        filename (string): Path to input data file (.sav format)
        instrument (string): The telescope/instrument (ie. IBIS, DKSIST).
        field (string): Band or data field to read.
    ----------
    Returns:
        data (numpy array): Array of data values.
    """

    try:
        data = sio.readsav(filename)
    except FileNotFoundError:
        raise FileNotFoundError('Cannot find ' + filename)
    except Exception:
        raise Exception('Data does not appear to be in correct .sav format')

    if instrument not in ['IBIS']:
        raise Exception('This functionality has so far only been ' +
                        'implemented for IBIS .sav data.')

    if field not in data.keys():
        raise Exception('Field ' + field + ' is not in file keys ',
                        data.keys())

    # TO DO: catch error if file is not in sav format
    file = sio.readsav(filename)
    data = file[field]

    return data


def segment(data_map, skimage_method, plot_intermed=True, out_dir='output/',
            res='DKIST'):
    """
    Segment optical image of the solar photosphere into tri-value maps
    with 0 = inter-granule, 0.5 = faculae, 1 = granule.
    ----------
    Parameters:
        data_map (SunPy map): SunPy map containing data to segment
        skimage_method (string): skimage thresholding method -
                                 options are 'otsu', 'li', 'isodata',
                                  'mean', 'minimum', 'yen', 'triangle'
        plot_intermed (True or False): whether to intermediate data product
                                       image
        out_dir (str): Desired directory in which to save intermediate data
                                  product image (if plot_intermed = True)
    ----------
    Returns:
        data_map (SunPy map): SunPy map containing segmentated image (with the
                              original header)
    """

    if type(data_map) != sunpy.map.mapbase.GenericMap:
        raise TypeError('Input must be sunpy map.')

    methods = ['otsu', 'li', 'isodata', 'mean', 'minimum', 'yen', 'triangle']
    if skimage_method not in methods:
        raise TypeError('Method must be one of: ' + str(methods))

    data = data_map.data
    header = data_map.meta

    # apply median filter
    median_filtered = scipy.ndimage.median_filter(data, size=3)

    # apply threshold
    threshold = get_threshold(median_filtered, skimage_method)

    # initial skimage segmentation
    segmented_image = np.uint8(median_filtered > threshold)

    # fix the extra IGM bits in the middle of granules
    segmented_image_fixed = trim_interganules(segmented_image)

    # mark faculae
    segmented_image_markfac = mark_faculae(segmented_image_fixed, data, res)

    if plot_intermed:
        # show pipeline process
        fig, ((ax0, ax1), (ax2, ax3)) = plt.subplots(2, 2, figsize=(15, 15))
        s1 = 18
        s2 = 24
        fig.suptitle('Intermediate Processesing Steps \n', fontsize=s2)

        im0 = ax0.imshow(data / np.max(data), origin='lower')
        ax0.set_title('scaled Intensity Data', fontsize=s1)
        plt.colorbar(im0, ax=ax0, shrink=0.8)

        im1 = ax1.imshow(segmented_image, cmap='gray', origin='lower')
        ax1.set_title('Initial Thresholding', fontsize=s1)
        plt.colorbar(im1, ax=ax1, shrink=0.8)

        im2 = ax2.imshow(segmented_image_fixed, cmap='gray', origin='lower')
        ax2.set_title('Extraneous IG Material Removed', fontsize=s1)
        plt.colorbar(im2, ax=ax2, shrink=0.8)

        im3 = ax3.imshow(segmented_image_markfac, cmap='gray', origin='lower')
        ax3.set_title('Faculae Identified', fontsize=s1)
        plt.colorbar(im3, ax=ax3, shrink=0.8)

        plt.tight_layout()

        outline = mpe.withStroke(linewidth=5, foreground='black')
        custom_lines = [lines.Line2D([0], [0], color='white', lw=4,
                        path_effects=[outline]),
                        lines.Line2D([0], [0], color='black', lw=4),
                        lines.Line2D([0], [0], color='grey', lw=4)]
        legax = plt.axes([0.1, 0.1, 0.8, 0.85], alpha=0)
        legax.axis('off')
        legax.legend(custom_lines, ['Granule', 'Intergranule', 'Faculae'],
                     loc='upper center', ncol=3, fontsize='x-large')

        if not os.path.exists(out_dir):
            try:
                os.mkdir(out_dir)
            except Exception:
                raise OSError('Could not make directory ' + out_dir)
        plt.savefig(out_dir + 'intermediate_outputs.png')

    # convert segmentated image back into SunPy map with original header
    segmented_map = sunpy.map.Map(segmented_image, header)

    return segmented_map


def get_threshold(data, method):
    """
    Get the threshold value using given skimage segmentation type.
    ----------
    Parameters:
        data (numpy array): data to threshold
        method (string): skimage thresholding method - options are 'otsu',
                        'li', 'isodata', 'mean', 'minimum', 'yen', 'triangle'
    ----------
    Returns:
        threshold (float): threshold
    """

    if not type(data) == np.ndarray:
        raise ValueError('Input data must be an array.')

    methods = ['otsu', 'li', 'isodata', 'mean', 'minimum', 'yen', 'triangle']
    if method not in methods:
        raise ValueError('Method must be one of: ' + str(methods))
    if method == 'otsu':  # works ok, but classifies low-flux ganules as IG
        threshold = skimage.filters.threshold_otsu(data)
    elif method == 'li':  # slightly better than otsu
        threshold = skimage.filters.threshold_li(data)
    elif method == 'yen':  # poor - overidentifies IG
        threshold = skimage.filters.threshold_yen(data)
    elif method == 'mean':  # similar to li
        threshold = skimage.filters.threshold_mean(data)
    elif method == 'minimum':  # terrible - identifies almost all as granule
        threshold = skimage.filters.threshold_minimum(data)
    elif method == 'triangle':  # overidentifies IG worse than yen
        threshold = skimage.filters.threshold_triangle(data)
    elif method == 'isodata':  # similar to otsu
        threshold = skimage.filters.threshold_isodata(data)

    return threshold


def trim_interganules(segmented_image):
    """
    Remove the erronous idenfication of intergranule material in the
    middle of granules that pure threshold segmentation produces.
    ----------
    Parameters:
        segmented_image (numpy array): the segmented image containing
                                       incorrect middles
    ----------
    Returns:
        segmented_image_fixed (numpy array): the segmented image without
                                             incorrect middles
    """

    if len(np.unique(segmented_image)) > 2:
        raise ValueError('segmented_image must have only values of 1 and 0')

    segmented_image_fixed = np.copy(segmented_image)
    labeled_seg = skimage.measure.label(segmented_image + 1, connectivity=2)
    values = np.unique(labeled_seg)
    # find value of the 0 region that is big continuous region
    size = 0
    for value in values:
        if len((labeled_seg[labeled_seg == value])) > size:
            real_IG_value = value
            size = len(labeled_seg[labeled_seg == value])
    # set all other 0 regions to 1
    for value in values:
        if value != real_IG_value:
            segmented_image_fixed[labeled_seg == value] = 1

    return segmented_image_fixed


def mark_faculae(segmented_image, data, res):
    """
    Mark faculae seperatly from granules - give them a value of 0.5 not 1
    ----------
    Parameters:
        data (numpy array): the original flux values
        segmented_image (numpy array): the segmented image containing
                                       incorrect middles
    ----------
    Returns:
        segmented_image_fixed (numpy array): the segmented image with faculae
                                             marked as 0.5
    """

    if res == 'DKIST':
        fac_size_limit = 250  # number of pixels criterion for faculae
        fac_brightness_limit = 5000  # flux/pix criterion for faculae
    if res == 'IBIS':
        fac_size_limit = 20  # number of pixels criterion for faculae
        fac_brightness_limit = 4000  # flux/pix criterion for faculae

    if len(np.unique(segmented_image)) > 2:
        raise ValueError('segmented_image must have only values of 1 and 0')

    segmented_image = segmented_image
    segmented_image_fixed = np.copy(segmented_image.astype(float))
    labeled_seg = skimage.measure.label(segmented_image + 1, connectivity=2)
    values = np.unique(labeled_seg)
    for value in values:
        mask = np.zeros_like(segmented_image)
        mask[labeled_seg == value] = 1
        # check that is a 1 (white) region
        if np.sum(np.multiply(mask, segmented_image)) > 0:
            region_size = len(segmented_image_fixed[mask == 1])
            tot_flux = np.sum(data[mask == 1])
            # check that region is small
            if region_size < fac_size_limit:
                # check that avg flux very high
                if tot_flux / region_size > fac_brightness_limit:
                    segmented_image_fixed[mask == 1] = 0.5

    return segmented_image_fixed


def find_data(filepath):
    """
    Given a master filepath, traverses through all embedded directories
    and files, and returns a list of all the files that end in
    .fits or .sav.

    Useful for raw data, which is usually in nested directory structures.
    ----------
    Parameters:
        filepath (string): the filepath to search for .fits or .sav files
    ----------
    Returns:
        files_to_be_segmented (list of strings): the list of files
                                                 to be segmented.
    """
    files_to_be_segmented = []
    files = glob.glob(filepath + '**', recursive=True)
    for file in files:
        if file.endswith('.fits') or file.endswith('.sav'):
            files_to_be_segmented.append(file)
    return files_to_be_segmented
