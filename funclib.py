import glob
import os
import sys
import astropy.io.fits as fits
import astropy.units as u
from astropy.coordinates import SkyCoord
import matplotlib.pyplot as plt
import matplotlib.lines as lines
import matplotlib.patheffects as mpe
import matplotlib.colors as colors
import numpy as np
import scipy.ndimage as sndi
import scipy.io as sio
import skimage
import sunpy
import sunpy.map
from sunpy.coordinates import frames
from sklearn.cluster import KMeans as KMeans
import warnings
from erfa import ErfaWarning


def save_to_fits(segmented_map,
                 data_map,
                 out_file,
                 out_dir,
                 header,
                 confidence):
    """
    Save input sunpy map and output segmented sunpy map to fits file.
    ----------
    Parameters:
        segmented_map (sunpy map): segmented data
        data_map (sunpy map): input data
        out_file (str): desired name of output fits file
        out_dir (str): filepath for the fits file
        header (fits header object or None): header to use for output fits file
        confidence (float): fraction out of 1 returned by cross
                            correlation function (1 = perfect agreement,
                            0 = no agreement).
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

    if not type(data_map) == sunpy.map.mapbase.GenericMap:
        raise TypeError('data_map must be a sunpy map')

    if not type(segmented_map) == sunpy.map.mapbase.GenericMap:
        raise TypeError('segmented_map must be a sunpy map')

    try:
        if header is not None:
            seg_hdu = fits.PrimaryHDU(segmented_map.data, header)
        else:
            seg_hdu = fits.PrimaryHDU(segmented_map.data)
        seg_hdu.header.append(('CONFIDEN', confidence))
        raw_hdu = fits.ImageHDU(data_map.data)
        hdu = fits.HDUList([seg_hdu, raw_hdu])
    except Exception:
        raise TypeError('Segmented_map must be a sunpy map')

    hdu.writeto(filename, overwrite=True)

    return None


def sav_to_map(filename, field):
    """
    Read .sav file data into a sunpy map.
    ----------
    Parameters:
        filename (string): Path to input data file (.sav format)
        field (string): Field of sav file to access
    ----------
    Returns:
        data: SunPy map containing the data and generic header
    """

    # Suppress harmless warning about 'bad' placeholder date in
    # placeholder header (see comment below)
    warnings.filterwarnings(action='ignore', category=ErfaWarning)

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

    # Sunpy absolutely requires that a header be present in all sunpy.Map
    # objects. Becuase .sav files do not contain header information, create
    # a generic header as a placeholder under the strict requirements
    # imposed by sunpy on header field types and formats. No empty or
    # None fields are permitted by sunpy.Map creation function.
    print('WARNING: .sav input file contains no header; generating ' +
          'placeholder header in sunpy.Map object.')
    coord = SkyCoord(np.nan * u.arcsec,
                     np.nan * u.arcsec,
                     obstime='1111-11-11 11:11',
                     observer='earth',
                     frame=frames.Helioprojective)
    header = \
        sunpy.map.make_fitswcs_header(data=np.empty((0, 0)),
                                      coordinate=coord,
                                      reference_pixel=[np.nan, np.nan]
                                      * u.pixel,
                                      scale=[np.nan, np.nan]
                                      * u.arcsec / u.pixel,
                                      telescope='Unknown',
                                      instrument='Unknown',
                                      wavelength=np.nan * u.angstrom)
    data_map = sunpy.map.Map(data[field], header)

    return data_map


def fits_to_map(filename):
    """
    Read .fits file data into a sunpy map.
    ----------
    Parameters:
        filename (string): Path to input data file (.fits format)
    ----------
    Returns:
        data_map: SunPy map containing the data and header
    """

    try:
        hdu = fits.open(filename)
        data = hdu[0].data
    except FileNotFoundError:
        raise FileNotFoundError('Cannot find ' + filename)
    except Exception:
        raise Exception('Data does not appear to be in correct .fits format')

    data_map = sunpy.map.Map(filename)

    return data_map


def sav_to_numpy(filename, instrument, field):
    """
    Read .sav file data into a numpy array.
    ----------
    Parameters:
        filename (string): Path to input data file (.sav format)
        instrument (string): The telescope/instrument (ie. IBIS, DKIST).
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

    file = sio.readsav(filename)
    data = file[field]

    return data


def segment(file_id, data_map, skimage_method, res, plot_intermed=True,
            mark_dim_centers=False, out_dir='output/'):

    """
    Segment optical image of the solar photosphere into tri-value maps
    with 0 = intergranule, 0.5 = faculae, 1 = granule.
    ----------
    Parameters:
        file_id (string): id name of file to be segmented
        data_map (SunPy map): SunPy map containing data to segment
        skimage_method (string): skimage thresholding method -
                                options are 'otsu', 'li', 'isodata',
                                'mean', 'minimum', 'yen', 'triangle'
        plot_intermed (True or False): whether to plot and save intermediate
                                data product image
        mark_dim_centers (True or False): whether to mark dim granule centers
                                as a seperate catagory for future exploration
        out_dir (str): Desired directory in which to save intermediate data
                                product image (if plot_intermed = True);
        res (float): Spatial resolution (arcsec/pixel) of the data
    ----------
    Returns:
        segmented_map (SunPy map): SunPy map containing segmentated image
                              (with the original header)
    """

    if type(data_map) != sunpy.map.mapbase.GenericMap:
        raise TypeError('Input must be sunpy map.')

    methods = ['otsu', 'li', 'isodata', 'mean', 'minimum', 'yen', 'triangle']
    if skimage_method not in methods:
        raise TypeError('Method must be one of: ' + str(methods))

    data = data_map.data
    header = data_map.meta

    # apply median filter
    median_filtered = sndi.median_filter(data, size=3)

    # apply threshold
    threshold = get_threshold(median_filtered, skimage_method)

    # initial skimage segmentation
    segmented_image = np.uint8(median_filtered > threshold)

    # fix the extra IGM bits in the middle of granules
    segmented_image_fixed = trim_intergranules(segmented_image,
                                               mark=mark_dim_centers)

    # mark faculae
    segmented_image_markfac = mark_faculae(segmented_image_fixed, data, res)

    if plot_intermed:
        # show pipeline process
        fig, ((ax0, ax1), (ax2, ax3)) = plt.subplots(2, 2, figsize=(14, 13))
        s1 = 16
        s2 = 22
        fig.suptitle('Intermediate Processesing Steps \n', fontsize=s2)

        # define colormap to bring out faculae and dim middles
        col_dict = {0: "black",
                    0.5: "blue",
                    1: "white",
                    1.5: "#ffc406"}
        cmap = colors.ListedColormap([col_dict[x] for x in col_dict.keys()])
        norm_bins = np.sort([*col_dict.keys()]) + 0.5
        norm_bins = np.insert(norm_bins, 0, np.min(norm_bins) - 1.0)
        norm = colors.BoundaryNorm(norm_bins, 4, clip=True)

        im0 = ax0.imshow(data / np.max(data), origin='lower')
        ax0.set_title('Scaled Intensity Data', fontsize=s1)
        plt.colorbar(im0, ax=ax0, shrink=0.8)

        im1 = ax1.imshow(segmented_image, norm=norm, cmap=cmap,
                         interpolation='none', origin='lower')
        ax1.set_title('Initial Thresholding', fontsize=s1)

        im2 = ax2.imshow(segmented_image_fixed, norm=norm, cmap=cmap,
                         interpolation='none', origin='lower')
        if mark_dim_centers:
            ax2.set_title('Dim IG Material Marked', fontsize=s1)
        if not mark_dim_centers:
            ax2.set_title('Extraneous IG Material Removed', fontsize=s1)

        im3 = ax3.imshow(segmented_image_markfac, norm=norm, cmap=cmap,
                         interpolation='none', origin='lower')
        ax3.set_title('Faculae Identified', fontsize=s1)

        plt.tight_layout()

        # rescale axis
        l0, b0, w0, h0 = ax0.get_position().bounds
        newpos = [l0, b0-0.01, w0, h0]
        ax0.set_position(newpos)
        l1, b1, w1, h1 = ax1.get_position().bounds
        newpos = [l1, b0-0.01, w0, h0]
        ax1.set_position(newpos)
        l2, b2, w2, h2 = ax2.get_position().bounds
        newpos = [l0, b2, w0, h0]
        ax2.set_position(newpos)
        l3, b3, w3, h3 = ax3.get_position().bounds
        newpos = [l3, b3, w0, h0]
        ax3.set_position(newpos)

        # add color bar at top
        outline = mpe.withStroke(linewidth=5, foreground='black')
        legax = plt.axes([0.1, 0.1, 0.8, 0.85], alpha=0)
        legax.axis('off')
        if mark_dim_centers:
            labels = ['Granule', 'Intergranule', 'Faculae', 'Dim Centers']
            custom_lines = [lines.Line2D([0], [0], color='white', lw=4,
                                         path_effects=[outline]),
                            lines.Line2D([0], [0], color='black', lw=4),
                            lines.Line2D([0], [0], color="#ffc406", lw=4),
                            lines.Line2D([0], [0], color='blue', lw=4)]
            ncol = 4
        if not mark_dim_centers:
            labels = ['Granule', 'Intergranule', 'Faculae']
            custom_lines = [lines.Line2D([0], [0], color='white', lw=4,
                                         path_effects=[outline]),
                            lines.Line2D([0], [0], color='black', lw=4),
                            lines.Line2D([0], [0], color="#ffc406", lw=4)]
            ncol = 3
        legax.legend(custom_lines, labels, loc='upper center', ncol=ncol,
                     fontsize='x-large')

        if not os.path.exists(out_dir):
            try:
                os.mkdir(out_dir)
            except Exception:
                raise OSError('Could not make directory ' + out_dir)

        plt.savefig(out_dir + 'segmentation_plots_' + file_id + '.png')

    # convert segmentated image back into SunPy map with original header
    segmented_map = sunpy.map.Map(segmented_image_markfac, header)

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


def trim_intergranules(segmented_image, mark=False):
    """
    Remove the erronous idenfication of intergranule material in the
    middle of granules that pure threshold segmentation produces.
    ----------
    Parameters:
        segmented_image (numpy array): the segmented image containing
                                       incorrect extra intergranules
        mark (bool): if False, remove erronous intergranules. If True,
                     mark them as 0.5 instead (for later examination).
    ----------
    Returns:
        segmented_image_fixed (numpy array): the segmented image without
                                             incorrect extra intergranules
    """

    if len(np.unique(segmented_image)) > 2:
        raise ValueError('segmented_image must have only values of 1 and 0')

    segmented_image_fixed = np.copy(segmented_image).astype(float)
    labeled_seg = skimage.measure.label(segmented_image + 1, connectivity=2)
    values = np.unique(labeled_seg)
    # find value of the large continuous 0-valued region
    size = 0
    for value in values:
        if len((labeled_seg[labeled_seg == value])) > size:
            real_IG_value = value
            size = len(labeled_seg[labeled_seg == value])

    # set all other 0 regions to mark value (1 or 0.5)
    for value in values:
        if np.sum(segmented_image[labeled_seg == value]) == 0:
            if value != real_IG_value:
                if not mark:
                    segmented_image_fixed[labeled_seg == value] = 1
                elif mark:
                    segmented_image_fixed[labeled_seg == value] = 0.5

    return segmented_image_fixed


def mark_faculae(segmented_image, data, res):
    """
    Mark faculae seperatley from granules - give them a value of 2 not 1.
    ----------
    Parameters:
        segmented_image (numpy array): the segmented image containing
                                incorrect middles
        data (numpy array): the original flux values
        res (float): Spatial resolution (arcsec/pixel) of the data
    ----------
    Returns:
        segmented_image_fixed (numpy array): the segmented image with faculae
                                             marked as 1.5
    """

    fac_size_limit = 2  # max size of a faculae in sqaure arcsec
    fac_pix_limit = fac_size_limit / res
    fac_brightness_limit = np.mean(data) + 0.5 * np.std(data)

    if len(np.unique(segmented_image)) > 3:
        raise ValueError('segmented_image must have only values of 1, 0, ' +
                         'an 0.5 (if dim centers marked)')

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
            if region_size < fac_pix_limit:
                # check that avg flux very high
                if tot_flux / region_size > fac_brightness_limit:
                    segmented_image_fixed[mask == 1] = 1.5

    return segmented_image_fixed


def find_data(filepath):
    """
    Given a master filepath, traverse through all embedded directories
    and files, and return a list of all the files that end in
    .fits or .sav.

    Useful for raw data, assumed to be in nested directory structures.
    ----------
    Parameters:
        filepath (string): the filepath to search for .fits or .sav files
    ----------
    Returns:
        files_to_be_segmented (list of strings): the list of files
                                                 to be segmented
    """

    if not os.path.isdir(filepath):
        raise Exception(filepath + ' is not a directory.')

    files_to_be_segmented = []
    files = os.listdir(filepath)
    for file in files:
        if file.endswith('.fits') or file.endswith('.sav'):
            files_to_be_segmented.append(file)
    if not files_to_be_segmented:
        raise OSError('No data found in ' + filepath)

    return files_to_be_segmented


def overplot_velocities(seg_map, input_file, output_path):
    """
    Overplot segmentation contours on velocity data as a check on accuracy.
    ----------
    Parameters:
        seg_map (sunpy map): Sunpy map containing segmented data
        input_file (sav data file): IBIS map to be segmented and overplotted
        output_path (string): path to save output plot
    ----------
    Returns:
        None; saves outplot plot
    """

    # Suppress harmless warning from plt.contour
    warnings.filterwarnings(action='ignore', category=UserWarning)

    if input_file.endswith('.sav'):
        file = sio.readsav(input_file)
        velocities = file['velocity_cln']
    else:
        raise Exception('Velocity data only reported for IBIS datasets.')

    segmented_data = seg_map.data

    plt.figure(figsize=[8, 8])
    plt.imshow(velocities)
    plt.colorbar(label='Normalized radial velocity')
    plt.contour(segmented_data, [2], colors='black')
    ax = plt.gca()
    # the next two lines add a legend entry for the contours
    patches = [lines.Line2D([0], [0],
               color='black',
               label='Segmentation \n contours',
               linewidth=1)]
    plt.legend(handles=patches,
               bbox_to_anchor=(0.65, 0.10),
               loc=2, borderaxespad=0.)

    plt.savefig(output_path)

    return None


def kmeans_cluster(data, llambda_axis=-1):
    """kmeans clustering: uses a kmeans algorithm to cluster data,
       in order to independently cross correlate the skimage clustering method.
        ----------
       Parameters:
            data (numpy array): data to be clustered
            llambda_axis (int): index for wavelength, -1 if scalar array.
        ----------
        Returns:
            labels (numpy array): an array of labels, with 0 = granules,
                                  2 = intergranules, 1 = in-between.
            """

    if llambda_axis not in [-1, 2]:
        raise Exception('Wrong data shape. \
        (either scalar or (x, y, llambda) )')
    n_clusters = 3
    n_init = 20
    x_size = np.shape(data)[0]
    y_size = np.shape(data)[1]
    if llambda_axis == -1:
        data_flat = np.reshape(data, (x_size * y_size, 1))
        labels_flat = KMeans(n_clusters).fit(data_flat).labels_
        labels = np.reshape(labels_flat, (x_size, y_size))
    else:
        llambda_size = np.shape(data)[llambda_axis]
        data = np.reshape(data, (x_size * y_size, llambda_size))
        labels = np.reshape(Kmeans(n_clusters, n_init).fit(data),
                            (x_size, y_size))

    # now, make intergranules 0, granules 1:
    group0_mean = np.mean(data[labels == 0])
    group1_mean = np.mean(data[labels == 1])
    group2_mean = np.mean(data[labels == 2])
    # intergranules
    min_index = np.argmin([group0_mean,
                           group1_mean,
                           group2_mean])
    return_labels = np.ones(labels.shape)

    return_labels[[labels[:, :] == min_index][0]] -= 1

    return return_labels


def cross_correlation(segment1, segment2):
    """
    Return -1 and print a message if the agreement between two
    arrays is low, 0 otherwise. Designed to be used with segment.py
    and funclib.kmeans function.
    ----------
    Parameters:
        segment1 (numpy array): 'main' array to compare the other input
                                 array against.
        segment2 (numpy array): 'other' array (i.e. data
                                 segmented using kmeans).
    ----------
    Returns:
        [label, confidence] (list): where
        label is a label to summarize the confidence metric (int):
            -1: if agreement is low (below 75%)
            0: otherwise
        confidence is the actual confidence metric (float):
            float between 0 and 1 (0 if no agreement,
                                   1 if completely agrees)
    """

    total_granules = np.count_nonzero(segment1 == 1)
    total_intergranules = np.count_nonzero(segment1 == 0)

    if total_granules == 0:
        raise Exception('clustering problematic (no granules found)')

    if total_intergranules == 0:
        raise Exception('clustering problematic (no intergranules found)')

    x_size = np.shape(segment1)[0]
    y_size = np.shape(segment1)[1]

    granule_agreement_count = 0
    intergranule_agreement_count = 0
    for i in range(x_size):
        for j in range(y_size):
            if segment1[i, j] == 1 and segment2[i, j] == 1:
                granule_agreement_count += 1
            elif segment1[i, j] == 0 and segment2[i, j] == 0:
                intergranule_agreement_count += 1

    percentage_agreement_granules = \
        granule_agreement_count / total_granules
    percentage_agreement_intergranules = \
        intergranule_agreement_count / total_intergranules
    try:
        confidence = np.mean(percentage_agreement_granules,
                             percentage_agreement_intergranules)
    except TypeError:
        confidence = 0

    if percentage_agreement_granules < 0.75 \
            or percentage_agreement_intergranules < 0.75:
        print('Low agreement with K-Means clustering. \
                         Saved output has low confidence.')
        return [-1, confidence]
    else:
        return [0, confidence]
