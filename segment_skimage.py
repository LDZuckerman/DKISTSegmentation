import sys
import argparse
import scipy
import matplotlib.pyplot as plt
import numpy as np
from skimage import data, img_as_float
import scipy.misc
import scipy.ndimage
import skimage.filters
import sklearn.metrics
import funclib

def main():
    """ Use skimage to segement solar data.

    Parameters
    ----------
    No input parameters; paramters are passed from
    command line call using argparse.

    Returns
    -------
    No returned output.

    Output
    ------
    Saves output plot.

    """

    parser = argparse.ArgumentParser(description='') 
    parser.add_argument('--input_file', dest='input_file',
                        type=str,
                        help='Input datafile path.',  # 'ie: ./IBIS.granulation.aligned.25Apr2019.seq56.sav'
                        required=True)    
    parser.add_argument('--method', dest='method',
                        type=str,
                        help='Method to use for segmentation',
                        required=True)
    parser.add_argument('--output_file', dest='output_file',
                        type=str,
                        help='Output file name to save plot.',
                        required=True)

    args = parser.parse_args()
    method = str(args.method)
    output_file = str(args.output_file)
    input_file = str(args.input_file)

    methods = ['otsu', 'li', 'isodata', 'mean', 'minimum', 'yen', 'triangle']
    if method not in methods:
        print('Method must be one of: ' + str(methods))
        sys.exit(1)
        
    # read in data
    data = funclib.sav_to_numpy(input_file, 'IBIS', 'rosa_gband')
    
    # apply median filter
    median_filtered = scipy.ndimage.median_filter(grayscale, size=3)

    # apply threshold
    if method == 'otsu':
        threshold == skimage.filters.threshold_otsu(median_filtered)
    elif method == 'li':
        threshold == skimage.filters.threshold_li(median_filtered)
    elif method == 'yen':
        threshold == skimage.filters.threshold_yen(median_filtered)
    elif method == 'mean':
        threshold == skimage.filters.threshold_mean(median_filtered)
    elif method == 'minimum':
        threshold = skimage.filters.threshold_minimum(median_filtered)
    elif method == 'triangle':
        threshold = skimage.filters.threshold_triangle(median_filtered)
    elif method == 'isodata':
        threshold == skimage.filters.threshold_isodata(median_filtered)
        
    plt.figure(figsize = [20,40])
    segmented_image = np.uint8(median_filtered > threshold) * 255
    plt.imshow(segmented, cmap='gray')
    plt.axis('off')
    plt.title('Image segmented using ' + method + 'method')
    plt.savefig(output_file)
    
    
if __name__ == "__main__":
    main()
