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
    """ Segement solar data.

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
    parser.add_argument('--skimage_method', dest='skimage_method',
                        type=str,
                        help='Skimage method to use for initial thresholding',
                        required=True)
    parser.add_argument('--output_file', dest='output_file',
                        type=str,
                        help='Output file name to save plot.',
                        required=True)

    args = parser.parse_args()
    skimage_method = str(args.skimage_method)
    output_file = str(args.output_file)
    input_file = str(args.input_file)
        
    # read in data
    data = funclib.sav_to_numpy(input_file, 'IBIS', 'rosa_gband')
    
    # apply median filter
    median_filtered = scipy.ndimage.median_filter(data, size=3)

    # apply threshold
    threshold = funclib.get_threshold(median_filtered, skimage_method)

    # initial skimage segmentation
    segmented_image = np.uint8(median_filtered > threshold) * 255  

    # fix the extra IGM bits in the middle of granules
    segmented_image_fixed = funclib.remove_middles(segmented_image)

    fig, axs = plt.subplots(3, 1, figsize=(15, 30))
    s1 = 20; s2 =26
    fig.suptitle(skimage_method + ' segmentation plus correction', fontsize=s2)    
    im0 = axs[0].imshow(data)
    axs[0].set_title('input image', fontsize=s1)
    plt.colorbar(im0, ax=axs[0])
    im1 = axs[1].imshow(segmented_image, cmap='gray')
    axs[1].set_title('direct skimage segmentation', fontsize=s1)
    plt.colorbar(im1, ax=axs[1])
    im2 = axs[2].imshow(segmented_image_fixed, cmap='gray')
    axs[2].set_title('correction segmented image', fontsize=s1)
    plt.colorbar(im2, ax=axs[2])
    plt.axis('off')
    plt.savefig(output_file)
    
    
if __name__ == "__main__":
    main()
