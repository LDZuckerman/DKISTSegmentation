import argparse
import funclib


def main():

    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--input_file', dest='input_file',
                        type=str,
                        help='Input datafile path.',
                        required=True)
    parser.add_argument('--skimage_method', dest='skimage_method',
                        type=str,
                        help='Skimage method to use for initial thresholding',
                        required=True)
    parser.add_argument('--plot_intermed',
                        type=str,
                        help='True/False - Whether or not to save an '
                             + 'intermediate data products image',
                        required=True)
    parser.add_argument('--out_file',
                        type=str,
                        help='(Optional) Desired name of output fits file '
                             + 'containing segmented map (extension 0) and '
                             + 'input map (extension 1)',
                        default=False,
                        required=False)
    parser.add_argument('--out_dir',
                        type=str,
                        help='(Optional) Desired directory in which to save '
                             + 'out_file',
                        default='output/',
                        required=False)
    parser.add_argument('--vel_comparison_file',
                        type=str,
                        help='(Optional) Desired name of output image file '
                             + 'containing segmented contours overplotted '
                             + 'on velocity data.',
                        default=None,
                        required=False)

    args = parser.parse_args()
    skimage_method = args.skimage_method
    input_file = args.input_file
    out_dir = args.out_dir
    plot_intermed = args.plot_intermed

    # read data into map to mimic use within SunPy
    if input_file.endswith('.sav'):
        data_map = funclib.sav_to_map(input_file, 'rosa_gband')
    if input_file.endswith('.fits'):
        data_map = funclib.fits_to_map(input_file)

    # check for res flags
    if "dkist" in input_file.lower():
        res = "DKIST"
    elif "ibis" in input_file.lower():
        res = "IBIS"
    else:
        print('Currently file name must contain "dkist" or "ibis", for ' +
              'proper faculae detection. This will be updated in a later ' +
              'verions. Defaulting to DKIST.')
        res = "DKIST"

    # apply segmentation pipeline
    segmented_map = funclib.segment(data_map,
                                    skimage_method,
                                    plot_intermed,
                                    out_dir,
                                    res)

    # create a visual comparison against velocity data
    if args.vel_comparison_file is not None:
        funclib.overplot_velocities(segmented_map,
                                    input_file,
                                    out_dir + '/' + args.vel_comparison_file)

    # save map as fits file
    if args.out_file is not None:
        funclib.save_to_fits(segmented_map,
                             data_map,
                             args.out_file,
                             out_dir)

        # check out put via kmeans:
        # still working on IBIS data:
        if input_file.endswith('.fits'):
            kmeans_labels = funclib.kmeans_cluster(data_map, llambda_axis=-1)
            funclib.cross_correlation(segmented_map, kmeans_labels)


if __name__ == "__main__":
    main()
