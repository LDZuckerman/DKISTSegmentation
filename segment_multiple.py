import argparse
import funclib


def main():

    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--filepath_to_data',
                        dest='filepath_to_data',
                        type=str,
                        help='filepath to directory containing data to ' +
                        'be segmented.',
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
    parser.add_argument('--vel_comparison_file',
                        type=str,
                        help='(Optional) Desired name of output image file '
                             + 'containing segmented contours overplotted '
                             + 'on velocity data.',
                        default=None,
                        required=False)
    parser.add_argument('--out_dir',
                        type=str,
                        help='(Optional) Desired directory in which to save '
                             + 'out_file',
                        default='output/',
                        required=False)

    args = parser.parse_args()
    skimage_method = args.skimage_method
    filepath_to_data = args.filepath_to_data
    plot_intermed = args.plot_intermed
    vel_comparison_file = args.vel_comparison_file
    out_dir = args.out_dir

    # get the list of data:
    data_to_be_segmented = funclib.find_data(filepath_to_data)
    for data in data_to_be_segmented:
        input_file = data

        # read data into map to mimic use within SunPy
        if input_file.endswith('.sav'):
            data_map = funclib.sav_to_map(input_file, 'rosa_gband')
        if input_file.endswith('.fits'):
            data_map = funclib.fits_to_map(input_file)

        # apply segmentation pipeline
        segmented_map = funclib.segment(data_map,
                                        skimage_method,
                                        plot_intermed,
                                        args.out_dir)

        # create a visual comparison against velocity data
        if vel_comparison_file is not None:
            funclib.overplot_velocities(segmented_map,
                                        input_file,
                                        out_dir + '/' + vel_comparison_file)

        # save map as fits file
        if args.out_file is not None:
            funclib.save_to_fits(segmented_map,
                                 data_map,
                                 args.out_file,
                                 args.out_dir)


if __name__ == "__main__":
    main()