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

    args = parser.parse_args()
    skimage_method = args.skimage_method
    input_file = args.input_file
    plot_intermed = args.plot_intermed
  
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

    # save map as fits file
    if args.out_file is not None:
        funclib.save_to_fits(segmented_map,
                             data_map,
                             args.out_file,
                             args.out_dir)


if __name__ == "__main__":
    main()
