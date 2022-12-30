import argparse
import astropy.io.fits as fits
import funclib


def main():

    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--data_path',
                        dest='data_path',
                        type=str,
                        help='filepath to directory containing data to ' +
                        'be segmented.',
                        required=True)
    parser.add_argument('--resolution',
                        dest='resolution',
                        type=float,
                        help='resolution (arcseconds/pixel) of data in ' +
                        'input directory.',
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
    parser.add_argument('--mark_dim_centers',
                        type=str,
                        help='True/False - Whether or not to mark dim '
                             + 'intergranule centers as a seperate class',
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
    data_path = args.data_path
    res = float(args.resolution)
    skimage_method = args.skimage_method
    out_dir = args.out_dir
    out_file = args.out_file
    plot_intermed = args.plot_intermed
    mark_dim_centers = args.mark_dim_centers
    vel_comp_file = args.vel_comparison_file

    # define wavelength band of interest for dkist
    dkist_band = 'rosa_gband'

    # get the list of data:
    data_to_be_segmented = funclib.find_data(data_path)

    for input_file in data_to_be_segmented:

        # get the name of the input file (minus extenstion) for bookkeeping
        file_id = input_file.rsplit('.', 1)[0]

        # read data into map to mimic use within SunPy
        if input_file.endswith('.sav'):
            data_map = funclib.sav_to_map(data_path + '/' + input_file,
                                          dkist_band)
            header = None
        elif input_file.endswith('.fits'):
            data_map = funclib.fits_to_map(data_path + '/' + input_file)
            header = fits.open(data_path + '/' + input_file)[0].header

        # apply segmentation pipeline
        segmented_map = funclib.segment(file_id,
                                        data_map,
                                        skimage_method,
                                        res,
                                        plot_intermed,
                                        mark_dim_centers,
                                        out_dir)

        # check out put via kmeans:
        # still working on IBIS data:
        if input_file.endswith('.fits'):
            kmeans_labels = funclib.kmeans_cluster(data_map.data,
                                                   llambda_axis=-1)
            kmeans_return = funclib.cross_correlation(segmented_map.data,
                                                      kmeans_labels)
            confidence = kmeans_return[1]
        else:
            confidence = 0

        # save map as fits file
        if args.out_file is not None:
            funclib.save_to_fits(segmented_map,
                                 data_map,
                                 out_file + '_' + file_id + '.fits',
                                 out_dir,
                                 header,
                                 confidence)


if __name__ == "__main__":
    main()
