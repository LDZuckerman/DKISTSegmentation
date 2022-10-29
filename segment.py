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

    args = parser.parse_args()
    skimage_method = args.skimage_method
    input_file = args.input_file
        
    # read data into map to mimic use within SunPy
    data_map = funclib.sav_to_map(input_file, 'rosa_gband')

    # apply segmentation pipeline
    segmented_map = funclib.segment(data_map, skimage_method)
    
    
if __name__ == "__main__":
    main()
