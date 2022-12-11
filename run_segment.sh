#!/usr/bin/bash
set -e  # stop on error
set -u  # raise error if variable is unset
set -o pipefail  # fail if any prior step failed

# to run on IBIS data
ibis_res=0.096
python segment.py --data_path 'data/IBIS' --resolution $ibis_res --skimage_method 'li' --plot_intermed True --mark_dim_centers True --out_file 'segmented_data' --vel_comparison_file 'velocity_comparison' --out_dir 'example_outputs/IBIS/'

# to run on DKIST data
dkist_res=0.016
python segment.py --data_path 'data/DKIST' --resolution $dkist_res --skimage_method 'li' --plot_intermed True  --mark_dim_centers True --out_file 'segmented_data' --out_dir 'example_outputs/DKIST/'
