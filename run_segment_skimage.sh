#!/usr/bin/bash

set -e  # stop on error
set -u  # raise error if variable is unset
set -o pipefail  # fail if any prior step failed

python segment_skimage.py  \
    --input_file 'IBIS.granulation.aligned.25Apr2019.seq56.sav' \
    --method 'li' \
    --output_file 'li_segmented.png' 