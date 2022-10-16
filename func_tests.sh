test -e ssshtest || wget -q https://raw.githubusercontent.com/ryanlayer/ssshtest/master/ssshtest
. ssshtest


# check that output plot directory is correct
#wget "https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"
run test_segment_skimage python segment_skimage.py  \
    --input_file 'IBIS.granulation.aligned.25Apr2019.seq56.sav' \
    --method 'li' \
    --output_file 'li_segmented.png' \

curr_dir=$(pwd)
file='/li_segmented.png'
path=$curr_dir$file
# check that the output file actually exists
if [ -f "$path" ]; then
   echo ' TEST SUCCEEDED: output plot found in expected location.'
else
   echo ' TEST FAILED: output plot not found in expected location!'
fi
assert_exit_code 0
