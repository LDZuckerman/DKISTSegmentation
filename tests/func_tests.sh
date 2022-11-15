test -e ssshtest || wget -q https://raw.githubusercontent.com/ryanlayer/ssshtest/master/ssshtest
. ssshtest


# check that output plot directory is correct
run test_segment python segment.py  \
    --input_file 'data/IBIS.granulation.aligned.25Apr2019.seq56.sav' \
    --skimage_method 'li' \
    --plot_intermed True \
    --out_file 'output.fits' \
    --out_dir 'output/'

curr_dir=$(pwd)
file='/output/intermediate_outputs.png'
path=$curr_dir$file
# check that the output file actually exists
if [ -f "$path" ]; then
   echo ' TEST SUCCEEDED: output plot found in '$path 
else
   echo ' TEST FAILED: output plot not found in '$path 
fi
rm $path
assert_exit_code 0