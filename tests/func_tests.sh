test -e ssshtest || wget -q https://raw.githubusercontent.com/ryanlayer/ssshtest/master/ssshtest
. ssshtest


# check that output plot directory is correct
run test_segment python segment.py  \
    --input_file 'IBIS.granulation.aligned.25Apr2019.seq56.sav' \
    --skimage_method 'li' \

curr_dir=$(pwd)
file='/intermediate_outputs.png'
path=$curr_dir$file
# check that the output file actually exists
if [ -f "$path" ]; then
   echo ' TEST SUCCEEDED: output plot found in expected location.'
else
   echo ' TEST FAILED: output plot not found in expected location!'
fi
rm $path
assert_exit_code 0