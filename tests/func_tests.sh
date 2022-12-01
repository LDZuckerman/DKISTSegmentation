test -e ssshtest || wget -q https://raw.githubusercontent.com/ryanlayer/ssshtest/master/ssshtest
. ssshtest

   
# --------------- Checks for running on IBIS data ------------------ #
# check run on IBIS data, with velocity comparison 
test_dir='/output_IBIS_test'
run test_segment_ibis python segment.py \
    --data_path 'data/IBIS' \
    --skimage_method 'li' \
    --plot_intermed True \
    --out_file 'output' \
    --vel_comparison_file 'velocity_comparison' \
    --out_dir 'output_IBIS_test/'

curr_dir=$(pwd)
file_id='IBIS_example'
file='/intermediate_outputs_'$file_id'.png'
path=$curr_dir$test_dir$file

# check that the output file actually exists
if [ -f "$path" ]; then
   echo ' TEST SUCCEEDED: intermediate output plot found in '$path 
else
   echo ' TEST FAILED: intermediate output plot not found in '$path 
fi
rm $path
assert_exit_code 0

file='/output_'$file_id'.fits'
path=$curr_dir$test_dir$file
# check that the output file actually exists
if [ -f "$path" ]; then
   echo ' TEST SUCCEEDED: output datafile found in '$path 
else
   echo ' TEST FAILED: output datafile not found in '$path 
fi
rm $path
assert_exit_code 0

file='/velocity_comparison_'$file_id'.png'
path=$curr_dir$test_dir$file
# check that the output file actually exists
if [ -f "$path" ]; then
   echo ' TEST SUCCEEDED: velocity comparison plot found in '$path 
else
   echo ' TEST FAILED: velocity comparison plot not found in '$path 
fi
rm $path
assert_exit_code 0

rm -r $curr_dir$test_dir


# --------------- Checks for running on DKIST data ------------------ #
# check run on IBIS data, with velocity comparison 
#test_dir='/output_DKIST_test'
#run test_segment_ibis python segment.py \
#    --data_path 'data/DKIST' \
#    --skimage_method 'li' \
#    --plot_intermed True \
#    --out_file 'output' \
#    --out_dir $test_dir

#curr_dir=$(pwd)
#file_id='DKIST_example'
#file='/intermediate_outputs_'$file_id'.png'
#path=$curr_dir$test_dir$file
## check that the output file actually exists
#if [ -f "$path" ]; then
#   echo ' TEST SUCCEEDED: intermediate output plot found in '$path 
#else
#   echo ' TEST FAILED: intermediate output plot not found in '$path 
#fi
#rm $path
#assert_exit_code 0

#file='/output_'$file_id'.fits'
#path=$curr_dir$test_dir$file
## check that the output file actually exists
#if [ -f "$path" ]; then
#   echo ' TEST SUCCEEDED: output datafile found in '$path 
#else
#   echo ' TEST FAILED: output datafile not found in '$path 
#fi
#rm $path
#assert_exit_code 0

#rm -r $curr_dir$test_dir