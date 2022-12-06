# DKISTSegmentation
Pipeline to segment DKIST images for granule detection

## Summary

This repository contains a pipeline for segmenting optical images of the solar
photosphere to identify granules. We segment into tri-valued images with 0 = 
intergranule, 0.5 = faculae, and 1 = granule. 

These segmentation pipelines are currently implemented for data from the IBIS
and DKIST telescopes.

For possible future integration with the SunPy package, we first convert test 
data to SunPy maps, and outputs are converted to SunPy maps as well. 

## Conda Environment

This repository will work correctly using the conda envrironment stored in the
environment.yml file. To create this environment from this file, use 
"conda env create -f environment.yml".

## Input data

This pipeline is currently implemented for data from the IBIS and DKIST telescopes,
though other data in .sav or .fits format should also produce good results. Input 
data must be enclosed in a directory containing only data from one telescope.
If you wish perform segmentation on data from several telescopes, please
create unique directories for each and run on each directory individually.

## Usage

The key functionality is contained in the `segment` function in funclib.py, which
takes in a SunPy map and a thresholding method and performs the segmentation, with 
the segmented image out as a SunPy map. Calls to segment look like:
```
segmented_data = segment(file_id, data, skimage_method, plot_intermed)
```
where 
* `file_id` is a string which identifies the input file for bookkeping purposes
* `data` is a SunPy map containing an optical image of the solar photosphere
* `skimage_method` is the base sci-kit image (skimage) segmentation technique 
   to use (see note on thresholds below)
* `plot_intermed` is a boleen flag indiacting whether the user would like to 
   save an image showing the intermediate data products (default True). If True, 
   will create a plot  showing the steps of segmentation processing and final
   results. 

We recomend running this function through the `segment.py` script, as per the
example below. `segment.py` also includes an option to save the output as a 
fits file `output.fits` for easy inspection. In this output file, the first 
HDU is the segmented map and the second HDU is the input data. 

## Example

The script segment.py shows an example use. It can be run on either DKIST or IBIS data. Example data files are located in the `data` directory. 

To perform the segmentation on the IBIS data, creating the optional output plots 
and fits file:

```
python segment.py
    --data_path 'data/DKIST'
    --skimage_method 'li'
    --plot_intermed True
    --out_file 'segmented_data.fits
    --vel_comparison_file 'velocity_comparison'
    --out_dir 'output_IBIS/'
```

To perform the segmentation on the DKIST data, creating the optional output plots 
and fits file:

```
python segment.py
    --data_path 'data/DKIST' 
    --skimage_method 'li' 
    --plot_intermed True 
    --out_file 'segmented_data.fits' 
    --out_dir 'output_DKIST/'
```
The outputs of this example call are located in the `example_outputs\IBIS` directory.

where the arugments are \
`data_path`: filepath to directory containing data to be segmented\
`skimage_method`: Skimage method to use for initial thresholding\
`plot_intermed`: True/False - whether or not to save an intermediate data products image\
`out_file`: (Optional) Desired name of output fits file containing segmented map (extension 0) and input map (extension 1)\
`vel_comparison_file`: (Optional) Desired name of output image file containing segmented contours overplotted on velocity data.\
`out_dir`: (Optional) Desired directory in which to save out_file.

The outputs of this example call are located in the `example_outputs\DKIST` directory.

## Warning messages

When running the above on IBIS data, you may see several warning messages about 'dubious year' from the sunpy ERFA function. 
*These can be safely ignored.* Sunpy absolutely requires that a header be present in all sunpy.Map objects. Becuase .sav files do
not contain header information, when the input file is in the .sav format we create as generic a header as possible given the strict
requirements imposed by sunpy on header field types and formats. No empty or None fields are permitted by sunpy.Map creation function. 
The warning messages simply note that the placeholder date in the header is not a true data value.

You may also see a UserWarning from matplotlib stating that no contour levels were found within the data range. *This can also be
safely ignored.* 


## Thresholding methods

The pipeline is built on existing sci-kit image (skimage) packages to determine
thresholding values. These are passed to the segment function as 'skimage-
method'. We recomend use of the 'Li' or 'mean' methods, although all should be
tested on a given dataset. We retain the functionality for the user to define the method used for generality to other datasets. We have determined that the 'Li' method performs best on both IBIS and DKIST data becuase the cross-correlation of velocity data with segementation labels is highest using this method.