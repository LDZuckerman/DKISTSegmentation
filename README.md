# DKISTSegmentation
Pipeline to segment DKIST images for granule detection

## Summary

This repository contains a pipeline for segmenting optical images of the solar
photosphere to identify granules. We segment into tri-valued images with 0 = 
inter-granule, 0.5 = faculae, and 1 = granule. 

These segmentation pipelines are currentlyy implemented for data from the IBIS
and DKIST telescopes.

For possible future integration with the SunPy package, we first convert test 
data to SunPy maps, and outputs are converted to SunPy maps as well. 

## Conda Environment

This repository will work correctly using the conda envrironment stored in the
environment.yml file. To create this environment from this file, use 
"conda env create -f environment.yml".

## Usage

The key functionality is contained in the `segment` function in funclib.py, which
takes in a SunPy map and a thresholding method and performs the segmentation, with 
the segmented image out as a SunPy map. Calls to segment look like:
```
segmented_data = segment(data, skimage_method, plot_intermed)
```
where 
* `data` is a SunPy map containing an optical image of the solar photosphere
* `skimage_method` is the base sci-kit image (skimage) segmentation technique 
   to use (see note on thresholds below)
* `plot_intermed` is a boleen flag indiacting whether the user would like to 
   save an image showing the intermediate data products (default True) If True, 
   will create a plot `intermediate_outputs.png` showing the steps of 
   segmentation processing. 

We recomend running this function through the `segment.py` script, as per the
example below. `segment.py` also includes an option to save the output as a 
fits file `output.fits` for easy inspection. In this output file, the first 
HDU is the segmented map and the second HDU is the input data. 

## Example

The script segment.py shows an example use. It can be run on either DKIST or IBIS data. Example data files are located in the `data` directory. 

To perform the segmentation on the IBIS data, creating the optional output plots 
and fits file:

```
python segment.py\
    --input_file 'data/IBIS.granulation.aligned.25Apr2019.seq56.sav'\
    --skimage_method 'li'\
    --plot_intermed True\
    --out_file 'output.fits\
    --vel_comparison_file 'velocity_comparison'\
    --out_dir 'output_IBIS/'
```

To perform the segmentation on the DKIST data, creating the optional output plots 
and fits file:

```
python segment.py\
    --input_file 'data/dkist.cont789nm.scaled.fits' \
    --skimage_method 'li' \
    --plot_intermed True \
    --out_file 'output.fits' \
    --out_dir 'output_DKIST/'
```

where the arugments are \
`input_file`: Input datafile path.\
`skimage_method`: Skimage method to use for initial thresholding\
`plot_intermed`: True/False - whether or not to save an intermediate data products image\
`out_file`: (Optional) Desired name of output fits file containing segmented map (extension 0) and input map (extension 1)\
`vel_comparison_file`: (Optional) Desired name of output image file containing segmented contours overplotted on velocity data.\
`out_dir`: (Optional) Desired directory in which to save out_file.

## Thresholding methods

The pipeline is built on existing sci-kit image (skimage) packages to determine
thresholding values. These are passed to the segment function as 'skimage-
method'. We recomend use of the 'Li' or 'mean' methods, although all should be
tested on a given dataset. We retain the functionality for the user to define the method used for generality to other datasets. We have determined that the 'Li' method performs best on both IBIS and DKIST data becuase the cross-correlation of velocity data with segementation labels is highest using this method.