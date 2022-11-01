# DKISTSegmentation
Pipeline to segment DKIST images for granule detection

## Summary

This repository contains a pipeline for segmenting optical images of the solar
photosphere to identify granules. We segment into tri-valued images with 0 = 
inter-granule, 0.5 = faculae, and 1 = granule. 

Ideally, this technique will be accurate enough to apply to high-resolution 
DKIST data. Currently, we only have access to test data from IBIS.

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
segmented_data = segment(data, skimage_method)
```
where 
* `data` is a SunPy map containing an optical image of the solar photosphere
* `skimage_method` is the base sci-kit image (skimage) segmentation technique 
   to use (see note on thresholds below)

This will also create a plot `intermediate_outputs.png` showing the steps of
segmentation processing. 

## Thresholding methods

The pipeline is built on existing sci-kit image (skimage) packages to determine
thresholding values. These are passed to the segment function as 'skimage-
method'. We recomend use of the 'Li' or 'mean' methods, although all should be
tested on a given dataset.  

## Example

The script segment.py shows an example use on IBIS data. This IBIS test data
file IBIS.granulation.aligned.25Apr2019.seq56.sav test can be downloaded from
[Leah Zuckerman's google drive](https://drive.google.com/file/d/1bEYOzn3KXJp3EDTsF5Rgya3x8ABxXP2D/view?usp=share_link)

To perform the segmentation on this data:

python segment.py\ 
        --input_file IBIS.granulation.aligned.25Apr2019.seq56.sav 
        --skimage_method li 
        --output_file test_pipeline.png
