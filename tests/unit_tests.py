import sys
sys.path.append('..')  # nopep8
import unittest
import funclib
import random
import numpy as np
import sunpy.map as sm
import scipy
import sunpy
from astropy.io import fits
import os
import pathlib as pl


class TestUtils(unittest.TestCase):

    @classmethod
    def setUpClass(cls):

        """ Set up for unit testing by creating toy data
        """
        cls.testfile = 'data/IBIS.granulation.aligned.25Apr2019.seq56.sav'
        cls.dkist_testfile = 'data/dkist.cont789nm.scaled.fits'
        cls.test_instrument = 'IBIS'
        cls.dkist_instrument = 'DKIST'
        cls.test_method = 'li'
        cls.test_band = 'rosa_gband'

    @classmethod
    def tearDownClass(cls):

        """ Tear down unit testing toy data
        """

        cls.testfile = None
        cls.test_instrument = None
        cls.test_method = None
        cls.test_band = None

    def test_sav_to_map(cls):

        """ Unit tests for sav_to_map() function
        """

        # -------- positive tests -------- :
        # check that the returned type is as expected
        returned_map = funclib.sav_to_map(cls.testfile, cls.test_band)
        test_type = type(returned_map)
        cls.assertEqual(test_type, sunpy.map.mapbase.GenericMap)
        # check that returned data is not empty
        cls.assertTrue(np.size(returned_map.data) > 0)

        # -------- negative tests -------- :
        # check that no other type is returned
        cls.assertFalse(test_type != sunpy.map.mapbase.GenericMap)
        cls.assertNotEqual(test_type, np.ndarray)

        # ------ error raising tests ------ :
        # check that errors are raised for incorrect inputs
        cls.assertRaises(Exception, funclib.sav_to_map, 'ABC.txt',
                         cls.test_band)
        cls.assertRaises(Exception, funclib.sav_to_map, '', cls.test_band)

    def test_fits_to_map(cls):

        """ Unit tests for fits_to_map() function
        """

        # -------- positive tests -------- :
        # check that the returned type is as expected
        returned_map = funclib.fits_to_map(cls.dkist_testfile)
        test_type = type(returned_map)
        cls.assertEqual(test_type, sunpy.map.mapbase.GenericMap)
        # check that returned data is not empty
        cls.assertTrue(np.size(returned_map.data) > 0)

        # -------- negative tests -------- :
        # check that no other type is returned
        cls.assertFalse(test_type != sunpy.map.mapbase.GenericMap)
        cls.assertNotEqual(test_type, np.ndarray)

        # ------ error raising tests ------ :
        # check that errors are raised for incorrect inputs
        cls.assertRaises(Exception, funclib.fits_to_map, 'ABC.txt')
        cls.assertRaises(Exception, funclib.fits_to_map, '')

    def test_sav_to_numpy(cls):

        """ Unit tests for sav_to_numpy() function
        """

        # -------- positive tests -------- :
        returned_array = funclib.sav_to_numpy(cls.testfile,
                                              cls.test_instrument,
                                              cls.test_band)
        # check that correct type is returned
        test_type = type(returned_array)
        cls.assertEqual(test_type, np.ndarray)
        # check that returned data is not empty
        cls.assertTrue(np.size(returned_array) > 0)

        # -------- negative tests -------- :
        # check that no other type is returned
        cls.assertFalse(test_type != np.ndarray)
        cls.assertNotEqual(test_type, sunpy.map.mapbase.GenericMap)

        # ------ error raising tests ------ :
        # check that errors are raised for incorrect inputs
        cls.assertRaises(Exception, funclib.sav_to_numpy, 'ABC.txt',
                         cls.test_instrument, cls.test_band)
        cls.assertRaises(Exception, funclib.sav_to_numpy, cls.testfile,
                         'telescope', cls.test_band)
        cls.assertRaises(Exception, funclib.sav_to_numpy, cls.testfile,
                         cls.test_instrument, 'visible')

    def test_segment(cls):
        """ Unit tests for segment() function
        """

        # -------- positive tests -------- :
        # check that the returned type is correct
        data_map = funclib.sav_to_map(cls.testfile, cls.test_band)
        segmented = funclib.segment(data_map, cls.test_method)
        test_type = type(segmented)
        cls.assertEqual(test_type, sunpy.map.mapbase.GenericMap)
        # get an array of pixel values and check it is not empty
        initial_pix = sm.all_pixel_indices_from_map(data_map).value
        seg_pixels = sm.all_pixel_indices_from_map(segmented).value
        cls.assertTrue(np.size(seg_pixels) > 0)
        # check that the returned shape is unchanged
        cls.assertEqual(seg_pixels.shape, initial_pix.shape)

        # -------- negative tests -------- :
        # check that the values in the array are changed
        random.seed(42)
        # pick 10 random indices to check
        x = random.sample(list(np.arange(0, data_map.data.shape[0], 1)), 10)
        y = random.sample(list(np.arange(0, data_map.data.shape[1], 1)), 10)
        for i in range(len(x)):
            cls.assertNotEqual(data_map.data[x[i], y[i]],
                               segmented.data[x[i], y[i]])

        # ------ error raising tests ------ :
        cls.assertRaises(TypeError, funclib.segment, data_map, 'skimage')
        cls.assertRaises(TypeError, funclib.segment,
                         funclib.sav_to_numpy(cls.testfile,
                                              cls.test_instrument,
                                              cls.test_band), 'skimage')

    def test_get_threshold(cls):
        """ Unit tests for get_threshold() function
        """

        # -------- positive tests -------- :
        # check type of output
        test_arr1 = np.array([[1, 2, 3, 4, 5], [1, 2, 3, 4, 5]])
        threshold1 = funclib.get_threshold(test_arr1, cls.test_method)
        cls.assertEqual(type(threshold1), np.float64)
        # check range of output
        cls.assertTrue((threshold1 > 0) * (threshold1 < np.max(test_arr1)))

        # -------- negative tests -------- :
        # check that different arrays return different thresholds
        test_arr2 = np.array([[2, 3, 4, 5, 6], [2, 3, 4, 5, 6]])
        threshold2 = funclib.get_threshold(test_arr2, cls.test_method)
        cls.assertNotEqual(threshold1, threshold2)

        # ------ error raising tests ------ :
        cls.assertRaises(ValueError, funclib.get_threshold, test_arr1,
                         'banana')
        cls.assertRaises(ValueError, funclib.get_threshold, [],
                         cls.test_method)

    def test_remove_middles(cls):
        """ Unit tests for remove_middles() function
        """

        # -------- positive tests -------- :
        data_map = funclib.sav_to_map(cls.testfile, cls.test_band)
        thresholded = np.uint8(data_map.data > np.nanmedian(data_map.data))
        # check that returned array is not empty
        cls.assertTrue(np.size(thresholded) > 0)
        # check that the correct dimensions are returned
        cls.assertEqual(thresholded.shape,
                        funclib.remove_middles(thresholded).shape)

        # -------- negative tests -------- :
        # check that the returned array has fewer (or same number)
        # 0-valued pixels as input array (for a data set which we
        # know by eye should have some middle sections removed)
        middles_removed = funclib.remove_middles(thresholded)
        cls.assertFalse(np.count_nonzero(middles_removed)
                        < np.count_nonzero(thresholded))

        # ------ error raising tests ------ :
        # check that raises error if passed array is not binary
        cls.assertRaises(ValueError, funclib.remove_middles, data_map.data)

    def test_mark_faculae(cls):
        """ Unit tests for mark_faculae() function
        """

        # -------- positive tests -------- :
        data_map = funclib.sav_to_map(cls.testfile, cls.test_band)
        thresholded = np.uint8(data_map.data > np.nanmedian(data_map.data))
        faculae_marked = funclib.mark_faculae(thresholded, data_map.data)
        # check that the correct dimensions are returned
        cls.assertEqual(thresholded.shape,
                        faculae_marked.shape)
        # check that returned array is not empty
        cls.assertTrue(np.size(thresholded) > 0)

        # -------- negative tests -------- :
        # check that the returned array has some 0.5 values (for a dataset
        # that we know has faculae by eye)
        cls.assertNotEqual(len(np.where(faculae_marked == 0.5)[0]), 0)

        # ------ error raising tests ------ :
        cls.assertRaises(ValueError, funclib.mark_faculae, data_map.data,
                         data_map.data)

    def test_save_to_fits(cls):
        """ Unit tests for save_to_fits() function
        """

        data_map = funclib.sav_to_map(cls.testfile, cls.test_band)
        segmented_map = funclib.segment(data_map, cls.test_method)
        funclib.save_to_fits(segmented_map, data_map,
                             'test_output.fits',
                             'output/')
        path = pl.Path('output/test_output.fits')
        read_seg_data = fits.open('output/test_output.fits')[0].data
        read_data = fits.open('output/test_output.fits')[1].data

        # positive tests
        # check that fits file gets created and two data is correct
        cls.assertTrue(os.path.exists(path))
        cls.assertTrue(np.array_equal(read_seg_data, segmented_map.data))
        cls.assertTrue(np.array_equal(read_data, data_map.data))

        # negtaive tests
        # check that that data maps are not switched
        cls.assertFalse(np.array_equal(read_seg_data, data_map.data))
        cls.assertFalse(np.array_equal(read_data, segmented_map.data))

        # error raising tests
        cls.assertRaises(TypeError,
                         funclib.save_to_fits,
                         data_map.data, segmented_map,
                         'test_output.fits',
                         'output/')
        cls.assertRaises(TypeError,
                         funclib.save_to_fits,
                         data_map,
                         segmented_map,
                         'test_output.fits',
                         4)

        os.remove('output/test_output.fits')

    def test_find_files(self):
        self.assertTrue(True)
        self.assertFalse(False)
