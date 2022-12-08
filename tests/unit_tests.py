import sys

sys.path.append('..')  # nopep8
import unittest
import funclib
import random
import numpy as np
import sunpy.map as sm
import sunpy
from astropy.io import fits
import os
import pathlib as pl
import scipy
import shutil


class TestUtils(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        """ Set up for unit testing by creating toy data
        """
        cls.ibis_testfile = 'data/IBIS/IBIS_example.sav'
        cls.ibis_res = 0.096
        cls.dkist_testfile = 'data/DKIST/DKIST_example.fits'
        cls.dkist_res = 0.016
        cls.dkist_fileid = 'DKIST_example'
        cls.ibis_fileid = 'IBIS_example'
        cls.ibis_instrument = 'IBIS'
        cls.dkist_instrument = 'DKIST'
        cls.test_method = 'li'
        cls.test_band = 'rosa_gband'
        cls.test_header = None

    @classmethod
    def tearDownClass(cls):
        """ Tear down unit testing toy data
        """

        cls.ibis_testfile = None
        cls.dkist_testfile = None
        cls.ibis_fileid = None
        cls.dkist_fileid = None
        cls.ibis_instrument = None
        cls.dkist_instrument = None
        cls.test_method = None
        cls.test_band = None
        shutil.rmtree('test_output')

    def test_save_to_fits(self):
        """ Unit tests for save_to_fits() function
        """

        data_map = funclib.sav_to_map(self.ibis_testfile, self.test_band)
        segmented_map = funclib.segment(self.ibis_fileid,
                                        data_map,
                                        self.test_method,
                                        True,
                                        'test_output/',
                                        self.ibis_res)
        confidence = 0
        funclib.save_to_fits(segmented_map, data_map, 'test_output.fits',
                             'output/', self.test_header, confidence)
        path = pl.Path('output/test_output.fits')
        read_seg_data = fits.open('output/test_output.fits')[0].data
        read_data = fits.open('output/test_output.fits')[1].data

        # -------- positive tests -------- :
        # Test 1: check that fits file gets created and two data fields are
        # correct
        self.assertTrue(os.path.exists(path))
        self.assertTrue(np.array_equal(read_seg_data, segmented_map.data))
        self.assertTrue(np.array_equal(read_data, data_map.data))

        # -------- negative tests -------- :
        # Test 2: check that that data maps are not switched
        self.assertFalse(np.array_equal(read_seg_data, data_map.data))
        self.assertFalse(np.array_equal(read_data, segmented_map.data))

        # ------ error raising tests ------ :
        # Test 4: check that errors are raised for incorrect inputs
        self.assertRaises(TypeError,
                          funclib.save_to_fits, data_map.data,
                          segmented_map, 'test_output.fits',
                          'output/', self.test_header, confidence)
        self.assertRaises(TypeError,
                          funclib.save_to_fits, data_map,
                          segmented_map, 'test_output.fits',
                          4, self.test_header, confidence)

        os.remove('output/test_output.fits')

    def test_sav_to_map(self):
        """ Unit tests for sav_to_map() function
        """

        # -------- positive tests -------- :
        # Test 1: check that the returned type is as expected
        returned_map = funclib.sav_to_map(self.ibis_testfile, self.test_band)
        test_type = type(returned_map)
        self.assertEqual(test_type, sunpy.map.mapbase.GenericMap)
        # Test 2: check that returned data is not empty
        self.assertTrue(np.size(returned_map.data) > 0)

        # -------- negative tests -------- :
        # Test 3: check that input type is not returned
        self.assertNotEqual(test_type, np.ndarray)
        # Test 4: check that no other type is returned
        self.assertFalse(test_type != sunpy.map.mapbase.GenericMap)

        # ------ error raising tests ------ :
        # Test 5: check that errors are raised for incorrect inputs
        self.assertRaises(Exception, funclib.sav_to_map, 'ABC.txt',
                          self.test_band)
        self.assertRaises(Exception, funclib.sav_to_map, '', self.test_band)

    def test_fits_to_map(self):
        """ Unit tests for fits_to_map() function
        """

        # -------- positive tests -------- :
        # Test 1: check that the returned type is as expected
        returned_map = funclib.fits_to_map(self.dkist_testfile)
        test_type = type(returned_map)
        self.assertEqual(test_type, sunpy.map.mapbase.GenericMap)
        # Test 2: check that returned data is not empty
        self.assertTrue(np.size(returned_map.data) > 0)

        # -------- negative tests -------- :
        # Test 3: check that input type is not returned
        self.assertFalse(test_type != sunpy.map.mapbase.GenericMap)
        # Test 4: check that no other type is returned
        self.assertNotEqual(test_type, np.ndarray)

        # ------ error raising tests ------ :
        # Test 5: check that errors are raised for incorrect inputs
        self.assertRaises(Exception, funclib.fits_to_map, 'ABC.txt')
        self.assertRaises(Exception, funclib.fits_to_map, '')

    def test_sav_to_numpy(self):
        """ Unit tests for sav_to_numpy() function
        """

        # -------- positive tests -------- :
        returned_array = funclib.sav_to_numpy(self.ibis_testfile,
                                              self.ibis_instrument,
                                              self.test_band)
        # Test 1: check that correct type is returned
        test_type = type(returned_array)
        self.assertEqual(test_type, np.ndarray)
        # Test 2: check that returned data is not empty
        self.assertTrue(np.size(returned_array) > 0)

        # -------- negative tests -------- :
        # Test 3: check that common output type is not returned
        self.assertFalse(test_type != np.ndarray)
        # Test 4: check that no other type is returned
        self.assertNotEqual(test_type, sunpy.map.mapbase.GenericMap)

        # ------ error raising tests ------ :
        # Test 5: check that errors are raised for incorrect inputs
        self.assertRaises(Exception,
                          funclib.sav_to_numpy,
                          'ABC.txt',
                          self.ibis_instrument,
                          self.test_band)
        self.assertRaises(Exception,
                          funclib.sav_to_numpy,
                          self.ibis_testfile,
                          'telescope',
                          self.test_band)
        self.assertRaises(Exception,
                          funclib.sav_to_numpy,
                          self.ibis_testfile,
                          self.ibis_instrument,
                          'visible')

    def test_segment(self):
        """ Unit tests for segment() function
        """

        # -------- positive tests -------- :
        data_map = funclib.sav_to_map(self.ibis_testfile, self.test_band)
        segmented = funclib.segment(self.ibis_fileid, data_map,
                                    self.test_method, True, 'test_output/',
                                    self.ibis_res)
        test_type = type(segmented)
        # Test 1: check that the returned type is correct
        self.assertEqual(test_type, sunpy.map.mapbase.GenericMap)
        # Test 2: get an array of pixel values and check it is not empty
        initial_pix = sm.all_pixel_indices_from_map(data_map).value
        seg_pixels = sm.all_pixel_indices_from_map(segmented).value
        self.assertTrue(np.size(seg_pixels) > 0)
        # Test 3: check that the returned shape is unchanged
        self.assertEqual(seg_pixels.shape, initial_pix.shape)

        # -------- negative tests -------- :
        # Test 4: check that the values in the array are changed
        random.seed(42)
        # pick 10 random indices to check
        x = random.sample(list(np.arange(0, data_map.data.shape[0], 1)), 10)
        y = random.sample(list(np.arange(0, data_map.data.shape[1], 1)), 10)
        for i in range(len(x)):
            self.assertNotEqual(data_map.data[x[i], y[i]],
                                segmented.data[x[i], y[i]])

        # ------ error raising tests ------ :
        # Test 5: check that errors are raised for incorrect inputs
        self.assertRaises(TypeError, funclib.segment, self.ibis_testfile,
                          data_map, 'skimage')
        self.assertRaises(TypeError, funclib.segment,
                          funclib.sav_to_numpy(self.ibis_testfile,
                                               self.ibis_instrument,
                                               self.test_band), 'skimage')

    def test_get_threshold(self):
        """ Unit tests for get_threshold() function
        """

        # -------- positive tests -------- :
        # Test 1: check type of output
        test_arr1 = np.array([[1, 2, 3, 4, 5], [1, 2, 3, 4, 5]])
        threshold1 = funclib.get_threshold(test_arr1, self.test_method)
        self.assertEqual(type(threshold1), np.float64)
        # Test 2: check range of output
        self.assertTrue((threshold1 > 0) * (threshold1 < np.max(test_arr1)))

        # -------- negative tests -------- :
        # Test 3: check that different arrays return different thresholds
        test_arr2 = np.array([[2, 3, 4, 5, 6], [2, 3, 4, 5, 6]])
        threshold2 = funclib.get_threshold(test_arr2, self.test_method)
        self.assertNotEqual(threshold1, threshold2)

        # ------ error raising tests ------ :
        # Test 4: check that errors are raised for incorrect inputs
        self.assertRaises(ValueError,
                          funclib.get_threshold,
                          test_arr1,
                          'banana')
        self.assertRaises(ValueError,
                          funclib.get_threshold,
                          [],
                          self.test_method)

    def test_trim_interganules(self):
        """ Unit tests for trim_intergranules() function
        """

        # -------- positive tests -------- :
        data_map = funclib.sav_to_map(self.ibis_testfile, self.test_band)
        thresholded = np.uint8(data_map.data > np.nanmedian(data_map.data))
        # Test 1: check that returned array is not empty
        self.assertTrue(np.size(thresholded) > 0)
        # Test 2: check that the correct dimensions are returned
        self.assertEqual(thresholded.shape,
                         funclib.trim_interganules(thresholded).shape)
        # Test 3: mark erronous material, not remove when flag is True.
        middles_marked = \
            funclib.trim_interganules(thresholded, mark=True)
        marked_erroneous = \
            np.count_nonzero(middles_marked[middles_marked == 2])
        self.assertNotEqual(marked_erroneous, 0)

        # -------- negative tests -------- :
        # Test 4: check that the returned array has fewer (or same number)
        # 0-valued pixels as input array (for a data set which we
        # know by eye should have some middle sections removed)
        middles_removed = funclib.trim_interganules(thresholded)
        self.assertFalse(np.count_nonzero(middles_removed)
                         < np.count_nonzero(thresholded))

        # ------ error raising tests ------ :
        # Test 5: check that raises error if passed array is not binary
        self.assertRaises(ValueError, funclib.trim_interganules, data_map.data)

    def test_mark_faculae(self):
        """ Unit tests for mark_faculae() function
        """

        # -------- positive tests -------- :
        data_map = funclib.sav_to_map(self.ibis_testfile, self.test_band)
        thresholded = np.uint8(data_map.data > np.nanmedian(data_map.data))
        faculae_marked = funclib.mark_faculae(thresholded,
                                              data_map.data,
                                              res=self.ibis_res)

        # Test 1: check that the correct dimensions are returned
        self.assertEqual(thresholded.shape,
                         faculae_marked.shape)
        # Test 2: check that returned array is not empty
        self.assertTrue(np.size(faculae_marked) > 0)

        # -------- negative tests -------- :
        # Test 3: check that the returned array has some 0.5 values (for a
        # dataset that we know has faculae by eye)
        self.assertNotEqual(len(np.where(faculae_marked == 0.5)[0]), 0)

        # ------ error raising tests ------ :
        # Test 4: check that errors are raised for incorrect inputs
        self.assertRaises(ValueError,
                          funclib.mark_faculae,
                          data_map.data,
                          data_map.data,
                          res=self.ibis_res)

    def test_find_data(self):
        """ Unit tests for find_data() function
        """

        fake_dir = './test_dir/'
        fake_dir_2 = './test_dir_2/'
        # NOTE: the yml workflow is OCCASIONALLY running out of time and
        # not removing these folders, which causes an error next time it
        # is run. For that reason, remove them first thing here .
        try:
            shutil.rmtree(fake_dir)
        except FileNotFoundError:
            pass
        try:
            shutil.rmtree(fake_dir_2)
        except FileNotFoundError:
            pass
        os.mkdir(fake_dir)
        fake_data = np.empty((1, 1))
        fake_fits_name = 'test.fits'
        fits.writeto(fake_dir + fake_fits_name, fake_data)
        found_fits = funclib.find_data(fake_dir)

        # -------- positive tests -------- :
        # Test 1: test that it can find a fits file in a directory
        self.assertEqual(fake_fits_name, found_fits[0])

        # -------- negative tests -------- :
        # Test 2: test that it doens't find files outside the scope
        # of the given dataset:
        fake_fits_name_2 = 'test2.fits'
        fits.writeto('./' + fake_fits_name_2, fake_data)
        self.assertNotIn(fake_fits_name_2, funclib.find_data(fake_dir))
        os.remove('./' + fake_fits_name_2)

        # ------ error raising tests ------ :
        # Test 3: test that it errors if filepath passed doesn't include data
        os.mkdir(fake_dir_2)
        self.assertRaises(OSError, funclib.find_data, fake_dir_2)

        shutil.rmtree(fake_dir)
        shutil.rmtree(fake_dir_2)

    def test_kmeans(self):
        """ Unit tests for test_kmeans() function
        """

        # -------- positive tests -------- :
        # Test 1: check that returns numpy array of same shape as input
        N = 10
        array_to_be_clustered = np.ones((N, N))
        # give some fake different values, so kmeans has something to cluster,
        array_to_be_clustered[0, 0] = 1
        array_to_be_clustered[0, 1] = 2
        clustered_array = funclib.kmeans_cluster(array_to_be_clustered,
                                                 llambda_axis=-1)
        self.assertEqual(np.shape(clustered_array)[0], N)

        # -------- negative tests -------- :
        # Test 2: test that the returned labels don't contian a
        # label they shouldn't (should only be 0 or 1)
        non_label = 3
        count_non_label_in_cluster =\
            np.count_nonzero(clustered_array[clustered_array == non_label])
        self.assertEqual(count_non_label_in_cluster, 0)

        # ------ error raising tests ------ :
        # Test 3: should error if passed in data of wrong shape
        self.assertRaises(Exception,
                          funclib.kmeans_cluster,
                          array_to_be_clustered,
                          3)

    def test_cross_correlation(self):
        """ Unit tests for cross_correlation() function
        """

        # -------- positive tests -------- :
        # Test 1: if arrays agree, return 0:
        test_size = 10
        test_array_1 = np.ones((test_size, test_size))
        test_array_2 = np.ones((test_size, test_size))
        test_array_1[0, 0] = 0
        test_array_2[0, 0] = 0
        self.assertEqual(0, funclib.cross_correlation(test_array_1,
                                                      test_array_2)[0])

        # Test 2: if cross correlation too low, return -1:
        test_array_1 = np.ones((test_size, test_size))
        test_array_1[0, 0] = 0
        test_array_2 = np.zeros((test_size, test_size))
        test_array_2[0, 0] = 1

        self.assertEqual(-1, funclib.cross_correlation(test_array_1,
                                                       test_array_2)[0])

        # -------- negative tests -------- :
        # Test 1: check cross correlation isn't greater than 100% or
        # less than 0%
        self.assertFalse(funclib.cross_correlation(test_array_1,
                                                   test_array_2)[1] > 1)
        self.assertFalse(funclib.cross_correlation(test_array_1,
                                                   test_array_2)[1] < 0)

        # ------ error raising tests ------ :
        # Test 1: error if no granules or intergranules in skimage cluster
        test_array_1 = np.ones((test_size, test_size))
        test_array_2 = np.ones((test_size, test_size))
        self.assertRaises(Exception, funclib.cross_correlation,
                          test_array_1, test_array_2)

    def overplot_velocities(self):
        """ Unit tests for overplot_velocities() function
        """

        data_map = funclib.sav_to_map(self.ibis_testfile, self.test_band)
        segmented_map = funclib.segment(self.ibis_fileid,
                                        data_map,
                                        self.test_method,
                                        True,
                                        'test_output/',
                                        self.ibis_res)
        out_file_path = 'test_outputs/velocity_comparison.png'
        funclib.overplot_velocities(segmented_map,
                                    self.ibis_testsfile,
                                    out_file_path)

        # -------- positive tests -------- :
        # Test 1: check that produces output plot in expected location
        self.assertTrue(os.path.exists(out_file_path))

        # -------- negative tests -------- :
        # Test 2: check that path exists
        self.assertFalse(not os.path.exists(out_file_path))

        # ------ error raising tests ------ :
        # Test 3: check that errors for incorrect input
        self.assertRaises(Exception, funclib.overplot_velocities,
                          segmented_map, cls.dkist_testfile, out_file_path)
        self.assertRaises(Exception, funclib.overplot_velocities,
                          [[1, 2, 3], [4, 5, 6]], out_file_path)
