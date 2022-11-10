import unittest
import funclib
import random
import numpy as np
import sunpy.map as sm
import scipy
import sunpy


class TestUtils(unittest.TestCase):

    @classmethod
    def setUpClass(cls):

        """ Set up for unit testing by creating toy data
        """
        cls.testfile = 'IBIS.granulation.aligned.25Apr2019.seq56.sav'
        cls.test_instrument = 'IBIS'
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

    def test_sav_to_map(self):

        """ Unit tests for sav_to_map() function
        """

        # positive tests:
        test_type = type(funclib.sav_to_map(self.testfile, self.test_band))
        self.assertEqual(test_type, sunpy.map.mapbase.GenericMap)

        # negative test:
        # TO DO: this is a very weak test right now
        self.assertNotEqual(test_type, np.ndarray)

        # error raising tests
        self.assertRaises(Exception, funclib.sav_to_map, 'ABC.txt',
                          self.test_band)
        self.assertRaises(Exception, funclib.sav_to_map, '', self.test_band)

    def test_sav_to_numpy(self):

        """ Unit tests for sav_to_numpy() function
        """

        # positive tests:
        test_type = type(funclib.sav_to_numpy(self.testfile,
                                              self.test_instrument,
                                              self.test_band))
        self.assertEqual(test_type, np.ndarray)

        # negative test:
        # TO DO: this is a very weak test right now
        self.assertEqual(test_type, np.ndarray)

        # error raising tests
        self.assertRaises(Exception, funclib.sav_to_numpy, 'ABC.txt',
                          self.test_instrument, self.test_band)

    def test_segment(self):
        """ Unit tests for segment() function
        """

        # positive tests
        # check that the returnde type is correct
        data_map = funclib.sav_to_map(self.testfile, self.test_band)
        segmented = funclib.segment(data_map, self.test_method)
        test_type = type(segmented)
        self.assertEqual(test_type, sunpy.map.mapbase.GenericMap)

        # check that the shape is unchanged
        # get an array of pixel values
        initial_pix = sm.all_pixel_indices_from_map(data_map).value
        seg_pixels = sm.all_pixel_indices_from_map(segmented).value
        self.assertEqual(seg_pixels.shape, initial_pix.shape)

        # negtaive tests
        # check that the values in the array are changed
        random.seed(42)
        # pick 10 random indices to check
        x = random.sample(list(np.arange(0, data_map.data.shape[0], 1)), 10)
        y = random.sample(list(np.arange(0, data_map.data.shape[1], 1)), 10)
        for i in range(len(x)):
            self.assertNotEqual(data_map.data[x[i], y[i]],
                                segmented.data[x[i], y[i]])

        # error raising tests
        self.assertRaises(TypeError,
                          funclib.segment,
                          data_map,
                         'skimage')

        self.assertRaises(TypeError,
                          funclib.segment,
                          funclib.sav_to_numpy(self.testfile,
                                               self.test_instrument,
                                               self.test_band), 'skimage')

    def test_get_threshold(self):
        """ Unit tests for get_threshold() function
        """

        # positive tests
        # check type of output
        test_arr1 = np.array([[1, 2, 3, 4, 5], [1, 2, 3, 4, 5]])
        threshold1 = funclib.get_threshold(test_arr1, self.test_method)
        self.assertEqual(type(threshold1), np.float64)
        # check range of output
        self.assertTrue((threshold1 > 0) * (threshold1 < np.max(test_arr1)))

        # negtaive tests
        # check that different arrays return different thresholds
        test_arr2 = np.array([[2, 3, 4, 5, 6], [2, 3, 4, 5, 6]])
        threshold2 = funclib.get_threshold(test_arr2, self.test_method)
        self.assertNotEqual(threshold1, threshold2)

        # error raising tests
        self.assertRaises(ValueError, funclib.get_threshold, test_arr1,
                         'banana')

    def test_remove_middles(self):
        """ Unit tests for remove_middles() function
        """

        # positive tests
        # check that the correct dimensions are returned
        data_map = funclib.sav_to_map(self.testfile, self.test_band)
        thresholded = np.uint8(data_map.data > np.nanmedian(data_map.data))
        self.assertEqual(thresholded.shape,
                         funclib.remove_middles(thresholded).shape)

        # negative tests
        # check that the returned array has fewer (or same number)
        # 0-valued pixels as input array (for a data set which we
        # know by eye should have some middle sections removed)
        middles_removed = funclib.remove_middles(thresholded)
        self.assertFalse(np.count_nonzero(middles_removed)
                         < np.count_nonzero(thresholded))

        # error raising tests
        # check that raises error if passed array is not binary
        self.assertRaises(ValueError, funclib.remove_middles, data_map.data)

    def test_mark_faculae(self):
        """ Unit tests for mark_faculae() function
        """

        # positive tests
        # check that the correct dimensions are returned
        data_map = funclib.sav_to_map(self.testfile, self.test_band)
        thresholded = np.uint8(data_map.data > np.nanmedian(data_map.data))
        faculae_marked = funclib.mark_faculae(thresholded, data_map.data)
        self.assertEqual(thresholded.shape,
                         faculae_marked.shape)

        # negtaive tests
        # check that the returned array has some 0.5 values (for a dataset
        # that we know has faculae by eye)
        self.assertNotEqual(len(np.where(faculae_marked == 0.5)[0]), 0)

        # error raising tests
        self.assertRaises(ValueError, funclib.mark_faculae, data_map.data,
                          data_map.data)
