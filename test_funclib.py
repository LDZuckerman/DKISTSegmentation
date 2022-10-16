import unittest
import funclib
import random
import numpy as np


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

    def test_sav_to_maps(cls):

        """ Unit tests for () function
        """

        # positive tests:
        # cls.assertEqual()

        # negative test:
        # cls.assertNotEqual()

        # error raising tests
        # cls.assertRaises(Exception, funclib.sav_to_maps('ABC.txt'))
        # cls.assertRaises(Exception, funclib.sav_to_maps(''))

    def test_sav_to_numpy(cls):

        """ Unit tests for () function
        """

        # positive tests:
        # cls.assertEqual()

        # negative test:
        # cls.assertNotEqual()

        # error raising tests
        cls.assertRaises(Exception, funclib.sav_to_numpy(cls.testfile,
                                                         cls.test_instrument,
                                                         cls.test_band),
                         'ABC.txt')
        cls.assertRaises(Exception, funclib.sav_to_numpy(cls.testfile,
                                                         cls.test_instrument,
                                                         cls.test_band),
                         '')
