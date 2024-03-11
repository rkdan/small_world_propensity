import os
import unittest

import numpy as np
import pandas as pd
import scipy.io as sio

from small_world_propensity import (
    get_average_paths,
    get_avg_rad_eff,
    get_clustering_coefficient,
    make_symmetric,
    randomize_matrix,
    regular_matrix_generator,
    small_world_propensity,
)

TESTDATA_FILENAME = os.path.join(os.path.dirname(__file__), 'cat.mat')

class TestSmallWorldPropensity(unittest.TestCase):
    # use setUp to get cat matrix
    def setUp(self):
        self.cat = sio.loadmat(TESTDATA_FILENAME)['CIJctx']
        self.cat_shape = self.cat.shape

    def test_make_symmetric(self):
        result = make_symmetric(self.cat)
        self.assertIsInstance(result, np.ndarray)
        self.assertEqual(result.shape, self.cat.shape)
        self.assertTrue(np.allclose(result, result.T))


    def test_small_world_propensity_single_matrix(self):
        W = make_symmetric(self.cat)
        # test regular and binary
        result = small_world_propensity(W)
        self.assertIsInstance(result, pd.DataFrame)
        self.assertEqual(result.shape, (1, 11))
        self.assertTrue(np.all(0 <= result['SWP'].values <= 1))

        result = small_world_propensity(W, bin=True)
        self.assertIsInstance(result, pd.DataFrame)
        self.assertEqual(result.shape, (1, 11))
        self.assertTrue(np.all(0 <= result['SWP'].values <= 1))


    def test_small_world_propensity_list_of_matrices(self):
        W1 = make_symmetric(self.cat)
        W2 = make_symmetric(self.cat)
        result = small_world_propensity([W1, W2], bin=[True, True])
        self.assertIsInstance(result, pd.DataFrame)
        self.assertEqual(result.shape, (2, 11))
        self.assertTrue(np.all(0 <= result['SWP']) and np.all(result['SWP'] <= 1))

class TestNetworkMeasures(unittest.TestCase):
    def test_get_avg_rad_eff(self):
        W = np.array([[0, 1, 0], [1, 0, 1], [0, 1, 0]])
        result = get_avg_rad_eff(W)
        self.assertIsInstance(result, int)
        self.assertEqual(result, 1)

    def test_get_average_paths(self):
        W = np.array([[0, 1, 0], [1, 0, 1], [0, 1, 0]])
        result = get_average_paths(W)
        self.assertIsInstance(result, float)
        self.assertGreater(result, 0)

    def test_get_clustering_coefficient(self):
        W = np.array([[0, 1, 0], [1, 0, 1], [0, 1, 0]])
        result = get_clustering_coefficient(W)
        self.assertIsInstance(result, float)
        self.assertTrue(0 <= result <= 1)

class TestMatrixOperations(unittest.TestCase):
    def test_randomize_matrix(self):
        W = np.array([[0, 1, 0], [1, 0, 1], [0, 1, 0]])
        result = randomize_matrix(W)
        self.assertIsInstance(result, np.ndarray)
        self.assertEqual(result.shape, W.shape)
        self.assertTrue(np.allclose(result, result.T))

    def test_regular_matrix_generator(self):
        W = np.array([[0, 1, 0], [1, 0, 1], [0, 1, 0]])
        result = regular_matrix_generator(W, 1)
        self.assertIsInstance(result, np.ndarray)
        self.assertEqual(result.shape, W.shape)
        self.assertTrue(np.allclose(result, result.T))

    def test_make_symmetric(self):
        W = np.array([[0, 1, 0], [0, 0, 1], [1, 0, 0]])
        result = make_symmetric(W)
        self.assertIsInstance(result, np.ndarray)
        self.assertEqual(result.shape, W.shape)
        self.assertTrue(np.allclose(result, result.T))

if __name__ == '__main__':
    unittest.main()