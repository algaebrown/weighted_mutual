from weighted_mutual.weight import *
import unittest
import numpy as np

# test file
t = '../test_file/small_test'

class TestWeight(unittest.TestCase):

    def test_weight(self):
        assert np.allclose(weight(t), np.array([16/5, 16/6, 16/5]))

    def test_probability(self):
        assert np.allclose(probability(weight(t)).sum(), 1)


if __name__ == '__main__':
        unittest.main()
