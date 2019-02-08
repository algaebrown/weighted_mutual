import unittest
import numpy as np
from weighted_mutual.entropy import *

class TestEntropy(unittest.TestCase):

    def test_entropy(self):
        p = np.array([0.5,0.5,0])
        row = np.array([1,1,1])

        assert np.allclose(weighted_entropy(row, p), np.log2(1))

        # test another case
        row = np.array([1,0,1])
        assert np.allclose(weighted_entropy(row, p), -np.log2(0.5))

    def test_joint_entropy(self):
        p = np.array([0.5, 0.2, 0.3])
        row1 = np.array([1,0,1])
        row2 = np.array([1,1,0])

        assert np.allclose(weighted_joint_entropy(row1, row2, p), -(0.5*np.log2(0.5) + 0.2*np.log2(0.2) + 0.3*np.log2(0.3)))

        # another case: exactly the same protein
        row2 = np.array([1,0,1])
        assert np.allclose(weighted_joint_entropy(row1, row2, p), -(0.8*np.log2(0.8) + 0.2*np.log2(0.2)))

if __name__ == '__main__':
        unittest.main()
