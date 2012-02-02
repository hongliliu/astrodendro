import unittest
from astrodendro import Dendrogram
from astrodendro.components import Leaf
import numpy as np
import sys

class TestRecursionLimit(unittest.TestCase):
    """
    Test that we can efficiently compute deep dendrogram trees
    without hitting the recursion limit.
    Note: plot() uses recursion but we should be able to compute
    dendrograms without using deep recursion, which will ensure
    deep dendrograms get computed quickly.
    """
    def setUp(self):
        self._oldlimit = sys.getrecursionlimit()
        sys.setrecursionlimit(100) # Reduce recursion limit dramatically (default is 1000)
        size = 10000 # number of leaves desired in the dendrogram
        data1 = np.arange(size*2) # first row
        data2 = np.arange(size*2) # second row
        data2[::2] += 2;
        data1[-1] = 0 # set the last pixel in the first row to zero, to trigger a deep ancestor search
        self.data = np.vstack((data1,data2))
        self.size = size
        # self.data now looks like this:
        # [[ 0, 1, 2, 3, 4, 5, ...],
        #  [ 2, 1, 4, 3, 6, 5, ...]]
        # Notice every second column has a local maximum
        # so there are [size] local maxima in the array
    
    def test_recursionlimit(self):
        d = Dendrogram(self.data, verbose=False)
        
        leaves = d.get_leaves()
        
        self.assertEqual(len(leaves), self.size, msg="We expect {n} leaves, not {a}.".format(n=self.size, a=len(leaves)))
    
    def tearDown(self):
        sys.setrecursionlimit(self._oldlimit)

if __name__ == '__main__':
    unittest.main()
