import unittest
from astrodendro import Dendrogram
import numpy as np
import os

class TestHDF5(unittest.TestCase):
    def setUp(self):
        n = np.nan
        self.data = np.array([[[n,n,n,n,n,n,n,n],
                               [n,4,n,n,n,n,n,n],
                               [n,n,n,1,n,n,0,5],
                               [3,n,n,2,3,2,0,n]],
                              [[n,n,n,n,n,n,n,n],
                               [1,n,n,n,n,n,n,n],
                               [1,n,1,1,0,n,0,1],
                               [2,n,n,1,3,1,n,1]],
                              [[n,2,3,4,n,n,n,n],
                               [1,1,n,n,n,n,n,n],
                               [n,n,n,n,n,n,0,1],
                               [n,n,n,1,0,1,0,n]]])
        self.test_filename = 'test.hdf5'
    
    def tearDown(self):
        os.remove(self.test_filename)
    
    def test_write(self):
        d = Dendrogram(self.data, verbose=False)
        d.to_hdf5(self.test_filename)
    
    def test_read(self):
        d = Dendrogram(self.data, verbose=False)
        d.to_hdf5(self.test_filename)
        d2 = Dendrogram()
        d2.from_hdf5(self.test_filename)
        
        self.assertEqual(len(d.get_leaves()), len(d2.get_leaves()))
        
        np.testing.assert_array_equal(d.data, d2.data) # Do we recover the data exactly?
                         
        for leaf2 in d2.get_leaves():
            matches = [leaf for leaf in d.get_leaves() if leaf.idx == leaf2.idx]
            self.assertEqual(len(matches), 1)
            leaf = matches[0]
            self.assertItemsEqual(leaf.coords, leaf2.coords)
            self.assertItemsEqual(leaf.f, leaf2.f)
        
        # TODO: Test that merge levels are preserved
                

if __name__ == '__main__':
    unittest.main()
