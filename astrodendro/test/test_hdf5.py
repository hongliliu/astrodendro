import unittest
from astrodendro import Dendrogram
from astrodendro.components import Branch
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
        
        self.assertEqual(len(d.items_dict), len(d2.items_dict))
        
        np.testing.assert_array_equal(d.data, d2.data) # Do we recover the data exactly?
                         
        for idx in d2.items_dict:
            item1, item2 = d.items_dict[idx], d2.items_dict[idx]
            self.assertItemsEqual(item1.coords, item2.coords)
            self.assertItemsEqual(item1.f, item2.f)
            self.assertEqual(type(item1), type(item2))
            if type(item2) == Branch:
                self.assertEqual(item1.merge_level, item2.merge_level)        

if __name__ == '__main__':
    unittest.main()
