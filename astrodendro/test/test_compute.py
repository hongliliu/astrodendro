import unittest
from astrodendro import Dendrogram
from astrodendro.components import Leaf
import numpy as np

class Test2DimensionalData(unittest.TestCase):
    def setUp(self):
        n = np.nan
        self.data = np.array([[n,n,n,n,n,n,n,n],
                              [n,4,n,n,n,n,n,n],
                              [n,n,n,1,n,n,0,5],
                              [3,n,n,2,3,2,0,n]])
    
    def test_dendrogramWithNan(self):
        d = Dendrogram(self.data, verbose=False)
        
        leaves = [item for item in d.trunk if type(item) == Leaf]
        branches = [item for item in d.trunk if item not in leaves]
        
        self.assertEqual(len(leaves), 2, msg="We expect two leaves among the lowest structures (the trunk)")
        self.assertEqual(len(branches), 1, msg="We expect one branch among the lowest structures (the trunk)")
        
        for leaf in leaves:
            self.assertEqual(len(leaf.f), 1, msg="Leaves in the trunk are only expected to contain one point")
            if leaf.f[0] == 4:
                self.assertEqual(leaf.coords[0], (1,1))
            elif leaf.f[0] == 3:
                self.assertEqual(leaf.coords[0], (3,0))
            else:
                self.fail("Invalid value of flux in one of the leaves")
                
class Test3DimensionalData(unittest.TestCase):
    def setUp(self):
        from astrodendro.test._testdata import data
        self.data = data
    
    def test_dendrogramComputation(self):
        d = Dendrogram(self.data, minimum_npix=8, minimum_delta=0.3, minimum_flux=1.4, verbose=False)
        
        num_leaves = len(d.get_leaves())
        self.assertEqual(num_leaves,55)

if __name__ == '__main__':
    unittest.main()