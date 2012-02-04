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
        
        ########################################
        # Check the trunk elements:
        
        leaves = [item for item in d.trunk if type(item) == Leaf]
        branches = [item for item in d.trunk if item not in leaves]
        
        self.assertEqual(len(leaves), 2, msg="We expect two leaves among the lowest structures (the trunk)")
        self.assertEqual(len(branches), 1, msg="We expect one branch among the lowest structures (the trunk)")
        
        for leaf in leaves:
            self.assertEqual(len(leaf.f), 1, msg="Leaves in the trunk are only expected to contain one point")
            self.assertIsNone(leaf.parent)
            self.assertEqual(leaf.ancestor, leaf)
            if leaf.f[0] == 4:
                self.assertEqual(leaf.coords[0], (1,1))
            elif leaf.f[0] == 3:
                self.assertEqual(leaf.coords[0], (3,0))
            else:
                self.fail("Invalid value of flux in one of the leaves")
        
        ########################################
        # Check all the structures on the branch:
        branch = branches[0]
        self.assertIsNone(branch.parent)
        self.assertEqual(branch.ancestor, branch)
        self.assertEqual(branch.merge_level, 0)
        
        self.assertEqual(len(branch.items), 2)
        for leaf in branch.items:
            self.assertIsInstance(leaf, Leaf)
            self.assertEqual(leaf.ancestor, branch)
            self.assertEqual(leaf.parent, branch)
            if 5 in leaf.f:
                self.assertEqual(leaf.f_sum, 5)
            elif 3 in leaf.f:
                self.assertEqual(leaf.f_sum, 1+2+3+2)
            else:
                self.fail("Invalid child of the branch")
        
                
class Test3DimensionalData(unittest.TestCase):
    def setUp(self):
        from astrodendro.test._testdata import data
        self.data = data
    
    def test_dendrogramComputation(self):
        d = Dendrogram(self.data, minimum_npix=8, minimum_delta=0.3, minimum_flux=1.4, verbose=False)
        
        num_leaves = len(d.get_leaves())
        self.assertEqual(num_leaves,55) # This data with these parameters should produce 55 leaves
        
        # Now check every pixel in the data cube (this takes a while).
        # The following loop construct may look crazy, but it is a more
        # efficient way of iterating through the array than using a regular
        # nditer with multi_index.
        for coord in [tuple(c) for c in np.array(np.unravel_index( np.arange(self.data.size), self.data.shape)).transpose()]:
            f = self.data[coord]
            if (f < 1.4):
                self.assertEqual(d.item_at(coord), None)
            else:
                item = d.item_at(coord)
                if item:
                    # The current pixel is associated with part of the dendrogram.
                    self.assertIn(coord, item.coords, "Pixel at {0} is claimed to be part of {1}, but that item does not contain the coordinate {0}!".format(coord, item))
                    fmax_item, _, fmax = item.get_peak_recursive()
                    if fmax_item is item:
                        # The current pixel is the peak pixel in this item
                        pass
                    else:
                        self.assertTrue(fmax >= f)

if __name__ == '__main__':
    unittest.main()