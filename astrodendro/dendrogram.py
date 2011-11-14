# Computing Astronomical Dendrograms
# Copyright (c) 2011 Thomas P. Robitaille
#
# Permission is hereby granted, free of charge, to any person obtaining a
# copy of this software and associated documentation files (the "Software"),
# to deal in the Software without restriction, including without limitation
# the rights to use, copy, modify, merge, publish, distribute, sublicense,
# and/or sell copies of the Software, and to permit persons to whom the
# Software is furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
# DEALINGS IN THE SOFTWARE.

# Notes:
# - An item is a leaf or a branch
# - An ancestor is the largest structure that an item is part of

import numpy as np

from astrodendro.components import Trunk, Branch, Leaf
from astrodendro.newick import parse_newick
try:
    import matplotlib
    import matplotlib.pylab
except ImportError:
    # The plot method won't work without matplotlib, but everything else will be fine
    pass


class Dendrogram(object):

    def __init__(self, *args, **kwargs):

        if len(args) == 1:
            self._compute(*args, **kwargs)

        self._reset_idx()

    def _reset_idx(self):
        self._idx_counter = 0

    def _next_idx(self):
        self._idx_counter += 1
        return self._idx_counter

    def _compute(self, data, minimum_flux=-np.inf, minimum_npix=0, minimum_delta=0, verbose=True):

        # Reset ID counter
        self._reset_idx()

        # Initialize list of ancestors
        ancestor = {}

        # If array is 2D, recast to 3D
        if len(data.shape) == 2:
            self.n_dim = 2
            self.data = data.reshape(1, data.shape[0], data.shape[1])
        elif len(data.shape) == 3:
            self.n_dim = 3
            self.data = data
        else:
            raise Exception("Invalid # of dimensions")

        # Record data shape
        nz, ny, nx = self.data.shape

        # Create a list of all points in the cube above minimum_flux
        keep = self.data.ravel() > minimum_flux
        flux_values = self.data.ravel()[keep]
        coords = np.array(np.unravel_index( np.arange(nz*ny*nx)[keep] , self.data.shape)).transpose()
        
        if verbose:
            print "Number of points above minimum: %i" % len(flux_values)
            
        # Define index array indicating what item each cell is part of
        # We expand each dimension by one, so the last value of each
        # index (accessed with e.g. [nx,#,#] or [-1,#,#]) is always zero
        # This permits an optimization below when finding adjacent items
        self.index_map = np.zeros((nz+1,ny+1,nx+1), dtype=np.int32)

        # Dictionary of currently-defined items:
        items = {}

        # Loop from largest to smallest flux value. Each time, check if the 
        # pixel connects to any existing leaf. Otherwise, create new leaf.
        
        count = 0

        for i in np.argsort(flux_values)[::-1]:
            
            flux = flux_values[i]
            coord = coords[i]
            z,y,x = coord
            
            # Print stats
            if verbose and count % 10000 == 0:
                print "%i..." % count
            count += 1

            # Check if point is adjacent to any leaf
            indices_adjacent = [(z,y,x-1),(z,y,x+1),(z,y-1,x),(z,y+1,x),(z-1,y,x),(z+1,y,x)]
            adjacent = [self.index_map[c] for c in indices_adjacent if self.index_map[c] != 0]
            
            # Replace adjacent elements by its ancestor
            for j in range(len(adjacent)):
                if ancestor[adjacent[j]] is not None:
                    adjacent[j] = ancestor[adjacent[j]]

            # Remove duplicates
            adjacent = list(set(adjacent))

            # Find how many unique adjacent structures there are
            n_adjacent = len(adjacent)

            if n_adjacent == 0:  # Create new leaf

                # Set absolute index of the new element
                idx = self._next_idx()

                # Create leaf
                leaf = Leaf(x, y, z, flux, idx=idx)

                # Add leaf to overall list
                items[idx] = leaf

                # Set absolute index of pixel in index map
                self.index_map[z, y, x] = idx

                # Create new entry for ancestor
                ancestor[idx] = None

            elif n_adjacent == 1:  # Add to existing leaf or branch

                # Get absolute index of adjacent element
                idx = adjacent[0]

                # Get adjacent item
                item = items[idx]

                # Add point to item
                item.add_point(x, y, z, flux)

                # Set absolute index of pixel in index map
                self.index_map[z, y, x] = idx

            else:  # Merge leaves

                # At this stage, the adjacent items might consist of an arbitrary
                # number of leaves and branches.

                # Find all leaves that are not important enough to be kept
                # separate. These leaves will now be treated the same as the pixel
                # under consideration
                merge = []
                for idx in adjacent:
                    if type(items[idx]) == Leaf:
                        leaf = items[idx]
                        if leaf.npix < minimum_npix or leaf.fmax - flux < minimum_delta:
                            merge.append(idx)

                # Remove merges from list of adjacent items
                for idx in merge:
                    adjacent.remove(idx)

                # If there is only one item left, then if it is a leaf, add to the
                # list to merge, and if it is a branch then add the merges to the
                # branch.

                if len(adjacent) == 0:

                    # There are no separate leaves left (and no branches), so pick the
                    # first one as the reference and merge all the others onto it

                    idx = merge[0]
                    leaf = items[idx]

                    # Add current point to the leaf
                    leaf.add_point(x, y, z, flux)

                    # Set absolute index of pixel in index map
                    self.index_map[z, y, x] = idx

                    for i in merge[1:]:

                        # print "Merging leaf %i onto leaf %i" % (i, idx)

                        # Remove leaf
                        removed = items.pop(i)

                        # Merge old leaf onto reference leaf
                        leaf.merge(removed)

                        # Update index map
                        removed.add_footprint(self.index_map, idx)

                elif len(adjacent) == 1:

                    if type(items[adjacent[0]]) == Leaf:

                        idx = adjacent[0]
                        leaf = items[idx]

                        # Add current point to the leaf
                        leaf.add_point(x, y, z, flux)

                        # Set absolute index of pixel in index map
                        self.index_map[z, y, x] = idx

                        for i in merge:

                            # print "Merging leaf %i onto leaf %i" % (i, idx)

                            # Remove leaf
                            removed = items.pop(i)

                            # Merge old leaf onto reference leaf
                            leaf.merge(removed)

                            # Update index map
                            removed.add_footprint(self.index_map, idx)

                    else:

                        idx = adjacent[0]
                        branch = items[idx]

                        # Add current point to the branch
                        branch.add_point(x, y, z, flux)

                        # Set absolute index of pixel in index map
                        self.index_map[z, y, x] = idx

                        for i in merge:

                            # print "Merging leaf %i onto branch %i" % (i, idx)

                            # Remove leaf
                            removed = items.pop(i)

                            # Merge old leaf onto reference leaf
                            branch.merge(removed)

                            # Update index map
                            removed.add_footprint(self.index_map, idx)

                else:

                    # Set absolute index of the new element
                    idx = self._next_idx()

                    # Create branch
                    branch = Branch([items[j] for j in adjacent], \
                                    x, y, z, flux, idx=idx)

                    # Add branch to overall list
                    items[idx] = branch

                    # Set absolute index of pixel in index map
                    self.index_map[z, y, x] = idx

                    # Create new entry for ancestor
                    ancestor[idx] = None

                    for i in merge:

                        # print "Merging leaf %i onto branch %i" % (i, idx)

                        # Remove leaf
                        removed = items.pop(i)

                        # Merge old leaf onto reference leaf
                        branch.merge(removed)

                        # Update index map
                        removed.add_footprint(self.index_map, idx)

                    for j in adjacent:
                        ancestor[j] = idx
                        for a in ancestor:
                            if ancestor[a] == j:
                                ancestor[a] = idx

        if verbose and not count % 10000 == 0:
            print "%i..." % count

        # Remove orphan leaves that aren't large enough
        remove = []
        for idx in items:
            item = items[idx]
            if type(item) == Leaf:
                if item.npix < minimum_npix or item.fmax - item.fmin < minimum_delta:
                    remove.append(idx)
        for idx in remove:
            items.pop(idx)

        # Create trunk from objects with no ancestors
        self.trunk = Trunk()
        for idx in items:
            if ancestor[idx] is None:
                self.trunk.append(items[idx])

        # Make map of leaves vs branches
        self.item_type_map = np.zeros(self.data.shape, dtype=np.uint8)
        for idx in items:
            item = items[idx]
            if type(item) == Leaf:
                item.add_footprint(self.item_type_map, 2)
            else:
                item.add_footprint(self.item_type_map, 1, recursive=False)

        # Re-cast to 2D if original dataset was 2D
        if self.n_dim == 2:
            self.data = self.data[0, :, :]
            self.index_map = self.index_map[0, :, :]
            self.item_type_map = self.item_type_map[0, :, :]

    def get_leaves(self):
        return self.trunk.get_leaves()

    def to_newick(self):
        return self.trunk.to_newick()

    def to_hdf5(self, filename):

        import h5py

        f = h5py.File(filename, 'w')

        f.attrs['n_dim'] = self.n_dim

        f.create_dataset('newick', data=self.to_newick())

        d = f.create_dataset('index_map', data=self.index_map, compression=True)
        d.attrs['CLASS'] = 'IMAGE'
        d.attrs['IMAGE_VERSION'] = '1.2'
        d.attrs['IMAGE_MINMAXRANGE'] = [self.index_map.min(), self.index_map.max()]

        d = f.create_dataset('item_type_map', data=self.item_type_map, compression=True)
        d.attrs['CLASS'] = 'IMAGE'
        d.attrs['IMAGE_VERSION'] = '1.2'
        d.attrs['IMAGE_MINMAXRANGE'] = [self.item_type_map.min(), self.item_type_map.max()]

        d = f.create_dataset('data', data=self.data, compression=True)
        d.attrs['CLASS'] = 'IMAGE'
        d.attrs['IMAGE_VERSION'] = '1.2'
        d.attrs['IMAGE_MINMAXRANGE'] = [self.data.min(), self.data.max()]

        f.close()

    def from_hdf5(self, filename):

        import h5py

        f = h5py.File(filename, 'r')

        self.n_dim = f.attrs['n_dim']

        # If array is 2D, reshape to 3D
        if self.n_dim == 2:
            self.data = f['data'].value.reshape(1, f['data'].shape[0], f['data'].shape[1])
            self.index_map = f['index_map'].value.reshape(1, f['data'].shape[0], f['data'].shape[1])
            self.item_type_map = f['item_type_map'].value.reshape(1, f['data'].shape[0], f['data'].shape[1])
        else:
            self.data = f['data'].value
            self.index_map = f['index_map'].value
            self.item_type_map = f['item_type_map'].value


        # Create arrays with pixel positions
        x = np.arange(self.data.shape[2], dtype=np.int32)
        y = np.arange(self.data.shape[1], dtype=np.int32)
        z = np.arange(self.data.shape[0], dtype=np.int32)
        X, Y, Z = meshgrid_nd(x, y, z)

        tree = parse_newick(f['newick'].value)

        def construct_tree(d):
            items = []
            for idx in d:
                x = X[self.index_map == idx]
                y = Y[self.index_map == idx]
                z = Z[self.index_map == idx]
                f = self.data[self.index_map == idx]
                if type(d[idx]) == tuple:
                    sub_items = construct_tree(d[idx][0])
                    b = Branch(sub_items, x[0], y[0], z[0], f[0], idx=idx)
                    for i in range(1, len(x)):
                        b.add_point(x[i], y[i], z[i], f[i])
                    items.append(b)
                else:
                    l = Leaf(x[0], y[0], z[0], f[0], idx=idx)
                    for i in range(1, len(x)):
                        l.add_point(x[i], y[i], z[i], f[i])
                    items.append(l)
            return items

        self.trunk = Trunk()
        for item in construct_tree(tree):
            self.trunk.append(item)

        # Re-cast to 2D if original dataset was 2D
        if self.n_dim == 2:
            self.data = self.data[0, :, :]
            self.index_map = self.index_map[0, :, :]
            self.item_type_map = self.item_type_map[0, :, :]
    
    def plot(self, line_width = 1, spacing = 5):
        axis = matplotlib.pylab.gca()
        plot = self.trunk.plot_dendrogram(line_width, spacing)
        axis.set_xlim([plot.xmin, plot.xmax]) 
        axis.set_ylim([plot.ymin, plot.ymax])
        axis.set_xticks([])
        axis.set_xticklabels([])
        if line_width > 1:
            # Y values will not be correct, so hide them:
            axis.set_yticks([])
            axis.set_yticklabels([])
        line_collection = matplotlib.collections.LineCollection(plot.lines, linewidths = line_width)
        axis.add_collection(line_collection)
        matplotlib.pylab.draw_if_interactive()
