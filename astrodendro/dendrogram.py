# Computing Astronomical Dendrograms
# Copyright (c) 2011 Thomas P. Robitaille and Braden MacDonald
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
from .components import Branch, Leaf
from .newick import parse_newick
from .progressbar import AnimatedProgressBar
try:
    import matplotlib.pylab
    from .plot import RecursiveSortPlot, SpatialCoordPlot, SpatialMeanCoordPlot
    _plotting_enabled = True
except ImportError:
    # The plot method won't work without matplotlib, but everything else will be fine
    _plotting_enabled = False
    pass


class Dendrogram(object):

    def __init__(self, data=None, minimum_flux=-np.inf, minimum_npix=0, minimum_delta=0, verbose=True, compute=True):

        if data is not None:
            self.data = data
            self.n_dim = len(data.shape)
            self._params_used = {'min_flux':None, 'min_npix': None, 'min_delta':None}
            if compute:
                self.compute(minimum_flux, minimum_npix, minimum_delta, verbose)


    def compute(self, minimum_flux, minimum_npix, minimum_delta, verbose=False):
        
        # For reference, store the parameters used:
        self._params_used = {'min_flux':minimum_flux, 'min_npix': minimum_npix, 'min_delta':minimum_delta}

        # Create a list of all points in the cube above minimum_flux
        keep = self.data.ravel() > minimum_flux
        flux_values = self.data.ravel()[keep]
        coords = np.array(np.unravel_index( np.arange(self.data.size)[keep] , self.data.shape)).transpose()
        
        if verbose:
            print("Generating dendrogram using {:,} of {:,} pixels ({}% of data)".format(flux_values.size, self.data.size, (100*flux_values.size/self.data.size)))
            progress_bar = AnimatedProgressBar(end=flux_values.size, width=40, fill='=', blank=' ')
            
        # Define index array indicating what item each cell is part of
        # We expand each dimension by one, so the last value of each
        # index (accessed with e.g. [nx,#,#] or [-1,#,#]) is always zero
        # This permits an optimization below when finding adjacent items
        self.index_map = np.zeros(np.add(self.data.shape, 1), dtype=np.int32)

        # Dictionary of currently-defined items:
        items = {}
        
        # Define a list of offsets we add to any coordinate to get the coords
        # of all neighbouring pixels
        if self.n_dim is 3:
            neighbour_offsets = np.array([(0,0,-1),(0,0,1),(0,-1,0),(0,1,0),(-1,0,0),(1,0,0)])
        else:
            neighbour_offsets = np.array([(0,-1),(0,1),(-1,0),(1,0)])

        # Loop from largest to smallest flux value. Each time, check if the 
        # pixel connects to any existing leaf. Otherwise, create new leaf.
        
        count = 0

        for i in np.argsort(flux_values)[::-1]:
            
            def next_idx():
                return i+1
                # Generate IDs index i. We add one to avoid ID 0
            
            flux = flux_values[i]
            coord = tuple(coords[i])
            
            # Print stats
            count += 1
            if verbose and (count % 100 == 0):
                progress_bar + 100
                progress_bar.show_progress()

            # Check if point is adjacent to any leaf
            # We don't worry about the edges, because overflow or underflow in 
            # any one dimension will always land on an extra "padding" cell 
            # with value zero added above when index_map was created
            indices_adjacent = [tuple(c) for c in np.add(neighbour_offsets, coord)]
            adjacent = [self.index_map[c] for c in indices_adjacent if self.index_map[c] != 0]
            
            # Replace adjacent elements by its ancestor
            adjacent = [ items[a].ancestor.idx for a in adjacent]

            # Remove duplicates
            adjacent = list(set(adjacent))

            # Find how many unique adjacent structures there are
            n_adjacent = len(adjacent)

            if n_adjacent == 0:  # Create new leaf

                # Set absolute index of the new element
                idx = next_idx()

                # Create leaf
                leaf = Leaf(coord, flux, idx=idx)

                # Add leaf to overall list
                items[idx] = leaf

                # Set absolute index of pixel in index map
                self.index_map[coord] = idx

            elif n_adjacent == 1:  # Add to existing leaf or branch

                # Get absolute index of adjacent element
                idx = adjacent[0]

                # Get adjacent item
                item = items[idx]

                # Add point to item
                item.add_point(coord, flux)

                # Set absolute index of pixel in index map
                self.index_map[coord] = idx

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

                # Now, how many significant adjacent items are left?

                if len(adjacent) == 0:

                    # There are no separate leaves left (and no branches), so pick the
                    # first one as the reference and merge all the others onto it

                    idx = merge[0]
                    leaf = items[idx]

                    # Add current point to the leaf
                    leaf.add_point(coord, flux)

                    # Set absolute index of pixel in index map
                    self.index_map[coord] = idx

                    for i in merge[1:]:

                        # print "Merging leaf %i onto leaf %i" % (i, idx)

                        # Remove leaf
                        removed = items.pop(i)

                        # Merge old leaf onto reference leaf
                        leaf.merge(removed)

                        # Update index map
                        removed.add_footprint(self.index_map, idx)

                elif len(adjacent) == 1:
                    
                    # There is one significant adjacent leaf/branch left.
                    # Add the point under consideration and all insignificant
                    # leaves in 'merge' to the adjacent leaf/branch

                    idx = adjacent[0]
                    item = items[idx]  # Could be a leaf or a branch

                    # Add current point to the leaf/branch
                    item.add_point(coord, flux)

                    # Set absolute index of pixel in index map
                    self.index_map[coord] = idx

                    for i in merge:

                        # print "Merging leaf %i onto leaf/branch %i" % (i, idx)

                        # Remove leaf
                        removed = items.pop(i)

                        # Merge insignificant leaves onto the leftover leaf/branch
                        item.merge(removed)

                        # Update index map
                        removed.add_footprint(self.index_map, idx)

                else:

                    # Set absolute index of the new element
                    idx = next_idx()

                    # Create branch
                    branch = Branch([items[j] for j in adjacent], \
                                    coord, flux, idx=idx)

                    # Add branch to overall list
                    items[idx] = branch

                    # Set absolute index of pixel in index map
                    self.index_map[coord] = idx

                    for i in merge:

                        # print "Merging leaf %i onto branch %i" % (i, idx)

                        # Remove leaf
                        removed = items.pop(i)

                        # Merge old leaf onto reference leaf
                        branch.merge(removed)

                        # Update index map
                        removed.add_footprint(self.index_map, idx)


        if verbose:
            progress_bar.progress = 100 # Done
            progress_bar.show_progress()
            print("") # newline

        # Remove orphan leaves that aren't large enough
        remove = [idx for idx,item in items.iteritems()
                  if item.parent is None and type(item) == Leaf
                  and (item.npix < minimum_npix or item.fmax - item.fmin < minimum_delta)]
        for idx in remove:
            removed = items.pop(idx)
            removed.add_footprint(self.index_map, 0)
    
        # Save a list of all items accessible by ID
        self.items_dict = items

        # Create trunk from objects with no ancestors
        self.trunk = [item for item in items.itervalues() if item.parent == None] 

    @property
    def all_items(self):
        " Return a flattened iterable containing all items in the dendrogram "
        return self.items_dict.itervalues()
    
    @property
    def main_branch(self):
        " Return the branch (or leaf!) in the trunk with the most children "
        if len(self.trunk) == 0:
            return None
        return sorted(self.trunk, key=lambda i: len(i.descendants) if type(i) == Branch else 1)[-1]

    def get_leaves(self):
        " Return a flattened list of all leaves in the dendrogram "
        return [i for i in self.items_dict.itervalues() if type(i) == Leaf]

    def to_newick(self):
        return "(%s);" % ','.join([item.newick for item in self.trunk])

    def to_hdf5(self, filename):

        import h5py

        f = h5py.File(filename, 'w')

        f.attrs['n_dim'] = self.n_dim

        f.create_dataset('newick', data=self.to_newick())

        d = f.create_dataset('index_map', data=self.index_map, compression=True)
        d.attrs['CLASS'] = 'IMAGE'
        d.attrs['IMAGE_VERSION'] = '1.2'
        d.attrs['IMAGE_MINMAXRANGE'] = [self.index_map.min(), self.index_map.max()]
        
        # Make map of leaves vs branches
        item_type_map = np.zeros(self.data.shape, dtype=np.uint8)
        for idx in self.items_dict:
            item = self.items_dict[idx]
            if type(item) == Leaf:
                item.add_footprint(item_type_map, 2)
            else:
                item.add_footprint(item_type_map, 1, recursive=False)

        d = f.create_dataset('item_type_map', data=item_type_map, compression=True)
        d.attrs['CLASS'] = 'IMAGE'
        d.attrs['IMAGE_VERSION'] = '1.2'
        d.attrs['IMAGE_MINMAXRANGE'] = [item_type_map.min(), item_type_map.max()]

        d = f.create_dataset('data', data=self.data, compression=True)
        d.attrs['CLASS'] = 'IMAGE'
        d.attrs['IMAGE_VERSION'] = '1.2'
        d.attrs['IMAGE_MINMAXRANGE'] = [self.data.min(), self.data.max()]

        f.close()

    def from_hdf5(self, filename):

        import h5py

        f = h5py.File(filename, 'r')

        self.n_dim = f.attrs['n_dim']

        self.data = f['data'].value
        self.index_map = f['index_map'].value
        self.item_type_map = f['item_type_map'].value
        self.items_dict = {}

        tree = parse_newick(f['newick'].value)

        def construct_tree(d):
            items = []
            for idx in d:
                item_coords = zip(*(np.where(self.index_map == idx)))
                f = self.data[self.index_map == idx]
                if type(d[idx]) == tuple:
                    sub_items_repr = d[idx][0] # Parsed representation of sub items
                    sub_items = construct_tree(sub_items_repr)
                    for i in sub_items:
                        self.items_dict[i.idx] = i
                    b = Branch(sub_items, item_coords[0], f[0], idx=idx)
                    for i in range(1, len(f)):
                        b.add_point(item_coords[i], f[i])
                    # Correct merge levels - complicated because of the
                    # order in which we are building the tree.
                    # What we do is look at the heights of this branch's
                    # 1st child as stored in the newick representation, and then
                    # work backwards to compute the merge level of this branch
                    first_child_repr = sub_items_repr.itervalues().next()
                    if type(first_child_repr) == tuple:
                        height = first_child_repr[1]
                    else:
                        height = first_child_repr
                    b.merge_level = sub_items[0].fmax - height
                    self.items_dict[idx] = b
                    items.append(b)
                else:
                    l = Leaf(item_coords[0], f[0], idx=idx)
                    for i in range(1, len(f)):
                        l.add_point(item_coords[i], f[i])
                    items.append(l)
                    self.items_dict[idx] = l
            return items

        self.trunk = construct_tree(tree)
    
    def item_at(self, coords):
        " Get the item at the given pixel coordinate, or None "
        idx = self.index_map[coords]
        if idx:
            return self.items_dict[idx]
        return None
    
    @property
    def min_flux(self):
        return self._params_used['min_flux']
    @property
    def min_npix(self):
        return self._params_used['min_npix']
    @property
    def min_delta(self):
        return self._params_used['min_delta']
    
    def plot(self, **kwargs):
        """
        Plot a Dendrogram using matplotlib.
        Works with IPython's pylab mode, so you can just type: dendrogram.plot()
        Returns a DendrogramPlot object 
        """
        if not _plotting_enabled:
            raise Exception("Plotting is not available without matplotlib.")
        axes = matplotlib.pylab.gca()
        plot = self.make_plot(axes, **kwargs)
        matplotlib.pylab.draw_if_interactive()
        return plot

    def make_plot(self, axes, style='rfmax', color=(0,0,0.2,1), line_width=1, autoscale=True, **kwargs):
        """
        Plot the dendrogram in 2D on the given matplotlib axes.
        
        style must be one of the following:
            'rfmax': sort the dendrogram by maximum flux found among all child items.
            'x_coord': plot items according to the X coordinate of the leaf item with peak flux value
            'y_coord': same, but uses Y coordinate
            'z_coord': same, but uses Y coordinate
        
        color can be:
            1. An RGBA tuple, e.g. for red color=(1,0,0,1)
            2. "npix" to color such that bigger items are red,
               smaller ones black.
            3. "fsum" to color such that bigger items are red,
               smaller ones black, but using the sum of flux values
               rather than just the number of pixels.
            4. Any lambda method or callable to compute the color for each item.
               The method must return an RGBA tuple,
        
        This returns a DendrogramPlot object which has two main methods:
            create_highlighter(color='red', alpha=1, line_width_extra=1):
                Create a highlighter object that can be used to color in items
                on the plot with different colors. The returned object will
                have a .highlight(item) method.
            
            item_at(x,y):
                Return the item (or None) associated with the given 2D coordinate on the axes.
        """
        if not _plotting_enabled:
            raise Exception("Plotting is not available without matplotlib.")
        if color == "npix":
            max_npix = max(i.npix_self for i in self.all_items)
            color_lambda = lambda item: (1*(float(item.npix_self)/max_npix),0,0,1)
        elif color == "fsum":
            max_fsum = max(i.f_sum_self for i in self.all_items)
            color_lambda = lambda item: (1*(float(item.f_sum_self)/max_fsum),0,0,1)
        elif not hasattr(color, '__call__'):
            color_lambda = lambda _: color
        else:
            color_lambda = color
        if style == 'rfmax':
            return RecursiveSortPlot(dendrogram=self, axes=axes, color_lambda=color_lambda, line_width=line_width, autoscale=autoscale, **kwargs)
        elif style[1:] == '_coord':
            coord_indices = {'x':0, 'y':1, 'z':2}
            coord_index = coord_indices[style[0]]
            return SpatialCoordPlot(dendrogram=self, axes=axes, color_lambda=color_lambda, line_width=line_width, coord_index=coord_index, **kwargs)
        elif style[1:] == '_mean':
            coord_indices = {'x':0, 'y':1, 'z':2}
            coord_index = coord_indices[style[0]]
            return SpatialMeanCoordPlot(dendrogram=self, axes=axes, color_lambda=color_lambda, line_width=line_width, coord_index=coord_index, **kwargs)
        else:
            raise Exception("Invalid plot style requested.")
