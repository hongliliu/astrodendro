import numpy as np
import matplotlib as mpl

def item_sort_key(item):
    " Given a Branch or Leaf, this returns a key used to sort them for plotting "
    return item.height # TODO: test recursive-fmax, npix, height, f_sum

class Leaf(object):

    def __init__(self, coord, f, idx=None):
        self.coords = [coord]
        self.f = [f]
        self.fmin, self.fmax = f, f
        self.idx = idx
        self.parent = None

    @property
    def npix(self):
        return len(self.f)
    @property
    def f_sum(self):
        return sum(self.f)
    @property
    def height(self):
        if self.parent == None:
            return self.fmax - self.fmin
        else:
            return self.fmax - self.parent.merge_level
    
    @property
    def ancestor(self):
        if self.parent == None:
            return self
        else:
            return self.parent.ancestor

    def add_point(self, coord, f):
        "Add point to current leaf"
        self.coords.append(coord)
        self.f.append(f)
        self.fmin, self.fmax = min(f, self.fmin), max(f, self.fmax)

    def merge(self, leaf):
        self.coords.extend(leaf.coords)
        self.f.extend(leaf.f)
        self.fmin, self.fmax = min(leaf.fmin, self.fmin), max(leaf.fmax, self.fmax)

    def add_footprint(self, image, level):
        "Fill in a map which shows the depth of the tree"
        for c in self.coords:
            image[c] = level

    def to_newick(self):
        return "%i:%.3f" % (self.idx, self.height)

    def get_peak(self):
        imax = np.argmax(self.f)
        return self.coords[imax], self.f[imax]
    
    def get_peak_recursive(self): # Preserve same interface as Branch, below
        c,f = self.get_peak()
        return self, c, f


class Branch(Leaf):

    def __init__(self, items, coord, f, idx=None):
        self.merge_level = f # Record the exact flux level that triggered creation of this branch
        self.items = sorted(items, key=item_sort_key)
        for item in items:
            item.parent = self
        Leaf.__init__(self, coord, f, idx=idx)

    @property
    def npix(self):
        " Number of pixels among this branch and it's children "
        return len(self.f) + self.npix_children
    @property
    def npix_children(self):
        " Number of pixels among this branch's children "
        return np.sum([item.npix for item in self.items])
    @property
    def npix_self(self):
        " Number of pixels in this Branch, not counting children "
        return len(self.f)
    
    @property
    def f_sum(self):
        return np.sum(self.f) + self.f_sum_children
    @property
    def f_sum_children(self):
        return np.sum([item.f_sum for item in self.items])
    @property
    def f_sum_self(self):
        return np.sum(self.f)


    def add_footprint(self, image, level, recursive=True):
        if recursive:
            for item in self.items:
                item.add_footprint(image, level + 1)
        Leaf.add_footprint(self, image, level)

    def to_newick(self):
        newick_items = []
        for item in self.items:
            newick_items.append(item.to_newick())
        return "(%s)%s:%.3f" % (','.join(newick_items), self.idx, self.height)
    
    def get_peak_recursive(self):
        """
        Return the item, coordinate, and flux value of the pixel with highest
        intensity among this branch and its children.
        """
        imax = np.argmax(self.f)
        item = self
        c, f = self.coords[imax], self.f[imax]
        for child in self.items:
            # get Child Item, Child peak Coords, Child Flux:
            ci, cc, cf = child.get_peak_recursive() 
            if cf > f:
                item, c, f = ci, cc, cf
        return item, c, f


class DendrogramPlot():
    def __init__(self, trunk, axes=None, color='blue', line_width=1, spacing=5):
        self.line_width = line_width
        self.spacing = spacing
        self.x_increment = line_width + spacing
        self.flux_mult = line_width # All y values (flux values) should be multiplied by this factor 
        self._next_x = 0
        self.lines = [] # list of lists of (x,y) tuples
        
        # the xmin, xmax, ymin, ymax properties give the size of the canvas onto which the lines are drawn
        self.ymin = np.min([item.fmin for item in trunk]) * self.flux_mult
        self.ymax = self.ymin # ymax will get calculated during _plot_item calls
        # Create a dictionary to find what item is found at a given point.
        # Keys to this dictionary are (x1,x2,y1,y2) tuples defining the rect
        # in which the Branch or Leaf is plotted. Values in the dictionary
        # are references to the actual Branch or Leaf
        self._item_map = {}
        self._x_map = {} # Map between item idx and the value of cur_x when the item was created

        # Do the plotting recursively:
        for item in trunk:
            self._plot_item(item, self.lines)
        # Add a bit of padding above & below the plot:
        vspace = (self.ymax - self.ymin) * 0.01
        self.ymin -= vspace
        self.ymax += vspace
        
        # Create the line collection:
        self.line_collection = mpl.collections.LineCollection(self.lines, linewidths = self.line_width)
        self.line_collection.set_color(color)
        
        # Now attach to the given matplotlib axes:
        self.axes = axes
        if axes:
            axes.set_xlim([self.xmin, self.xmax]) 
            axes.set_ylim([self.ymin, self.ymax])
            axes.set_xticks([])
            axes.set_xticklabels([])
            if self.line_width > 1:
                # Y values will not be correct, so hide them:
                axes.set_yticks([])
                axes.set_yticklabels([])
            axes.add_collection(self.line_collection)
    
    @property
    def xmin(self): return 0 - self.x_increment
    @property
    def xmax(self): return self._next_x + self.x_increment
    
    def item_at(self,x,y):
        """ Returns the item at the given point on the plot, or None """
        if x is None or y is None: # a convenience check that makes None an acceptable parameter
            return None
        for (xmin,xmax,ymin,ymax), item in self._item_map.iteritems():
            if x >= xmin and x <= xmax and y > ymin and y < ymax:
                return item
        return None
    
    def _plot_item(self, item, line_list):
        if item.parent:
            parent_y_level = item.parent.merge_level*self.flux_mult
        else:
            parent_y_level = self.ymin
            
        if type(item) == Leaf:
            # Draw the vertical line for this leaf:
            top_y = item.fmax*self.flux_mult
            if item.idx in self._x_map:
                # We've plotted this before, so we know
                # what x value to use
                x = self._x_map[item.idx]
            else: # This is an item never before plotted
                x = self._next_x
                self._next_x += self.x_increment
                xmin,xmax,ymin,ymax = (x-self.x_increment/2,x+self.x_increment/2, parent_y_level, top_y)
                # Add the item to our maps:
                self._item_map[(xmin,xmax,ymin,ymax)] = item
                self._x_map[item.idx] = x
            line_list.append([(x, top_y), (x, parent_y_level)])
            if top_y > self.ymax:
                self.ymax = top_y
            return x
        elif type(item) == Branch:
            # Plot child branches & leaves, recording the x values at which they are plotted:
            xvalues = [self._plot_item(child, line_list) for child in item.items]
            mean_x = int(sum(xvalues) / len(xvalues))
            y_level = item.merge_level * self.flux_mult
            # Add vertical line for this branch:
            line_list.append( [(mean_x, parent_y_level), (mean_x, y_level)] )
            # Add horizontal line for this branch:
            line_list.append( [(xvalues[0], y_level), (xvalues[-1], y_level)] )
            if item.idx not in self._x_map: # If this is the main plot:
                self._item_map[(xvalues[0], xvalues[-1], parent_y_level,y_level)] = item
                self._x_map[item.idx] = xvalues[0]
            return mean_x
    
    def create_highlighter(self, color='red', alpha=1):
        """
        Set up a plot for interactive highlighting.
        This will return a Highlighter object with a highlight(item)
        method; simply call that method to highlight any given item in the 
        color specified when the highlighter was created.
        """
        h = self._Highlighter(self, color, alpha, self.line_width+1)
        self.axes.add_collection(h.line_collection)
        return h
    
    
    class _Highlighter:
        def __init__(self, plot, color, alpha, line_width):
            """ Do not call this yourself, but rather use plot.create_highlighter() """
            self.item = None
            self.line_collection = mpl.collections.LineCollection([], linewidths = line_width)
            self.line_collection.set_color(color)
            self.line_collection.set_alpha(alpha)
            self.plot = plot
        def highlight(self, item):
            """
            Highlight the given item (a Branch or Leaf object)
            Item can be None, which will remove highlighting.
            Returns True if the highlight changed at all.
            """
            if self.item == item:
                return False # No change
            elif item is None:
                self.clear()
            else:
                self.item = item
                item_lines = []
                self.plot._plot_item(item, item_lines)
                self.line_collection.set_segments(item_lines)
            return True
        def clear(self):
            self.item = None
            self.line_collection.set_segments([])