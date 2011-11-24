import string
import numpy as np

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

    def plot_dendrogram(self, plot, parent_y_level):
        # Draw the vertical line for this leaf:
        x = plot.cur_x
        plot.add_line( (x, self.fmax*plot.flux_mult), (x, parent_y_level) )
        plot.cur_x += plot.x_increment
        return x

    def to_newick(self):
        return "%i:%.3f" % (self.idx, self.height)

    def get_peak(self):
        imax = np.argmax(self.f)
        return self.coords[imax], self.f[imax]


class Branch(Leaf):

    def __init__(self, items, coord, f, idx=None):
        self.merge_level = f # Record the exact flux level that triggered creation of this branch
        self.items = items
        for item in items:
            item.parent = self
        Leaf.__init__(self, coord, f, idx=idx)

    @property
    def npix(self):
        return len(self.f) + self.npix_children
    @property
    def npix_children(self):
        return np.sum([item.npix for item in self.items])
    
    @property
    def f_sum(self):
        return np.sum(self.f) + self.f_sum_children
    @property
    def f_sum_children(self):
        return np.sum([item.f_sum for item in self.items])


    def add_footprint(self, image, level, recursive=True):
        if recursive:
            for item in self.items:
                item.add_footprint(image, level + 1)
        Leaf.add_footprint(self, image, level)

    def plot_dendrogram(self, plot, parent_y_level):
        y_level = self.merge_level * plot.flux_mult
        # Plot child branches & leaves, recording the x values at which they are plotted:
        xvalues = [item.plot_dendrogram(plot, y_level) for item in self.items]
        mean_x = int(sum(xvalues) / len(xvalues))
        # Add vertical line for this branch:
        plot.add_line( (mean_x, parent_y_level), (mean_x, y_level) )
        # Add horizontal line for this branch:
        plot.add_line( (xvalues[0], y_level), (xvalues[-1], y_level) )
        return mean_x

    def to_newick(self):
        newick_items = []
        for item in self.items:
            newick_items.append(item.to_newick())
        return "(%s)%s:%.3f" % (string.join(newick_items, ','), self.idx, self.height)

    def get_leaves(self):
        leaves = []
        for item in self.items:
            if type(item) == Leaf:
                leaves.append(item)
            else:
                leaves += item.get_leaves()
        return leaves


class Trunk(list):

    def plot_dendrogram(self, line_width, spacing):
        # Find the minimum flux among all root branches:
        min_f = np.min([item.fmin for item in self])
        # Set up variables needed for plotting:
        plot = DendrogramPlotInfo(line_width, spacing, min_f)
        # recursively generate the necessary lines:
        for item in self:
            item.plot_dendrogram(plot, plot.ymin)
        # Add a bit of padding above & below the plot:
        plot_vspace = (plot.ymax - plot.ymin) * 0.01
        plot.ymin -= plot_vspace
        plot.ymax += plot_vspace
        return plot

    def to_newick(self):
        newick_items = []
        for item in self:
            newick_items.append(item.to_newick())
        return "(%s);" % string.join(newick_items, ',')

    def get_leaves(self):
        leaves = []
        for item in self:
            if type(item) == Leaf:
                leaves.append(item)
            else:
                leaves += item.get_leaves()
        return leaves

class DendrogramPlotInfo():
    def __init__(self, line_width, spacing, min_flux):
        self.line_width = line_width
        self.spacing = spacing
        self.x_increment = line_width + spacing
        self.flux_mult = line_width # All y values (flux values) should be multiplied by this factor 
        self.cur_x = 0
        self.lines = [] # list of lists of (x,y) tuples
        
        # the xmin, xmax, ymin, ymax properties give the size of the canvas onto which the lines are drawn
        self.ymin = min_flux * self.flux_mult
        self.ymax = self.ymin
    @property
    def xmin(self): return 0 - self.x_increment
    @property
    def xmax(self): return self.cur_x + self.x_increment
    def add_line(self, point1, point2):
        self.lines.append([point1, point2])
        if point1[1] > self.ymax: self.ymax = point1[1]
        if point2[1] > self.ymax: self.ymax = point2[1]
