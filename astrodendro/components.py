import string
import numpy as np


class Leaf(object):

    def __init__(self, x, y, z, f, idx=None):
        self.x = np.array([x], dtype=int)
        self.y = np.array([y], dtype=int)
        self.z = np.array([z], dtype=int)
        self.f = np.array([f], dtype=float)
        self.fmin, self.fmax = f, f
        self.idx = idx
        self.parent = None

    def __getattr__(self, attribute):
        if attribute == 'npix':
            return len(self.x)
        else:
            raise Exception("Attribute not found: %s" % attribute)

    def add_point(self, x, y, z, f):
        "Add point to current leaf"
        self.x = np.hstack([self.x, x])
        self.y = np.hstack([self.y, y])
        self.z = np.hstack([self.z, z])
        self.f = np.hstack([self.f, f])
        self.fmin, self.fmax = min(f, self.fmin), max(f, self.fmax)

    def merge(self, leaf):
        self.x = np.hstack([self.x, leaf.x])
        self.y = np.hstack([self.y, leaf.y])
        self.z = np.hstack([self.z, leaf.z])
        self.f = np.hstack([self.f, leaf.f])
        self.fmin, self.fmax = min(np.min(leaf.f), self.fmin), max(np.max(leaf.f), self.fmax)

    def add_footprint(self, image, level):
        "Fill in a map which shows the depth of the tree"
        image[self.z, self.y, self.x] = level
        return image

    def plot_dendrogram(self, plot, parent_y_level):
        # Draw the vertical line for this leaf:
        x = plot.cur_x
        plot.add_line( (x, self.fmax*plot.flux_mult), (x, parent_y_level) )
        plot.cur_x += plot.x_increment
        return x

    def to_newick(self):
        return "%i:%.3f" % (self.idx, self.fmax - self.fmin)

    def get_peak(self):
        imax = np.argmax(self.f)
        return self.x[imax], self.y[imax], self.z[imax], self.f[imax]


class Branch(Leaf):

    def __init__(self, items, x, y, z, f, idx=None):
        self.items = items
        for item in items:
            item.parent = self
        Leaf.__init__(self, x, y, z, f, idx=idx)

    def __getattr__(self, attribute):
        if attribute == 'npix':
            npix = len(self.x)
            for item in self.items:
                npix += item.npix
            return npix
        else:
            raise AttributeError("Attribute not found: %s" % attribute)

    def add_footprint(self, image, level, recursive=True):
        if recursive:
            for item in self.items:
                image = item.add_footprint(image, level + 1)
        return Leaf.add_footprint(self, image, level)

    def plot_dendrogram(self, plot, parent_y_level):
        y_level = self.fmin * plot.flux_mult
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
        return "(%s)%s:%.3f" % (string.join(newick_items, ','), self.idx, self.fmax - self.fmin)

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
