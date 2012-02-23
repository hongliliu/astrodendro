import matplotlib as mpl
from .components import Leaf, Branch

class DendrogramPlot:
    """
    A simple virtual class used for setting up more specific types of plots.
    This class is only responsible for setting up matplotlib according to
    the wishes of (i.e. the canvas size specified by) a specific plot type.
    
    Specific plot types need to subclass this an implement:
    -> a constructor that plots the dendrogram and then calls __init__
       to prepare axes and highlighters for plotting, and create a
       matplotlib LineCollection to store the dendrogram plot lines
    -> item_at(x,y)
    -> _plot_item(item, line_collection, include_children)
    """

    def __init__(self, axes, line_width):
        """
        Set up the given matplotlib axes for plotting.
        This will create a line_collection object containing the lines and colors
        returned by _plot_trunk().
        
        Because this calls _plot_trunk(), you should do the initialization
        specific to each type of plot in the subclass before calling __init__()
        on this superclass.
        """
        # Create the line collection:
        self.line_width =  line_width
        lines, colors = self._plot_trunk()
        self.line_collection = mpl.collections.LineCollection(lines, colors=colors, linewidths = self.line_width)
        
        # Now attach to the given matplotlib axes:
        self.axes = axes
        if axes: # may be None if we are testing plot computation speed
            axes.set_xticks([])
            axes.set_xticklabels([])
            if self.line_width > 1:
                # Y values will not be correct, so hide them:
                axes.set_yticks([])
                axes.set_yticklabels([])
            axes.add_collection(self.line_collection)
            
            # A list of highlighters created on this plot:
            self.highlighters_active = []
            # Methods can be registered that will get called whenever a highlight changes:
            self._highlight_change_notify = []

            axes.set_xmargin(0.05)
            axes.set_ymargin(0.05)
            axes.autoscale()

    ########## Subclasses need to implement these: #################
    def _plot_trunk(self):
        """ Plot the whole dendrogram. Returns (lines, colors). """
        raise Exception("This type of plot has not implemented _plot_trunk() !")
    def _plot_item(self, item):
        """ Plot an item and its children. Returns (lines, colors). """
        raise Exception("This type of plot has not implemented _plot_item() !")
    def item_at(self,x,y):
        """ Returns the item at the given point on the plot, or None """
        raise Exception("This type of plot has not implemented item_at() !")
    ################################################################
    
    def on_highlight_change(self, func):
        """ Register a function to be called whenever any highlight changes """
        if func not in self._highlight_change_notify:
            self._highlight_change_notify.append(func)
    
    def create_highlighter(self, color='red', alpha=1, line_width_extra=1):
        """
        Set up a plot for interactive highlighting.
        This will return a Highlighter object with a highlight(item)
        method; simply call that method to highlight any given item in the 
        color specified when the highlighter was created.
        """
        h = self._Highlighter(self, color, alpha, self.line_width + line_width_extra)
        self.highlighters_active.append(h)
        self.axes.add_collection(h.line_collection)
        return h

    class _Highlighter:
        def __init__(self, plot, color, alpha, line_width):
            """ Do not call this yourself, but rather use plot.create_highlighter() """
            self.item = None
            self.line_collection = mpl.collections.LineCollection([], linewidths=line_width)
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
                item_lines, _ = self.plot._plot_item(item)
                self.line_collection.set_segments(item_lines)
                for func in self.plot._highlight_change_notify:
                    func()
            return True

        def clear(self):
            self.item = None
            self.line_collection.set_segments([])
            for func in self.plot._highlight_change_notify:
                func()

class RecursiveSortPlot(DendrogramPlot):
    """
    A plot that sorts the dendrogram according to the max. flux value found
    among that branch and its children
    """
    def __init__(self, dendrogram, color_lambda, axes, line_width, spacing=5):
        
        self._dendrogram = dendrogram
        self._color_lambda = color_lambda
        
        self.spacing = spacing
        self.x_increment = line_width + spacing
        self._flux_mult = line_width # All y values (flux values) should be multiplied by this factor 
        self._next_x = 0
        
        # item_x_map: A map of the positions of each item
        # recursing through all items in the dendrogram.
        # Keys are item.idx
        # Values are (type, xmin, xmax, ymin, ymax, color)
        # For Leaf items, xmin=xmax.
        self._item_rect_map = {} # map of x index
        
        # Build our map of the positions of each item:
        self._build_rect_map(dendrogram.trunk, sort_by=lambda item: item.get_peak_recursive()[2])
        
        DendrogramPlot.__init__(self, axes=axes, line_width=line_width)

    def _build_rect_map(self, items, sort_by):
        items_sorted = sorted(items, key=sort_by)
        for item in items_sorted:
            if item.parent:
                ymin = item.parent.merge_level * self._flux_mult
            else:
                ymin = 0
                
            color = self._color_lambda(item)
            
            if type(item) == Leaf:
                ymax = item.fmax * self._flux_mult
                self._item_rect_map[item.idx] = (Leaf, self._next_x, self._next_x, ymin, ymax, color)
                self._next_x += self.x_increment
            elif type(item) == Branch:
                xmin, xmax = self._build_rect_map(item.items, sort_by) # This is where the recursion happens
                ymax = item.merge_level*self._flux_mult
                self._item_rect_map[item.idx] = (Branch, xmin, xmax, ymin, ymax, color)
        # Return the position of the vertical line representing the first and last item:
        _,first_xmin,first_xmax,_,_,_ = self._item_rect_map[items_sorted[0].idx]
        _,last_xmin,last_xmax,_,_,_ = self._item_rect_map[items_sorted[-1].idx]
        return (first_xmin+first_xmax)//2, (last_xmin+last_xmax)//2
        

    def _plot_trunk(self):
        # Now do the plot:
        return self._plot_item_rects(self._item_rect_map.itervalues())
        
    def _plot_item(self, item):
        """ Plot an item and its children. """
        items = [item]
        if type(item) is Branch:
            items.extend(item.descendants)
        return self._plot_item_rects([self._item_rect_map[i.idx] for i in items])
    
    def _plot_item_rects(self, item_rects):
        """
        Plot the list of item rects, which can be unordered.
        ** Child item rects must be included in the list or they won't get plotted. **
        """
        lines = []
        colors = []
        for (itemtype, xmin, xmax, ymin, ymax, color) in item_rects:
            if itemtype == Leaf:
                # Draw a vertical line:
                lines.append([(xmin,ymin), (xmax, ymax)]) # Note for leaves, xmin = xmax
                colors.append(color)
            elif itemtype == Branch:
                # Draw a vertical line:
                xmean = (xmin+xmax)//2
                lines.append([(xmean,ymin), (xmean, ymax)])
                # Draw a horizontal line:
                lines.append([(xmin, ymax), (xmax, ymax)])
                colors.extend([color, color]) # Add two colors, one for horiz. line, one for vertical
        return lines, colors

    def item_at(self,x,y):
        """ Returns the item at the given point on the plot, or None """
        if x is None or y is None: # a convenience check that makes None an acceptable parameter
            return None
        s = self.spacing//2
        for idx, (t, xmin,xmax,ymin,ymax,_) in self._item_rect_map.iteritems():
            # If it's a branch:
            if t is Branch and x >= xmin and x <= xmax and y > ymin and y < ymax:
                return self._dendrogram.items_dict[idx]
            elif x >= (xmin-s) and x <= (xmax+s) and y > ymin and y < ymax:
                return self._dendrogram.items_dict[idx]
        return None

class SpatialCoordPlot(DendrogramPlot):
    """
    Create a plot that places leaves at the position of the X coordinate of the
    pixel with maximum flux value in the leaf. (Or use Y coordinate, or Z, ...)
    
    coord_axis: set to 0 for X, 1 for Y, etc.
    """

    def __init__(self, dendrogram, color_lambda, axes, line_width, coord_index = 0):
        self.dendrogram = dendrogram
        self.coord_index = coord_index
        self._color_lambda = color_lambda
        DendrogramPlot.__init__(self, axes=axes, line_width=line_width)

    def _plot_trunk(self):
        (lines, colors,) = ([], [])
        for item in self.dendrogram.trunk:
            self._plot_item_recursive(item, lines, colors)
        return (lines, colors)

    def _plot_item(self, item):
        (lines, colors,) = ([], [])
        self._plot_item_recursive(item, lines, colors)
        return (lines, colors)

    def _plot_item_recursive(self, item, lines, colors):
        """
        Plots an item, and all children recursively.
        Adds the item's lines and colors to the lists given.
        Returns the mean x value of the item plotted.
        """
        color = self._color_lambda(item)
        ymin = item.parent.merge_level if item.parent else 0
        if type(item) is Leaf:
            (peak_coords, ymax,) = item.peak
            x = peak_coords[self.coord_index]
            lines.append([(x, ymin), (x, ymax)])
            colors.extend([color])
            return x
        else:
            y = item.merge_level
            x_values = [ self._plot_item_recursive(i, lines, colors) for i in item.items ]
            (xmin, xmax,) = (min(x_values), max(x_values))
            lines.append([(xmin, y), (xmax, y)])
            xmean = (xmin + xmax) // 2
            lines.append([(xmean, ymin), (xmean, y)])
            colors.extend([color, color])
            return xmean

    def item_at(self, x, y):
        """ Returns the item at the given point on the plot, or None """
        raise Exception("This type of plot has not implemented item_at() !")

class FuturePlot(DendrogramPlot):
    """ A template for writing new styles of dendrogram plots """
    def __init__(self, axes, line_width, other_args_to_add):
        # Set up plot here
        
        # The following will in turn set up the plotting axes, then call _plot_trunk
        DendrogramPlot.__init__(self, axes=axes, line_width=line_width)
    def _plot_trunk(self):
        """ Plot the whole dendrogram. Returns (lines, colors). """
        raise Exception("This type of plot has not implemented _plot_trunk() !")
    def _plot_item(self, item):
        """ Plot an item and its children. Returns (lines, colors). """
        raise Exception("This type of plot has not implemented _plot_item() !")
    def item_at(self,x,y):
        """ Returns the item at the given point on the plot, or None """
        raise Exception("This type of plot has not implemented item_at() !")
