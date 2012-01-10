'''
@author: Braden MacDonald
'''

import gtk
import matplotlib
import numpy as np
import astrodendro
from astrodendro.components import Branch, Leaf
from matplotlib.backends.backend_gtkagg import FigureCanvasGTKAgg, NavigationToolbar2GTKAgg



class DendrogramViewWidget(gtk.VBox):

    def __init__(self, cube, parent_window):
        """
        A Dendrogram Plot Widget. Designed as a companion to an 
        astrocube.cubeview.CubeViewWidget which should be passed as cube_view
        """
        gtk.VBox.__init__(self, False)

        # Main data structures:         
        self.cube = cube # The cube that this widget's dendrogram represents
        self.data = cube.data # Save a reference to the cube's original data, as we may replace it...
        self.dendrogram = None # The Dendrogram for the cube. To be built later

        # The plot:
        self.fig = matplotlib.figure.Figure()
        self.axes = self.fig.add_subplot(111)
        self.dendro_plot = None # Gets set to a DendrogramPlot object
        self.highlighter_clicked = None # a Highlighter object to color the clicked item
        self.highlighter_hover = None # a Highlighter object to color an item on mouseover
        
        # UI structure:
        canvas = FigureCanvasGTKAgg(self.fig)  # a gtk.DrawingArea
        self.pack_start(canvas)
        self.pack_start(gtk.HSeparator(), False,False)
        
        # # The matplotlib plot toolbar:
        self._plot_toolbar = self._NavigationToolbar(self.cube, canvas, parent_window) # parent_window is needed for the "Save image..." file chooser dialog
        self._plot_toolbar.remove(self._plot_toolbar.get_nth_item(6)) # Remove the "Configure subplots" button which causes rendering issues if used
        self.pack_start(self._plot_toolbar, False, False)
        
        self.pack_start(gtk.HSeparator(), False,False)
        
        # # the dendrogram parameters toolbar:
        self._dendro_toolbar = gtk.HBox()
        
        # # # Minimum flux:
        self._dendro_toolbar.pack_start(gtk.Label("Min flux: "))
        self._min_flux_widget = gtk.SpinButton(digits=1)
        self._min_flux_widget.set_range(np.nanmin(self.data)-0.1,np.nanmax(self.data))
        self._min_flux_widget.set_increments(0.1,1)
        self._min_flux_widget.set_value(np.nanmin(self.data))
        self._dendro_toolbar.pack_start(self._min_flux_widget)
        
        # # # Minimum npix:
        self._dendro_toolbar.pack_start(gtk.Label("Min # pix: "))
        self._min_npix_widget = gtk.SpinButton()
        self._min_npix_widget.set_range(0,999)
        self._min_npix_widget.set_increments(1,5)
        self._min_npix_widget.set_value(0)
        self._dendro_toolbar.pack_start(self._min_npix_widget)
        
        # # # Minimum delta:
        self._dendro_toolbar.pack_start(gtk.Label("Min delta: "))
        self._min_delta_widget = gtk.SpinButton(digits=2)
        self._min_delta_widget.set_range(0,np.nanmax(self.data))
        self._min_delta_widget.set_increments(0.05,0.3)
        self._min_delta_widget.set_value(0)
        self._dendro_toolbar.pack_start(self._min_delta_widget)
        
        # # # Compute button:
        self._compute_button = gtk.Button("Compute")
        self._compute_button.connect("button-press-event", self._compute_btn_clicked)
        self._dendro_toolbar.pack_start(self._compute_button)
        
        self.pack_start(self._dendro_toolbar, False, False)
        
        
        # Set up event handling:
        self._is_mouse_down = False # is the mouse button currently pressed?
        self._redraw_all = False # Set this to True to trigger a complete re-draw of the canvas when idle
        self._redraw_highlights = False # Set this True to  re-draw only the highlighted (clicked or hovered) leaves
        canvas.mpl_connect('button_press_event', self._figure_mousedown)
        canvas.mpl_connect('button_release_event', self._figure_mouseup)
        canvas.mpl_connect('motion_notify_event', self._figure_mousemoved)
        canvas.mpl_connect('resize_event', self._figure_resized)
        gtk.idle_add(DendrogramViewWidget._check_redraw, self) # we only want to re re-drawing when the GUI is idle, for maximum interactivity
        
        self._click_notify = []
        self._compute_notify = []
    
    def make_dendrogram(self, min_flux, min_npix, min_delta):
        " Call this to manually set the parameters and create a new dendrogram"
        self._min_flux_widget.set_value(min_flux)
        self._min_npix_widget.set_value(min_npix)
        self._min_delta_widget.set_value(min_delta)
        self._compute_btn_clicked(self._compute_button, None)
    
    def on_click(self, func):
        """
        Register a function to be called when the user clicks on an item.
        Passes the item (or None) to the event handler.
        """
        if not func in self._click_notify:
            self._click_notify.append(func) 
    
    def on_compute(self, func):
        """
        Register a function to be called when the user [re]computes a dendrogram
        Passes the new dendrogram object to the event handler
        """
        if not func in self._compute_notify:
            self._compute_notify.append(func)
    
    def _figure_mousedown(self, event):
        if event.xdata != None and event.ydata != None: # If we're in the canvas:
            self._is_mouse_down = True
    def _figure_mouseup(self, event):
        self._is_mouse_down = False
        if self.highlighter_clicked:
            item = self.dendro_plot.item_at(event.xdata, event.ydata)
            if self.highlighter_clicked.highlight(item):
                self._redraw_highlights = True
            for handler in self._click_notify:
                handler(item)
    def _figure_mousemoved(self, event):
        if self.highlighter_hover:
            item = self.dendro_plot.item_at(event.xdata, event.ydata)
            if self.highlighter_hover.highlight(item): # return true if changed:
                self._redraw_highlights = True
        # Note other mouse motion updates get processed below in _NavigationToolbar.mouse_move

    def _figure_resized(self, event):
        self._redraw_all = True
    
    def _highlights_changed(self):
        self._redraw_highlights = True

    def _compute_btn_clicked(self, btn, event):
        if not self.dendrogram:
            self.dendrogram = astrodendro.Dendrogram(self.data, compute=False)
        
        min_flux = self._min_flux_widget.get_value()
        min_npix = self._min_npix_widget.get_value()
        min_delta= self._min_delta_widget.get_value()
        self.dendrogram.compute(minimum_flux=min_flux, minimum_npix=min_npix, minimum_delta=min_delta)
        self.axes.clear()
        self.dendro_plot = self.dendrogram.make_plot(self.axes)
        self.highlighter_clicked = self.dendro_plot.create_highlighter('red', alpha=1)
        self.highlighter_hover = self.dendro_plot.create_highlighter('green', alpha=0.7)
        self.dendro_plot.on_highlight_change(self._highlights_changed)
        # The first time we render this diagram, we call the draw() method:
        self.fig.canvas.draw()
        # From now on we use our more efficient rendered triggered by _redraw_all: 
        self._redraw_all = True
        for handler in self._compute_notify:
            handler(self.dendrogram)

    def _check_redraw(self):
        '''
        Update this widget's plot if needed.
        Called only when the main event loop is idle.
        Does everything in a very controlled manner for max. performance & 
        interactivity.
        '''
        if self._redraw_all:
            # Render just the dendrogram:
            if self.dendro_plot:
                self.fig.canvas.renderer.clear()
                self.axes.draw_artist(self.dendro_plot.line_collection)
                # Now cache the result:
                self._mpl_plot_bitmap_cache = self.fig.canvas.copy_from_bbox(self.axes.bbox)
                self._redraw_highlights = True
            self._redraw_all = False
        if self._redraw_highlights:
            # Restore the previously cached axes lines and dendrogram plot:
            self.fig.canvas.restore_region(self._mpl_plot_bitmap_cache)
            # Render the highlights:
            if self.dendro_plot:
                for highlighter in self.dendro_plot.highlighters_active:
                    self.axes.draw_artist(highlighter.line_collection)
                    # This is the big bottleneck slowing interactivity
            # Redraw the axes:
            self.fig.canvas.blit(self.axes.bbox)
            self._redraw_highlights = False
        return True
    
    def set_clicked_item_by_coords(self, coords):
        " Highlight whatever item is at the specified coordinates "
        if self.dendrogram and self.highlighter_clicked:
            item = self.dendrogram.item_at(coords)
            if self.highlighter_clicked.highlight(item):
                self._redraw_highlights = True
            return item
        return None

    class _NavigationToolbar(NavigationToolbar2GTKAgg):
        def __init__(self, cube, canvas, parent_window):
            self.cube = cube
            NavigationToolbar2GTKAgg.__init__(self, canvas, parent_window)
        def mouse_move(self, event):
            #print 'mouse_move', event.button
    
            cursors = matplotlib.backend_bases.cursors
            if not event.inaxes or not self._active:
                if self._lastCursor != cursors.POINTER:
                    self.set_cursor(cursors.POINTER)
                    self._lastCursor = cursors.POINTER
            else:
                if self._active=='ZOOM':
                    if self._lastCursor != cursors.SELECT_REGION:
                        self.set_cursor(cursors.SELECT_REGION)
                        self._lastCursor = cursors.SELECT_REGION
                elif (self._active=='PAN' and self._lastCursor != cursors.MOVE):
                    self.set_cursor(cursors.MOVE)
                    self._lastCursor = cursors.MOVE
    
            dendro = self.get_parent().dendro_plot
            if dendro:
                item = dendro.item_at(event.xdata, event.ydata)
                if type(item) == Branch:
                    msg = "Branch with {0} children; has {1:,} ({2:,}) pixels".format(len(item.items), item.npix_self, item.npix)
                elif type(item) == Leaf:
                    msg = "Leaf with {0} pixels".format(item.npix)
                else:
                    msg = ""
            else:
                msg = "No dendrogram computed yet."
            
            self.set_message(msg)
        def draw(self):
            """
            This gets called on zoom or pan. It updates the locators then tells
            the canvas to draw.
            We override the default implementation and set our 
            DendrogramViewWidget redraw flag instead of calling the slow canvas
            draw() method. 
            """
            for a in self.canvas.figure.get_axes():
                xaxis = getattr(a, 'xaxis', None)
                yaxis = getattr(a, 'yaxis', None)
                locators = []
                if xaxis is not None:
                    locators.append(xaxis.get_major_locator())
                    locators.append(xaxis.get_minor_locator())
                if yaxis is not None:
                    locators.append(yaxis.get_major_locator())
                    locators.append(yaxis.get_minor_locator())
    
                for loc in locators:
                    loc.refresh()
            self.get_parent()._redraw_all = True