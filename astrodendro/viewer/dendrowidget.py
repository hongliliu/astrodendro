'''
@author: Braden MacDonald
'''

import gtk
import matplotlib
import numpy as np
import astrodendro
from matplotlib.backends.backend_gtkagg import FigureCanvasGTKAgg, NavigationToolbar2GTKAgg

class DendrogramViewWidget(gtk.VBox):

    def __init__(self, cube_view, parent_window):
        """
        A Dendrogram Plot Widget. Designed as a companion to an 
        astrocube.cubeview.CubeViewWidget which should be passed as cube_view
        """
        gtk.VBox.__init__(self, False)

        # Main data structures:         
        self.cube = cube_view.cube # The cube that this widget's dendrogram represents
        self.data = cube_view.cube.data # Save a reference to the cube's original data, as we may replace it...
        self.dendrogram = None # The Dendrogram for the cube. To be built later

        # The plot:
        self.fig = matplotlib.figure.Figure()
        self.axes = self.fig.add_subplot(111)
        self.dendro_plot = None # Gets set to a DendrogramPlot object
        
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
        self.needs_redraw = False # Set this to True if you want the canvas to be repainted
        canvas.mpl_connect('button_press_event', self._figure_mousedown)
        canvas.mpl_connect('button_release_event', self._figure_mouseup)
        canvas.mpl_connect('motion_notify_event', self._figure_mousemoved)
        gtk.idle_add(DendrogramViewWidget._check_redraw, self) # we only want to re re-drawing when the GUI is idle, for maximum interactivity
        
        # Finally, put the default message on the status bar:
        self._plot_toolbar.update_mouseout_message()
    
    def _figure_mousedown(self, event):
        if event.xdata != None and event.ydata != None: # If we're in the canvas:
            self._is_mouse_down = True
    def _figure_mouseup(self, event):
        self._is_mouse_down = False
    def _figure_mousemoved(self, event):
        if self._is_mouse_down and (event.xdata != None and event.ydata != None): # If we're in the canvas:
            print("motion ({0},{1})".format(int(event.xdata),int(event.ydata)))
        # Note other mouse motion updates get processed below in _NavigationToolbar.mouse_move
        if (event.xdata != None and event.ydata != None and self.dendro_plot):
            item = self.dendro_plot.item_at(event.xdata, event.ydata)
            if item:
                self.dendro_plot.highlight(item=item)
                self.needs_redraw = True

    def _compute_btn_clicked(self, btn, event):
        if not self.dendrogram:
            self.dendrogram = astrodendro.Dendrogram(self.data, compute=False)
        
        min_flux = self._min_flux_widget.get_value()
        min_npix = self._min_npix_widget.get_value()
        min_delta= self._min_delta_widget.get_value()
        self.dendrogram.compute(minimum_flux=min_flux, minimum_npix=min_npix, minimum_delta=min_delta)
        self.axes.clear()
        self.dendro_plot = self.dendrogram.make_plot()
        self.dendro_plot.add_to_axes(self.axes)
        self.needs_redraw = True

    def _check_redraw(self):
        ''' Update this widget's display if needed. Called only when the main event loop is idle '''
        if self.needs_redraw:
            self.fig.canvas.draw()
            self.needs_redraw = False
        return True

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
    
            if event.inaxes and event.inaxes.get_navigate():
                # We are hovering over the plot, so display the data coordinates and real coordinates:
                x,y = event.xdata, event.ydata
                self.set_message(u"({x}, {y})".format(x=int(x), y=int(y)))
            else:
                self.update_mouseout_message()
        def update_mouseout_message(self):
            ''' Set the message shown when the user's cursor is not over the plot '''
            if not self.get_parent().dendrogram:
                self.set_message(u"No dendrogram.")
