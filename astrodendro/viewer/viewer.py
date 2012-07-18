#!/usr/bin/env python
# coding: utf-8
'''
Standalone data cube dendrogram viewer.
Uses GTK as a widget library

@author: Braden MacDonald
'''
import sys
import gtk

from astrocube import DataCube
from astrocube.cubeview import CubeViewWidget
from astrodendro.viewer.dendrowidget import DendrogramViewWidget
from astrodendro.viewer.ipython_view import IPythonView
from matplotlib.backends.backend_pdf import PdfPages

import numpy as np


class DendroViewer:
    
    def __init__(self, filename):
        self.cube = DataCube(filename, calc_noise_dev=False)
        
        win = gtk.Window()
        win.connect("destroy", lambda x: gtk.main_quit())
        win.set_default_size(800,600)
        win.set_title("{ln} spectral line map of {o}".format(o=self.cube.object_name,ln=self.cube.line_name))
        
        # Create our three main widgets:
        
        # A CubeViewWidget and a connected DendrogramViewWidget in the top row:
        self.cube_view = CubeViewWidget(self.cube, win)
        self.dendro_view = DendrogramViewWidget(self.cube, win)
        first_row = gtk.HPaned()
        first_row.pack1(self.cube_view, resize=True, shrink=True)
        first_row.pack2(self.dendro_view, resize=True, shrink=True)
        
        # and an IPython command widget in the bottom row:
        
        ipython_namespace = {'cube': self.cube,
                             'cube_view': self.cube_view,
                             'dendro_view': self.dendro_view,
                             'create_highlighter': self.create_highlighter,
                             'set_color_map': self.set_color_map,
                             'set_data': self.set_data,
                             'make_dendrogram': self.dendro_view.make_dendrogram,
                             'export_pdf': self.export_pdf,
                             'export_png': self.export_png }
        
        self.ipython_widget = IPythonView()
        for key,val in ipython_namespace.items():
            self.ipython_widget.IP.user_ns[key] = val 
        # An IPython widget does not support scrolling by itself...
        self.ipython_widget.set_wrap_mode(gtk.WRAP_CHAR)
        ipython_scroller = gtk.ScrolledWindow() # so we wrap it in a scrolled window
        ipython_scroller.set_size_request(-1,75)
        ipython_scroller.set_policy(gtk.POLICY_NEVER, gtk.POLICY_ALWAYS)
        ipython_scroller.add(self.ipython_widget)
       
        
        main_pane = gtk.VPaned()
        main_pane.pack1(first_row, resize=True, shrink=True)
        main_pane.pack2(ipython_scroller, resize=True, shrink=False)
        main_pane.set_position(9000) # Expand the figure to be as big as possible
        
        win.add(main_pane)
        self.win = win
        
        # Set up events:
        self.cube_view.on_click(self.cube_clicked)
        self.cube_highlight = self.cube_view.create_highlighter()
        self.dendro_view.on_click(self.dendro_clicked)
        self.dendro_view.on_compute(self.dendro_computed)
        self.highlighted_item = None
        
    def run(self):
        self.win.show_all()
        gtk.main()
    
    def dendro_computed(self, dendrogram):
        self.ipython_widget.IP.user_ns['dendro'] = dendrogram
    
    def cube_clicked(self, coords, flux):
        item_clicked = self.dendro_view.set_clicked_item_by_coords(coords)
        if (item_clicked):
            self.highlight_item_in_cube(item_clicked)
    
    def dendro_clicked(self, item):
        if item:
            if item != self.highlighted_item:
                self.highlighted_item = item
                self.highlight_item_in_cube(item)
            else:
                # User has clicked on the item already highlighted.
                # Mark the point of greatest flux in the view:
                _,max_coords,_ = item.get_peak_recursive()
                self.cube_view.x, self.cube_view.y, self.cube_view.z = max_coords    
        else:
            self.highlighted_item = None
            self.cube_highlight.clear()
    
    def highlight_item_in_cube(self, item):
        mapdata = np.zeros(self.cube.data.shape)
        item.add_footprint(mapdata, 1)
        mapdata[mapdata>1] = 0.75 # Set the child items to be semi-transparent
        self.cube_highlight.highlight(mapdata)
    
    ### Convenience methods for the user to enter via IPython shell:
    
    def set_data(self, data):
        """
        Call this method to change the data shown in the cube view and used
        to compute the dendrogram.
        """
        self.cube.data = data
        self.cube_view.update()
        self.dendro_view.data = data # This will also clear the dendrogram
    
    def set_color_map(self, cmap = CubeViewWidget.default_cmap):
        """
        Sets the color map used by the cube view (left hand side).
        Use "bone" (greyscale), "spectral", or any other matplotlib color map.
        Pass no argument for the default color map.
        """
        self.cube_view.imgplot.set_cmap(cmap)
        self.cube_view.axes.figure.canvas.draw()
    
    def export_pdf(self, filename, title="default", subtitle="default", include_cube_slice=True):
        """
        Create a standardized one or two -page PDF of the current figures
        Filename is e.g. "export1.pdf"
        Title gets put on every figure in the PDF
        Set include_cube_slice to False to export only the dendrogram
        If include_cube_slide is true, the cube view will be exported,
        at the Z value currently visible in the viewer.
        """
        pp = PdfPages(filename)
        
        # Figure out the titles and subtitles
        cube_subtitle = None # This gets prepended to the subitle of the cube only, not the dendrogram
        if title is "default":
            # Default title:
            title = "Dendrogram from {obj} {line}".format(obj=self.cube.object_name, line=self.cube.line_name)
        if subtitle is "default":
            subtitle = "Dendrogram parameters: min. flux {f}, min. npix {n:.0f}, min. delta {d}".format(
                       f=self.dendro_view.dendrogram.min_flux,
                       n=self.dendro_view.dendrogram.min_npix,
                       d=self.dendro_view.dendrogram.min_delta,
                       )
            if include_cube_slice:
                if self.cube.has_coords:
                    cube_subtitle = "@ {v:.2f} km/s (z={z}). {rest}".format(
                        v=self.cube.velocity_at(self.cube_view.z),
                        z=self.cube_view.z,
                        rest=subtitle,
                        )
                else:
                    cube_subtitle = "@ z={z}. {rest}".format(z=self.cube_view.z, rest=subtitle)
        
        export_size = (10,7.5) # in inches
        
        if include_cube_slice:
            cfig = self.cube_view.fig
            self.cube_view._check_redraw()
            old_size_cfig = (cfig.get_figwidth(), cfig.get_figheight())
            cfig.set_size_inches(*export_size)
            if title: ctitle = cfig.suptitle(title, weight='heavy', size='large')
            if cube_subtitle: ctitle2 = cfig.suptitle(cube_subtitle, y=0.95)
            pp.savefig(cfig)
            # Reset the figure to how it was before:
            if title: ctitle.set_visible(False) # Closest we can come to deleting these titles
            if cube_subtitle: ctitle2.set_visible(False)
            cfig.set_figwidth(old_size_cfig[0])
            cfig.set_figheight(old_size_cfig[1])
        
        # Now export the dendrogram figure:
        dfig = self.dendro_view.fig
        #old_size_dfig = dfig.get_size_inches()
        old_size_dfig = (dfig.get_figwidth(), dfig.get_figheight())
        dfig.set_size_inches(*export_size)
        if title: dtitle = dfig.suptitle(title, weight='heavy', size='large')
        if subtitle: dtitle2 = dfig.suptitle(subtitle, y=0.95)
        pp.savefig(self.dendro_view.fig)
        # Reset the figure to how it was before:
        if title: dtitle.set_visible(False)
        if subtitle: dtitle2.set_visible(False) 
        #dfig.set_size_inches(old_size_dfig)
        dfig.set_figwidth(old_size_dfig[0])
        dfig.set_figheight(old_size_dfig[1])
        
        # Finalize the PDF:
        if title: 
            pp.infodict()['Title'] = title
        pp.close()
        pp = None
        
        # Now the draw() method must be called, or some of the dendrogram 
        # plot artists will cache the PDF renderer and cause bugs, due to
        # the dendro_view's special handling of rendering
        self.dendro_view.fig.canvas.draw()
    
    def export_png(self, filename_cube, filename_dendro):
        from matplotlib.backends.backend_agg import FigureCanvasAgg
        if filename_cube:
            #self.cube_view.fig.canvas.print_png(filename_cube)
            canvas = self.cube_view.fig.canvas
            agg = canvas.switch_backends(FigureCanvasAgg)
            agg.print_png(filename_cube)
            agg = None
            self.cube_view.fig.set_canvas(canvas)
        if filename_dendro:
            #self.dendro_view.fig.canvas.print_png(filename_dendro)
            canvas = self.dendro_view.fig.canvas
            agg = canvas.switch_backends(FigureCanvasAgg)
            agg.print_png(filename_dendro)
            agg = None
            self.dendro_view.fig.set_canvas(canvas)
            # Now the draw() method must be called, or some of the dendrogram 
            # plot artists will cache the PDF renderer and cause bugs, due to
            # the dendro_view's special handling of rendering
            self.dendro_view.fig.canvas.draw()
    
    def create_highlighter(self, color):
        return self._CombinedHighlighter(self.cube_view, self.dendro_view, color)
    
    class _CombinedHighlighter:
        """
        A highlighter that highlights the data cube and the dendrogram
        Do not initialize directly, but get from create_highlighter()
        """
        def __init__(self, cube_view, dendro_view, color):
            self.color = color
            self.cube_view = cube_view
            self.dendro_view = dendro_view
            self._highlighter_cube = cube_view.create_highlighter(color)
            self._highlighter_dend = dendro_view.dendro_plot.create_highlighter(color) if dendro_view.dendro_plot else None
        def highlight_coords(self, coords):
            if self.dendro_view.dendrogram:
                self.highlight_item(self.dendro_view.dendrogram.item_at(coords))
            else:
                print("Cannot highlight until a dendrogram is plotted")
        def highlight_item(self, item):
            if not item:
                self._highlighter_cube.clear()
                if self._highlighter_dend:
                    self._highlighter_dend.clear()
                return
            mapdata = np.zeros(self.cube_view.cube.data.shape)
            item.add_footprint(mapdata, 1)
            mapdata[mapdata>1] = 0.75 # Set the child items to be semi-transparent
            if not self._highlighter_dend:
                # There was no dendrogram plot when this combined highlighter was created
                if self.dendro_view.dendro_plot:
                    # But one exists now:
                    self._highlighter_dend = self.dendro_view.dendro_plot.create_highlighter(self.color)
                else:
                    # And there still isn't, so we're done
                    return
            elif self._highlighter_dend.plot != self.dendro_view.dendro_plot:
                # The dendrogram has changed (been re-plotted), so re-create the highlighter:
                self._highlighter_dend.clear()
                self._highlighter_dend = self.dendro_view.dendro_plot.create_highlighter(self.color)
            # Now highlight on the dendrogram plot:
            self._highlighter_cube.highlight(mapdata)
            self._highlighter_dend.highlight(item)    

