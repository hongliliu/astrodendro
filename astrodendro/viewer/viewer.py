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
    
        ipython_widget = IPythonView() # An IPython terminal. Does not support scrolling by itself...
        ipython_widget.set_wrap_mode(gtk.WRAP_CHAR)
        ipython_scroller = gtk.ScrolledWindow() # so we wrap it in a scrolled window
        ipython_scroller.set_size_request(-1,75)
        ipython_scroller.set_policy(gtk.POLICY_NEVER, gtk.POLICY_ALWAYS)
        ipython_scroller.add(ipython_widget)
       
        
        main_pane = gtk.VPaned()
        main_pane.pack1(first_row, resize=True, shrink=True)
        main_pane.pack2(ipython_scroller, resize=True, shrink=False)
        main_pane.set_position(9000) # Expand the figure to be as big as possible
        
        win.add(main_pane)
        self.win = win
        
        # Set up events:
        self.cube_view.on_click(self.cube_clicked)
        self.cube_highlight = self.cube_view.create_highlighter()
        
    def run(self):
        self.win.show_all()
        gtk.main()
    
    def cube_clicked(self, coords, flux):
        item_clicked = self.dendro_view.set_clicked_item_by_coords(coords)
        if (item_clicked):
            self.highlight_item_in_cube(item_clicked)
    
    def highlight_item_in_cube(self, item):
        mapdata = np.zeros(self.cube.data.shape)
        item.add_footprint(mapdata, 1)
        mapdata[mapdata>1] = 0.75 # Set the child items to be semi-transparent
        self.cube_highlight.highlight(mapdata)

if __name__ == "__main__":
    filename = False
    if len(sys.argv) == 2:
        filename = sys.argv[1]
    v = DendroViewer(filename)
    v.run()