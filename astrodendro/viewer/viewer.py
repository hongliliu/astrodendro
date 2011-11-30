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




if __name__ == "__main__":
    
    filename = False
    if len(sys.argv) == 2:
        filename = sys.argv[1]
    
    cube = DataCube(filename, calc_noise_dev=False)
    
    win = gtk.Window()
    win.connect("destroy", lambda x: gtk.main_quit())
    win.set_default_size(800,600)
    win.set_title("{ln} spectral line map of {o}".format(o=cube.object_name,ln=cube.line_name))
    
    # Create our three main widgets:
    
    # A CubeViewWidget and a connected DendrogramViewWidget in the top row:
    cube_view = CubeViewWidget(cube, win)
    dendro_view = DendrogramViewWidget(cube_view, win)
    first_row = gtk.HPaned()
    first_row.pack1(cube_view, resize=True, shrink=True)
    first_row.pack2(dendro_view, resize=True, shrink=True)
    
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
    
    win.show_all()
    gtk.main()