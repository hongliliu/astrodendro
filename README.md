About
=====

The aim of this Python module is to make it easy to compute dendrograms
of Astronomical data. The use of dendrograms to represent Astronomical
data is described in detail in [Goodman, A. (2009,
Nature)](http://adsabs.harvard.edu/abs/2009Natur.457...63G).

**DISCLAIMER**: The code has not yet been thoroughly tested.

Screenshot
==========

In this screenshot, you can see what the generated dendrograms look like, as
returned by the plot() method described below. A branch with five leaves is
highlighted in red.

![Screenshot](http://i.imgur.com/QDePB.png)

Installing
==========

To install the ``astrodendro`` module, simply do:

    pip install git+git://github.com/bradenmacdonald/astrodendro.git#egg=astrodendro

Using
=====

Dendrograms can be computed from any 2- or 3-D Numpy array:

    >>> from astrodendro import Dendrogram
    >>> d = Dendrogram(array)

where ``array`` can be read in from a FITS file for example:

    >>> import pyfits
    >>> array = pyfits.getdata('observations.fits')

(but ``array`` can also be generated in memory, or read in from other
files, e.g. HDF5)

The tree can be explored from the trunk:

    >>> d.trunk
    [<astrodendro.components.Branch object at 0x10279b250>]

This is a list of the lowest structures. We can select one of these:

    >>> d.trunk[0]
    <astrodendro.components.Branch object at 0x10279b250>

Sub-structures of branches can then be retrieved as a list with the ``items`` attribute:

    >>> d.trunk[0].items
    [<astrodendro.components.Branch object at 0x10279b1d0>, <astrodendro.components.Leaf object at 0x10279b210>]

    >>> d.trunk[0].items[0].items
    [<astrodendro.components.Leaf object at 0x10279b150>, <astrodendro.components.Branch object at 0x10279b190>]

The pixel positions and fluxes of all the pixels in a leaf can then be retrieved with the ``coords`` and ``f`` attributes:

    >>> In [7]: d.trunk[0].items[0].coords
    [(230, 50, 75),
    (230, 50, 74),
    (229, 50, 74),
    (229, 51, 74)]
    >>> d.trunk[0].items[0].items[0].f
    [1.4287809,
    1.4096074,
    1.4536692,
    1.4319911]

It is also possible to directly retrieve all the leaves:

    >>> leaves = d.get_leaves()

Options
=======

There are several options that can be used when initializing a
``Dendrogram`` object:

* ``minimum_flux`` - the minimum flux of pixels to be used in the
  dendrogram (default is -infinity)
* ``minimum_npix`` - the minimum number of pixels necessary to create a
  new leaf (default is 0)
* ``minimum_delta`` - the minimum flux difference for two structures to
  be treated as being separate (minimum is 0)

For example:

    d = Dendrogram(array, minimum_flux=1.2, minimum_npix=10, minimum_delta=0.1)
    
Plotting
========

Plotting a Dendrogram is easy when this module is used within IPython
in pylab mode:

	d = Dendrogram(array)
	d.plot()

Viewer
======

An simple application is included which uses PyGTK, astrodendro and 
[astrocube](https://github.com/bradenmacdonald/astrocube) to allow you to
interactively generate dendrograms and see the correspondence between pixels in
the data cube and items in the dendrogram. 

![Screenshot](http://i.imgur.com/GIUwf.png)

The viewer currently only loads FITS files.

To launch the viewer, run the following command:

	astrodendro-viewer fits_file.fits

To use the viewer, specify the desired `min_flux`, `min_npix`, and `min_delta`
values, and click on the __Compute__ button. A dendrogram will be generated.

Once a dendrogram has been created, you can click on any pixel in the cube to 
highlight the dendrogram node that contains it, or conversely you can click on
any item in the dendrogram to see which pixels it includes highlighted in the
data cube.  Double-click on any item in the dendrogram to automatically jump to
the pixel with the highest flux value found among the selected item and its 
children.  

The viewer also contains an integrated IPython shell that can be used to
manipulate the dendrogram and associated data cube. Any Python commands can be
used. The following data objects are automatically added to the local scope:

* `cube`: The astrocube data cube (a `DataCube` object). You can access the cube
  data using `cube.data`

* `cube_view`: The `CubeViewWidget` used to display the data cube. You can read
  and write to `cube_view.x`, `cube_view.y`, `cube_view.z`. The `z` value gets
  or sets the z coordinate of the X-Y slice currently visible on the left side
  of the viewer, and the `x` and `y` properties get set to the coordinates of
  any data cube pixel that you click on.

* `dendro_view`: The `DendrogramViewWidget` object used to display the 
  dendrogram

* `dendrogram`: The `Dendrogram` object itself.

In addition, the following commands are defined for convenience. These commands
should be preferred over direct manipulation of the objects above.

* `create_highlighter(color)`: Creates a highlighter useful for coloring in
  both dendrogram objects and the data pixels associated with that object in
  the specified color. The returned object has the following methods:
  * `highlight_coords(coords)`: Highlight whatever dendrogram item contains
    the pixel at `coords`, which must be an (x,y,z) tuple.
  * `highlight_item(item)`: Highlight the given `Leaf` or `Branch` and the 
    associated pixels. Pass item=None to clear this highlighting.

* `set_color_map(cmap = CubeViewWidget.default_cmap)`: Use to change the color
  map that defines how the data cube is shown. For example, use
  `set_color_map("bone")` to show the data cube in greyscale so that colored
  highlights are more visible. Accepts any matplotlib
  [color map](http://www.scipy.org/Cookbook/Matplotlib/Show_colormaps).

* `set_data(data)`: Replace the data in the data cube with the given numpy
  ndarray, and clear the dendrogram. Useful for adding noise to the data, etc.

* `make_dendrogram(min_flux, min_npix, min_delta)`: Compute a new dendrogram.

* `export_pdf(filename, title="default", subtitle="default",
  include_cube_slice=True)`: Export a two-page PDF, with the currently visible
  slice of the data cube on the first page, and the current dendrogram on the
  second page. A title and subtitle including the dendrogram parameters will
  automatically added, but you can change them using the `title` and `subtitle`
  arguments. Set `include_cube_slice=False` if you only want to export the
  dendrogram.

* `export_png(filename_cube, filename_dendro)`: Essentially exports a
  screenshot of the data cube and/or the dendrogram. Size of the exported
  PNG images depends on the current window size of the viewer program.
  Either .png filename may be `None` if you don't want to export one or the
  other. In general, PDF export is recommended over .PNG export.

Import/Export
=============

Dendrograms can be written out and read in from a portable HDF5-based format:

    # Create dendrogram and write it out
    d = Dendrogram(array)
    d.to_hdf5('observations_dendrogram.hdf5')

    # Read dendrogram into new instance
    d2 = Dendrogram()
    d2.from_hdf5('observations_dendrogram.hdf5')

Unit Tests and Benchmarks
=========================

Several unit tests are included, and are also installed with the package.
To run the unit tests, simply run the command

	python -m astrodendro.test

A benchmark is also included that uses realistic data to determine how fast the
dendrograms are being generated. To run the benchmark, use the command

	python -m astrodendro.test.benchmark
