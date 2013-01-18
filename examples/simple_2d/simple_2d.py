import pyfits
from astrodendro import Dendrogram

# Read in the data
image = pyfits.getdata('data.fits.gz')

# Compute the dendrogram
d = Dendrogram(image)

# Output to HDF5 file
d.to_hdf5('simple_2d_dendrogram.hdf5')

# plot in some different ways
import matplotlib as mpl
from matplotlib import pyplot as pl
max_npix = max(i.npix for i in d.all_items)
min_npix = min(i.npix for i in d.all_items)
color_lambda = lambda value: mpl.cm.jet(((value-min_npix)/float(max_npix-min_npix)))

pl.figure(1)
pl.clf()
p1=d.plot(style='fluxpix',brightest=True,color=color_lambda)
p1.axes.set_xscale('log')
p1.axes.set_yscale('log')

pl.figure(2)
pl.clf()
p2=d.plot(style='fluxpix',brightest=False,color=color_lambda)
p2.axes.set_xscale('log')
p2.axes.set_yscale('log')

pl.show()
