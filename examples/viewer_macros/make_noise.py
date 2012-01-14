# Example macro for astrodendro-viewer by Braden MacDonald
#
# This will add random noise to the data cube some desired number of times,
# exporting a PDF figure each time.
# 
# Simply paste this whole file into the IPython shell built into 
# astrodendro-viewer and press enter twice.

import numpy as np

##############################################
# Parameters (min flux, min npix, min delta):
min_flux = 1.4
min_npix = 10
min_delta = 0.3
noise_dev = 0.1
n = 10 # n+1 PDF files will be output in the current directory
cube_view.z = 212 # Set this to an interesting Z coordinate
# Set up highlighting of four specific pixels in the data:
highlighters = [create_highlighter(c) for c in ["purple", "green", "orange", "brown"]]
coords = [(51,69,213), (69,7,214), (3,14,212), (96,32,212)] # First coord is purple, etc.
##############################################

set_color_map("bone")

fig_num = 1

def generate(min_flux,min_npix,min_delta, desc=""):
    " Make a new dendrogram, and export it as a PDF "
    global fig_num
    make_dendrogram(min_flux=min_flux, min_npix=min_npix, min_delta=min_delta)
    for i in range(0,len(highlighters)):
        highlighters[i].highlight_coords(coords[i])
    title="Fig. {num}".format(num=fig_num)
    if desc:
        desc += " "
    desc=desc+"flux {f} npix {n} delta {d}".format(f=min_flux, n=min_npix, d=min_delta)
    filename="fig{fig:02}_{desc}.pdf".format(fig=fig_num, desc=desc.replace(" ", "_"))
    export_pdf(filename, title=title, subtitle=desc)
    fig_num += 1

# Here is where we actually create the dendrograms:
orig_data = cube.data

# First create a dendrogram with no added noise:
generate(min_flux, min_npix, min_delta, desc="No extra noise");

# Now, add noise, create the dendrogram, export a PDF, and repeat n times.
for i in range(1,n):
    set_data(orig_data + np.random.normal(scale=noise_dev, size=cube.data.shape))
    generate(min_flux, min_npix, min_delta, desc="Noise {0}".format(noise_dev));