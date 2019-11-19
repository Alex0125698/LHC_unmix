from functools import partial
import logging
from collections import OrderedDict

import numpy as np
import scipy.sparse
# from osgeo import gdal
# matplotlib inline
import matplotlib
import matplotlib.pyplot as plt
matplotlib.rc('font',size=15)

import proxmin
from proxmin import nmf
from proxmin.operators import prox_soft, get_gradient_x, get_gradient_y
from proxmin.utils import l2
# import dc

plogger = logging.getLogger("proxmin")
plogger.setLevel(logging.DEBUG)

points = OrderedDict([
    ("concrete", (162, 507)), # statue
    #("road", (175, 546)), # road, degenerate with concrete
    #("roof2", (114, 575)), # roof, degenerate with dirt
    ("dirt", (172, 597)), # dirt
    ("grass", (125, 575)), #grass
    #("trees", (183, 561)), #tree
    ("bkg", (None,None))
])

#color_cycle = [dc.ref_colors[obj] for obj in points.keys()]
#matplotlib.rcParams['axes.color_cycle'] = color_cycle

ds = gdal.Open('/Users/fred/Downloads/Hyperspectral_Project/dc.tif')
data_shape = ds.GetRasterBand(1).ReadAsArray().shape

# only use a subset of the image
shape = (ds.RasterCount, data_shape[0]*data_shape[1])

# Get hyperspectral data
data = np.zeros(shape)
for bidx in range(shape[0]):
    band = ds.GetRasterBand(bidx + 1).ReadAsArray()
    data[bidx] = band.flatten()
dc.plot_color_img(data, data_shape, figsize=(8,8), show=False);

# only use a subset of the image
xmin = 50
xmax = 250
ymin = 400
ymax = 600
img_shape = (ymax-ymin, xmax-xmin)
img = data.reshape(data.shape[0], data_shape[0], data_shape[1])[:,ymin:ymax, xmin:xmax]
img = img.reshape(data.shape[0], img_shape[0]*img_shape[1])

spectra = OrderedDict()
for obj, pt in points.items():
    if obj!="bkg":
        plt.plot(pt[0], pt[1],'rx', label=obj)
        spectra[obj] = dc.get_point_spec(pt[0], pt[1], data, data_shape)
    else:
        spectra[obj] = np.min(img, axis=1)

for obj, spec in spectra.items():
    if obj!="bkg":
        spectra[obj] = spec-spectra["bkg"]
plt.legend()
plt.xlim([50,250])
plt.ylim([600,400])
plt.show()

for obj, pt in points.items():
    plt.plot(spectra[obj], label=obj)
plt.legend()
plt.show()

# Get wavelengths used in hyperspectral data
wavelength_data = np.recfromcsv('/Users/fred/Downloads/Hyperspectral_Project/wavelengths.txt', delimiter=" ")
wavelength = wavelength_data["wavelength"]
idx = wavelength_data["idx"]