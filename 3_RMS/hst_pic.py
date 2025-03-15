import astropy as ap
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from astropy.utils.data import get_pkg_data_filename
import matplotlib.colors as clr
import pandas as pd
import aplpy


fig = plt.figure()
figure = aplpy.FITSFigure('ngc3256_hst.tif', figure=fig)
figure.show_rgb()
figure.ticks.hide()
figure.tick_labels.hide()
figure.axis_labels.hide()
figure.show_regions('NGC3256_hst_reg')
figure.add_scalebar(length=5000/4.85/41/3600, label='5 kpc', color='white', corner='top')
x = 15*(10+27/60+46.043/3600); y = -(43+55/60+16.43/3600)
figure.show_arrows(x, y, 0, 20/3600, color='white')
figure.show_arrows(x, y, 20/3600, 0, color='white')
figure.savefig('ngc3256_hst.png', dpi=300)