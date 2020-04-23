# /usr/bin/python
# -*- coding: UTF-8 -*-
import matplotlib
import matplotlib.colors as mcolors
from matplotlib import colors
from matplotlib.colors import LogNorm,Normalize,ListedColormap
from matplotlib import colors
#####
# these colors were seleced by Hashini
# made by Menaka
# colorbars defined for Revel et al 2020.
# created by Menaka@IIS
#####
##################
def make_colormap(seq):
    """Return a LinearSegmentedColormap
    seq: a sequence of floats and RGB-tuples. The floats should be increasing
    and in the interval (0,1).
    """
    seq = [(None,) * 3, 0.0] + list(seq) + [1.0, (None,) * 3]
    cdict = {'red': [], 'green': [], 'blue': []}
    for i, item in enumerate(seq):
        if isinstance(item, float):
            r1, g1, b1 = seq[i - 1]
            r2, g2, b2 = seq[i + 1]
            cdict['red'].append([item, r1, r2])
            cdict['green'].append([item, g1, g2])
            cdict['blue'].append([item, b1, b2])
    return mcolors.LinearSegmentedColormap('CustomMap', cdict)
###################
def colormap(name):
    #cmap--
    c = matplotlib.colors.ColorConverter().to_rgb
    if name=="H01":
        color_list=[c('xkcd:blood red'),c('xkcd:darkish red'),c('xkcd:orangered'),c('xkcd:orangeish'),c('xkcd:yellowish orange'),0.33,c('xkcd:yellowish orange'),c('xkcd:sunflower'),c('xkcd:butter'),c('xkcd:light beige'),0.67,c('xkcd:light beige'),c('xkcd:seaweed'),c('xkcd:cool green'),c('xkcd:viridian')]
    if name=="H02":
        color_list=[c('xkcd:sapphire'),c('xkcd:nice blue'),c('xkcd:cool blue'),0.2,c('xkcd:cool blue'),c('xkcd:indigo blue'),c('xkcd:heather'),c('xkcd:merlot'),0.4,c('xkcd:merlot'),c('xkcd:lipstick red'),c('xkcd:pink red'),c('xkcd:watermelon'),0.6,c('xkcd:watermelon'),c('xkcd:pale salmon'),c('xkcd:pale'),c('xkcd:very pale green'),0.8,c('xkcd:very pale green'),c('xkcd:pale teal'),c('xkcd:dark seafoam'),c('xkcd:deep teal')]
    return make_colormap(color_list)
###################
