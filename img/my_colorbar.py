# /usr/bin/python
# -*- coding: UTF-8 -*-
import matplotlib
import matplotlib.colors as mcolors
from matplotlib import colors
from matplotlib.colors import LogNorm,Normalize,ListedColormap,LinearSegmentedColormap
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
    if name=="H03":
        color_list=[c('#005a74'), c('#3a8494'), c('#68b0b5'), c('#a5d5d8'), c('#ffffff'), 0.5, c('#ffffff'), c('#ffbcaf'), c('#f4777f'), c('#cf3759'), c('#93003a')]
    if name=="H04":
        color_list=[c('#fd7028'), c('#ff8d45'), c('#ffa762'), c('#ffbf80'), c('#ffd6a2'), c('#ffecc8'), c('#ffffff'), 0.5,c('#ffffff'), c('#d4e3e7'), c('#aac8cf'), c('#81acb9'), c('#5891a3'), c('#2f768f'), c('#005a7d')]
    if name=="H05":
        color_list=[c('#008300'), c('#6e9d6e'), c('#9dbd9d'), c('#cedece'), c('#ffffff'), 0.5, c('#ffffff'), c('#ffe56f'), c('#ffc15c'), c('#ff9a62'), c('#ff6b6a')]
    if name=="H06":
        color_list=[c('#005300'), c('#008300'),0.4,c('#008300'), c('white'), 0.5, c('white'),c("#ffe156"),0.6,c("#ffe156"), c("#ff6b6a")]
    if name=="H07":
        color_list=[c('#5a389e') , c('#ffffff'), 0.5, c('#ffffff'), c('#008300')]
    if name=="H08":
        color_list=[c('#fffc6f'),0.167,c('#fffc6f'), c('#36eaa5'),0.333,c('#36eaa5'), c('#53bd77'),0.5,c('#53bd77') ,c('#2b928e'),0.667,c('#2b928e'), c('#0165a7'),0.833,c('#0165a7'), c('#850031')]
    return make_colormap(color_list)
###################
def diverging_colormap(colorL='#052f61',colorR='#67001f'):
    """
    Create diverging colormap by defined colors at left and right
    """
    nodes = [0.0, 0.5, 1.0] 
    color_list = [colorL,'white',colorR] # Left color _ white _ Righ color
    cmap = LinearSegmentedColormap.from_list("", list(zip(nodes, color_list)))
    return cmap
###################