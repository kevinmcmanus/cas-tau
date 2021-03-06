import pickle
from os import path
import re
import numpy as np

import astropy.units as u
import astropy.coordinates as coord
from astropy.coordinates.sky_coordinate import SkyCoord
from astropy.units import Quantity

import matplotlib.pyplot as plt
from matplotlib.pyplot import cm

class cluster():

    def __init__(self, clustname:str, rootdir:str, coldefsname=None):
        
        pathname = path.join(rootdir, clustname+'.pkl')
        with open(pathname,'rb') as f:
            self.objs = pickle.load(f)

        defname = path.join(rootdir,
            coldefsname+'_definition.pkl' if coldefsname is not None else clustname+ '_definition.pkl')
        with open(defname,'rb') as f:
            self.coldefs = pickle.load(f)

        self.cluster_name = clustname

    def findcolumns(self, rex):
        target = '.*' + rex + '.*'
        cols =   [(c, d['desc'])for c,d in self.coldefs.items() if re.match(target,d['desc'],flags=re.IGNORECASE)]
        return cols

    def plot_radec(self, ax=None):
        if ax is not None:
            yax = ax
        else:
            yax = plt.subplot(111)
        
        yax.scatter(self.objs.RAdeg, self.objs.DEdeg, s=1, label=self.cluster_name, color='blue')
        if ax is None:
            yax.legend()

    def plot_hrdiagram(self, **kwargs):
        ax = kwargs.get('ax')
        title = kwargs.get('title')
        color = kwargs.get('color')

        if color is None:
            color='blue'

        if ax is not None:
            yax = ax
        else:
            yax = plt.subplot(111)

        #distmod = 5*np.log10(self.objs.rest)-5
        distance = coord.Distance(parallax=u.Quantity(np.array(self.objs.Plx)*u.mas),allow_negative=True)

        abs_mag = self.objs.Gmag - distance.distmod.value
        BP_RP = self.objs.BPmag - self.objs.RPmag

        yax.scatter(BP_RP,abs_mag, s=1,  label=self.cluster_name, color=color)
        if not yax.yaxis_inverted():
            yax.invert_yaxis()
        yax.set_xlim(-1,5)
        #yax.set_ylim(20, -1)

        yax.set_title(title)
        if ax is None:
            yax.legend()



         


if __name__ == '__main__':

    alpha_per = cluster('alpha_per','./data/J_A+A_628_A66')

    print(f'Alpha-Persei has {len(alpha_per.objs)} members')

    cols = alpha_per.findcolumns('Gaia')
    for c in cols:
        print(f'Column: {c[0]}, Description: {c[1]}')    
