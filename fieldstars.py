import numpy as np
import pandas as pd

import astropy.units as u
import astropy.coordinates as coord
from astropy.coordinates.sky_coordinate import SkyCoord
from astropy.units import Quantity
from pyvo.dal import TAPService

import matplotlib.pyplot as plt
from matplotlib.pyplot import cm

class fieldstars():

    __tap_service_url = "http://gaia.ari.uni-heidelberg.de/tap"
    __tap_service = TAPService(__tap_service_url)

    #default query columns:
    __column_list = ['source_id', 'ra','dec','parallax','pmra','pmdec','radial_velocity',
                     'phot_g_mean_mag','phot_bp_mean_mag', 'phot_rp_mean_mag','r_est']

    def __init__(self, name:str):
        self.name=name


    @u.quantity_input(ra='angle', dec='angle', rad='angle')
    def conesearch(self, ra, dec, radius, **kwargs):

        # input parameters in degrees:
        ra_ = ra.to(u.degree); dec_= dec.to(u.degree); rad_=radius.to(u.degree)

        r_est = kwargs.get('r_est')
        plx_error_thresh = kwargs.get('plx_error_thresh')

        columnlist = ' '+ '\n\t\t,'.join(self.__column_list)
        dbsource = '\n\t'.join(['\nFROM gaiadr2_complements.geometric_distance gd',
                    'INNER JOIN gaiadr2.gaia_source gs using (source_id) '])

        constraints =  '\n\t'.join(['\nWHERE ', 
                'CONTAINS(POINT(\'\', gs.ra, gs.dec), ',
                '\tCIRCLE(\'\', {ra}, {dec}, {rad})) = 1 '.format(ra=ra_.value, dec=dec_.value, rad=rad_.value)])

        if r_est is not None:
            constraints = constraints + '\n\tAND gd.r_est BETWEEN {r_lo} AND {r_hi}'.format(r_lo=r_est[0],r_hi=r_est[1])

        if plx_error_thresh is not None:
            constraints = constraints + '\n\tAND parallax_over_error >= {thresh}'.format(thresh=plx_error_thresh)

        self.tap_query_string = 'SELECT \n\t\t'+ columnlist + dbsource + constraints

        tap_service = TAPService(self.__tap_service_url)

        tap_results = tap_service.search(self.tap_query_string)

        self.objs = tap_results.to_table().to_pandas()
        self.objs.set_index('source_id', inplace=True)

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

        distmod = 5*np.log10(self.objs.r_est)-5
        #distance = coord.Distance(parallax=u.Quantity(np.array(self.objs.Plx)*u.mas),allow_negative=True)

        abs_mag = self.objs.phot_g_mean_mag - distmod
        BP_RP = self.objs.phot_bp_mean_mag - self.objs.phot_rp_mean_mag

        yax.scatter(BP_RP,abs_mag, s=1,  label=self.name, color=color)
        if not yax.yaxis_inverted():
            yax.invert_yaxis()
        yax.set_xlim(-1,5)
        #yax.set_ylim(20, -1)

        yax.set_title(title)
        if ax is None:
            yax.legend()

if __name__ == '__main__':

    fs = fieldstars(52.074625695066345*u.degree, 48.932707471347136*u.degree, 1.0*u.degree, plx_error_thresh=5, r_est=(175,185))

    print(f'There are {len(fs.objs)} field stars')

    print(fs.tap_query_string)