import numpy as np
import pandas as pd

import astropy.units as u
import astropy.coordinates as coord
from astropy.coordinates.sky_coordinate import SkyCoord
from astropy.units import Quantity
from pyvo.dal import TAPService

import matplotlib.pyplot as plt
from matplotlib.pyplot import cm


__tap_service_url = "http://gaia.ari.uni-heidelberg.de/tap"
__tap_service = TAPService(__tap_service_url)

#default query columns:
__column_list = ['source_id', 'ra','dec','parallax','pmra','pmdec','radial_velocity',
                    'phot_g_mean_mag','phot_bp_mean_mag', 'phot_rp_mean_mag','r_est']

class fieldstars():


    def __init__(self, name:str):
        self.name=name
        self.coords = None


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

    def __get_coords__(self, recalc=False):
        #computes and caches sky coordinates for the objects
        #set recalc=True to force recalculation
        if self.coords is None or recalc:

            self.coords = coord.SkyCoord(ra=np.array(self.objs.ra)*u.degree,
                   dec=np.array(self.objs.dec)*u.degree,
                   distance=np.array(self.objs.r_est)*u.pc,
                   pm_ra_cosdec=np.array(self.objs.pmra)*u.mas/u.yr,
                   pm_dec=np.array(self.objs.pmdec)*u.mas/u.yr,
                   radial_velocity=np.array(self.objs.radial_velocity)*u.km/u.s)

    def get_coords(self):
        #returns sky coordinates for the objects
        self.__get_coords__(self)
        return self.coords

    def maxsep(self):
        #computes maximum separation from mean of objects
        ra_mean = self.objs.ra.mean()*u.degree
        dec_mean = self.objs.dec.mean()*u.degree
        c_mean=coord.SkyCoord(ra=ra_mean, dec=dec_mean)
        seps = c_mean.separation(self.get_coords())
        return seps.max()

def from_pandas(df, colmapper, name=None):

    my_fs = fieldstars(name=name)
    src_cols = [c for c in colmapper.values()]
    dest_cols = [c for c in colmapper.keys()]

    arg_df = df.reset_index()

    # source columns in argument dataframe?
    inv_cols = set(src_cols).difference(arg_df.columns)
    if len(inv_cols) != 0:
        raise ValueError('invalid source column(s): '+str(inv_cols))

    # get the right destination columns?
    # case 1: too few dest columns given
    missing_cols = set(__column_list).difference(dest_cols)
    if len(missing_cols) != 0:
        raise ValueError('Missing column mapping for: '+str(missing_cols))
    # case 2: too many dest columns given:
    extra_cols = set(dest_cols).difference(__column_list)
    if len(extra_cols) != 0:
        raise ValueError('Invalid destination column supplied: '+str(extra_cols))

    #swap the keys and values for purposes of renaming:
    col_renamer = {s:d for s,d in zip(src_cols, dest_cols)}
    my_fs.objs = arg_df[src_cols].rename(columns = col_renamer, copy=True)

    my_fs.objs.set_index('source_id', inplace=True)
    return my_fs

if __name__ == '__main__':

    fs = fieldstars(52.074625695066345*u.degree, 48.932707471347136*u.degree, 1.0*u.degree, plx_error_thresh=5, r_est=(175,185))

    print(f'There are {len(fs.objs)} field stars')

    print(fs.tap_query_string)