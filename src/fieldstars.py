import numpy as np
import pandas as pd
import pickle

import astropy.units as u
import astropy.coordinates as coord
from astropy.coordinates.sky_coordinate import SkyCoord
from astropy.table import Table
from astropy.units import Quantity
from astroquery.gaia import Gaia
#from pyvo.dal import TAPService

import matplotlib.pyplot as plt
from matplotlib.pyplot import cm





#default query columns:


class fieldstars():

    column_list = ['source_id', 'ra','dec','parallax','pmra','pmdec','radial_velocity',
                    'phot_g_mean_mag','phot_bp_mean_mag', 'phot_rp_mean_mag', 'e_bp_min_rp_val', 'a_g_val'] #,'r_est'] # r_est not in gaia archive
    tap_service_url = "http://gaia.ari.uni-heidelberg.de/tap" #need to change to use just gaia archive
    
    source_constraints = ' AND '.join([
        'parallax_over_error > 10',
        'phot_g_mean_flux_over_error>50',
        'phot_rp_mean_flux_over_error>20',
        'phot_bp_mean_flux_over_error>20',
        'phot_bp_rp_excess_factor < 1.3+0.06*power(phot_bp_mean_mag-phot_rp_mean_mag,2)',
        'phot_bp_rp_excess_factor > 1.0+0.015*power(phot_bp_mean_mag-phot_rp_mean_mag,2)',
        'visibility_periods_used>8',
        'astrometric_chi2_al/(astrometric_n_good_obs_al-5)<1.44*greatest(1,exp(-0.4*(phot_g_mean_mag-19.5)))'])

    def __init__(self, name:str):
        self.name=name
        self.coords = None
        self.tap_query_string = None

    def _get_col_list(self, prefix='gs'):
        collist = f'{prefix}.' +  f'\n\t\t,{prefix}.'.join(self.column_list)
        return collist
    
    @u.quantity_input(ra='angle', dec='angle', rad='angle')
    def conesearch(self, ra, dec, radius, **kwargs):
        ## ToDo: rewrite to use Gaia Archive, not tap service.  Need to ditch distance param and just use parallax
        # input parameters in degrees:
        ra_ = ra.to(u.degree); dec_= dec.to(u.degree); rad_=radius.to(u.degree)

        maxrec = kwargs.get('maxrec', 20000)

        columnlist = self._get_col_list() + ', gd.r_est'
        dbsource = '\n\t'.join(['\nFROM external.gaiadr2_geometric_distance gd',
                    'INNER JOIN gaiadr2.gaia_source gs using (source_id) '])

        constraints =  '\n\t'.join(['\nWHERE ', 
                'CONTAINS(POINT(\'\', gs.ra, gs.dec), ',
                '\tCIRCLE(\'\', {ra}, {dec}, {rad})) = 1 '.format(ra=ra_.value, dec=dec_.value, rad=rad_.value)])

        if self.source_constraints is not None:
            constraints = constraints + ' AND '+ self.source_constraints

        self.tap_query_string = 'SELECT \n\t\t'+ columnlist + dbsource + constraints
        
        #tap_service = TAPService(self.tap_service_url)
        #tap_results = tap_service.search(self.tap_query_string, maxrec=maxrec)
        #self.objs = tap_results.to_table().to_pandas()
        #fetch the data
        
        job = Gaia.launch_job_async(query=self.tap_query_string)
        self.objs = job.get_results().to_pandas()

        self.objs.set_index('source_id', inplace=True)
        
    def from_source_idlist(self, source_idlist, source_idcol=None, filters=False):
        #xml-ify the source_idlist to a file
        if isinstance(source_idlist, Table):
            #guess which column contains the source ids
            if source_idcol is None:
                if 'source_id' in source_idlist.colnames:
                    sidcol = 'source_id'
                elif 'source' in source_idlist.colnames:
                    sidcol = 'source'
                else:
                    raise ValueError('no column to use as source_id')
            else:
                if source_idcol not in source_idlist.colnames:
                    raise ValueError(f'invalid column specified as source id column: {source_idcol}')
                sidcol = source_idcol
                
            tbl = source_idlist
        elif isinstance(source_idlist, np.ndarray) or isinstance(source_idlist, list):
            sidcol = 'source_id'
            tbl = Table({sidcol:source_idlist})
        else:
            raise ValueError(f'invalid source_idlist type: {type(source_idlist)}')
        xml_path = 'source_idlist.xml'
        tbl.write(xml_path, table_id='source_idlist', format='votable', overwrite=True)
        
        #build the query:
        col_list = self._get_col_list() + ', gd.r_est'
        
        dbsource =  ''.join([' FROM tap_upload.source_idlist sidl',
                            f' LEFT JOIN gaiadr2.gaia_source gs ON gs.source_id = sidl.{sidcol}',
                            f' LEFT JOIN external.gaiadr2_geometric_distance gd on gs.source_id = gd.source_id' ])
        
        query_str = f'SELECT sidl.{sidcol} as "source", '+col_list+dbsource
        if filters:
            query_str = query_str + ' WHERE '+ self.source_constraints
        self.tap_query_string = query_str
        
        #fetch the data
        job = Gaia.launch_job_async(query=query_str, upload_resource=xml_path,upload_table_name='source_idlist')
        self.objs = job.get_results().to_pandas()
        
        #fetch data via tap query
        #tap_service = TAPService(self.tap_service_url)
        #tap_results = tap_service.search(self.tap_query_string, maxrec=len(tbl),uploads={'source_idlist':tbl})
        #self.objs = tap_results.to_table().to_pandas()
        
        
        self.objs.set_index('source', inplace=True)
                
    
    def merge(self, right):
        # joins right fieldstars to self; returns result
        consol_df = self.objs.merge(right.objs,left_index=True, right_index=True, how='outer', indicator=True)

        consol_df['which'] = consol_df._merge.apply(lambda s: right.name if s == 'right_only' else self.name if s == 'left_only' else 'both')
        
        #fix up columns; preferential treatment for self's columns
        mycols = self.objs.columns
        for c in mycols:
            consol_df[c]=np.where(np.isnan(consol_df[c+'_x']), consol_df[c+'_y'], consol_df[c+'_x'])

        consol_df.drop(columns = [s+'_x' for s in mycols]+[s+'_y' for s in mycols], inplace=True)
        
        my_fs = fieldstars(name = self.name + ' merged with ' + right.name)
        my_fs.objs = consol_df
        my_fs.tap_query_string = [self.tap_query_string, right.tap_query_string]
        
        return my_fs
    
    def to_pickle(self, picklepath):
        """
        pickles and dumps the object to file designated by picklepath
        """
        with open(picklepath,'wb') as pkl:
            pickle.dump(self, pkl)
        

    def plot_hrdiagram(self, **kwargs):
        ax = kwargs.get('ax')
        title = kwargs.get('title', 'HR Diagram')
        color = kwargs.get('color', 'blue')
        alpha = kwargs.get('alpha', 1.0)      
        absmag = kwargs.get('absmag', True)
        
        if ax is None:
            yax = plt.subplot(111)
        else:
            yax = ax
            
        distmod = 5*np.log10(self.objs.r_est)-5
        #distmod = (10.0 - 5.0*np.log10(self.objs.parallax)) if absmag else 0.0
        #distance = coord.Distance(parallax=u.Quantity(np.array(self.objs.Plx)*u.mas),allow_negative=True)

        abs_mag = self.objs.phot_g_mean_mag - distmod
        BP_RP = self.objs.phot_bp_mean_mag - self.objs.phot_rp_mean_mag

        yax.scatter(BP_RP,abs_mag, s=1,  label=self.name, color=color)
        if not yax.yaxis_inverted():
            yax.invert_yaxis()
        yax.set_xlim(-1,5)
        #yax.set_ylim(20, -1)

        yax.set_title(title)
        yax.set_ylabel('$M_G$',fontsize=14)
        yax.set_xlabel('$G_{BP}\ -\ G_{RP}$', fontsize=14)
        if ax is None:
            yax.legend()
            
    def plot_motions(self, **kwargs):
        from scipy.stats import kde
        ax = kwargs.get('ax')
        title = kwargs.get('title', 'Proper Motions')
        cmap = kwargs.get('cmap', 'viridis')
        nbins = kwargs.get('nbins', 300)
   
        if ax is None:
            fig, yax = plt.subplots()
        else:
            yax = ax
            
        x = np.array(self.objs.pmra)
        y = np.array(self.objs.pmdec)
        
        # Evaluate a gaussian kde on a regular grid of nbins x nbins over data extents
        k = kde.gaussian_kde([x,y])
        xi, yi = np.mgrid[x.min():x.max():nbins*1j, y.min():y.max():nbins*1j]
        zi = k(np.vstack([xi.flatten(), yi.flatten()]))
        
        pcm = yax.pcolormesh(xi, yi, zi.reshape(xi.shape),cmap=cmap)
        
        yax.set_title(title)
        yax.set_ylabel('PM Dec (mas/yr)')
        yax.set_xlabel('PM RA (mas/yr)')
        
        if ax is None:
            fig.colorbar(pcm)
        
        return pcm
    
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

    def from_pandas(self, df, colmapper):

        src_cols = [c for c in colmapper.values()]
        dest_cols = [c for c in colmapper.keys()]

        arg_df = df.reset_index()

        # source columns in argument dataframe?
        inv_cols = set(src_cols).difference(arg_df.columns)
        if len(inv_cols) != 0:
            raise ValueError('invalid source column(s): '+str(inv_cols))

        # get the right destination columns?
        # case 1: too few dest columns given
        missing_cols = set(self.column_list).difference(dest_cols)
        if len(missing_cols) != 0:
            raise ValueError('Missing column mapping for: '+str(missing_cols))
        # case 2: too many dest columns given:
        extra_cols = set(dest_cols).difference(self.column_list)
        if len(extra_cols) != 0:
            raise ValueError('Invalid destination column supplied: '+str(extra_cols))

        #swap the keys and values for purposes of renaming:
        col_renamer = {s:d for s,d in zip(src_cols, dest_cols)}
        self.objs = arg_df[src_cols].rename(columns = col_renamer, copy=True)

        self.objs.set_index('source_id', inplace=True)
        self.coords = None
        
        return


def from_pickle(picklefile):
    """
    reads up a pickle file and hopefully it's a pickled fieldstars object
    """
    with open(picklefile,'rb') as pkl:
        my_fs = pickle.load(pkl)
        
    if not isinstance(my_fs, fieldstars):
        raise ValueError(f'pickle file contained a {type(my_fs)}')
        
    return my_fs

if __name__ == '__main__':

    fs = fieldstars(52.074625695066345*u.degree, 48.932707471347136*u.degree, 1.0*u.degree, plx_error_thresh=5, r_est=(175,185))

    print(f'There are {len(fs.objs)} field stars')

    print(fs.tap_query_string)