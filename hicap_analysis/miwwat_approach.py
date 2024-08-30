from .geoprocessing import Geoprocess, WellGeometry
import geopandas as gpd
import numpy as np
import os
import pandas as pd
import shapely
from pathlib import Path
import sys
sys.path.append('../')
from .wells import Well

# these were needed to pass the testing, reprojection was giving (Inf, Inf)
import pyproj
pyproj.network.set_network_enabled(False)

class MiWWAT(Geoprocess):
    ''' Subclass to mimic the steps of the Michigan 
        Water Withdraw Assessment Tool (WWAT) online
        screening tool.  Inherits from Geoprocess class
        but recognizes that shapefiles supporting the
        MI WWAT have initial estimates of transmissivity
        and streambed conductance.  The subclass
        simplifies making a stand-alone script that
        mimics the on-line tool.
    '''
    def __init__(self, 
                 catchment_shp=None, 
                 streamline_shp=None, 
                 catch_idx=None,
                 stream_idx=None):
        ''' Constructor for WWAT approach removes
            unused columns from the WWAT shapefiles
            to clean up later processing.
        '''
        super().__init__(catchment_shp, 
                        streamline_shp, 
                        catch_idx,
                        stream_idx)
        
        # limit returned dataframe to desired columns
        self.catchment_df = self.catchment_df[['BASIN', 'SUBBASIN', 'EST_DPH_BR',
       'BDRK_FLAG', 'BDRK_T', 'MEDIAN_T', 'EST_Kv_W', 'TYPE',  'geometry']]
        
        self.streamline_df = self.streamline_df[['BASIN', 'SUBBASIN', 'EST_DPH_BR',
       'BDRK_FLAG', 'BDRK_T', 'MEDIAN_T', 'EST_Kv_W', 'TYPE',  'geometry']]
        
        return
        

    def get_geometries(self, well, wellnamefield='name',
                       latitudefield='lat', longitudefield='long'):
        ''' given a well location or list of well locations
            DataFrame, or GeoDataFrame; 
            get the home and neighboring catchments.
            Also get the distances to the closest stream
            in each catchment and the inverse distance
            and inverse distance squared approtionments.

            Parameters
            ----------
            well: lat/long location as list, list of lat/long lists, 
                  dataframe with well names and locations,
                  or geodataframe with that information.
            wellnamefield: string
                if a DataFrame or GeoDataFrame is passed, the name
                of the column with well name or ID to be put in
                WellGeometry object.  Defaults to 'name'.
            latitudefield: string
                if a DataFrame or GeoDataFrame is passed, the name
                of the column with latitude of the well.
                Defaults to 'lat'.
            longitudefield: string
                if a DataFrame or GeoDataFrame is passed, the name
                of the column with latitude of the well.
                Defaults to 'long'.
            

            Returns
            -------
            wellgeom: list of WellGeometry objects

        
        '''
        if isinstance(well, list):
            well = np.array(well)

        if isinstance(well, np.ndarray):
            if not isinstance(well[0], np.ndarray):   # list of lists
                well = np.array([well])
            wells = []
            count=0
            for w in well:
                welldict = {'name': f'testwell_{count}',
                            'lat': w[0],
                            'long': w[1]}
                count = count + 1
                wells.append(welldict)

            welltemp = pd.DataFrame(wells)
            well_df = gpd.GeoDataFrame(welltemp,
                                       geometry=gpd.points_from_xy(welltemp['long'],
                                                                   welltemp['lat']),
                                        crs='EPSG:4326')
            
            well_df = well_df.to_crs(self.streamline_df.crs)

        elif isinstance(well, gpd.GeoDataFrame):
            well_df = well.to_crs(self.streamline_df.crs)


        elif isinstance(well, pd.DataFrame):
            well_df = gpd.GeoDataFrame(well,
                                       geometry=gpd.points_from_xy(well[longitudefield],
                                                                   well[latitudefield]),
                                        crs='EPSG:4326')
            well_df = well_df.to_crs(self.streamline_df.crs)

        else:
            sys.exit('unknown parameter passed to get_geometries')

        # given the geodataframe with well locations, get the home catchment
        # and neighbors for each, get the stream reaches and closest points,
        # compute the apportionments.  Put the information into a list
        # of WellGeometry objects and return the list.
            
        well_df.set_index(wellnamefield, inplace=True)

        # check that well_df has depth information, needs to have it for WWAT
        # implementation

        if not 'depth' in well_df.columns:
            sys.exit('WWAT class needs to have depth defined for wells')

        well_geometries = []

        for n, w in well_df.iterrows():
            w_df = gpd.GeoDataFrame(w.to_frame().T, crs=well_df.crs)
            home_df = self.get_home(w_df)
            neighbors_df = self.get_neighbors(home_df)

            target_streams = neighbors_df.join(self.streamline_df, lsuffix='l_')
            target_streams.drop(columns='geometryl_', inplace=True) # drop catchment geometry
            
            nearest = []
            for k, v in target_streams.iterrows():
                d = dict()
                d['ADJ_SEGMNT'] = k
                d['distance'] = w_df.iloc[0].geometry.distance(v.geometry) * 3.28 # convert to ft
                p1, p2, = shapely.ops.nearest_points(v.geometry, 
                                                     w_df.iloc[0].geometry)
                d['geometry'] = p1
                nearest.append(d)
            nearest_df = gpd.GeoDataFrame(nearest)
            nearest_df.set_index('ADJ_SEGMNT', inplace=True)
            nearest_df['inv_distance'] = 1./nearest_df['distance']
            nearest_df['inv_dist2'] = nearest_df['inv_distance'].apply(lambda x: x**2)
            nearest_df['apportionment'] = nearest_df['inv_distance']/nearest_df['inv_distance'].sum()
            nearest_df['apport2'] = nearest_df['inv_dist2']/nearest_df['inv_dist2'].sum()

            # for WWAT include transmissivity, storativity and streambed conductance
            # check bedrock flag and depth when deciding bedrock or glacial
            # well, include aquifer type in WellGeometry

            # see if depth puts well into bedrock

            if (home_df['BDRK_FLAG'].values[0] == 2) and (well_df['depth'] > home_df['EST_DPH_BR'].values[0]):
                # bedrock well
                aq_type = 'bedrock'
                stor = 0.0004
                trans = home_df['BDRK_T'].values[0]
                # take conductance as Kv*W for drift and stream divided by the depth to bedrock
                vert_conductance = home_df['EST_Kv_W'].values[0]/home_df['EST_DPH_BR'].values[0]
            else:
                # drift well or using drift parameters
                aq_type = 'drift'
                stor = 0.01
                trans = home_df['MEDIAN_T'].values[0]
                vert_conductance = home_df['EST_Kv_W'].values[0]/w_df['depth'].values[0]

            well_geometries.append(WellGeometry(name= n, 
                                                well_df = w_df,
                                                home_df =  home_df,
                                                neighbors_df =  neighbors_df,
                                                streams_df =  target_streams,
                                                close_points_df =  nearest_df,
                                                transmissivity= trans,
                                                storativity= stor,
                                                streambed_cond= vert_conductance,
                                                aq_type= aq_type))
        return well_geometries

        






