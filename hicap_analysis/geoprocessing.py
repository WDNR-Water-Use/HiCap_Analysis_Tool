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

class Geoprocess():
    ''' Object to find catchment containing a point, neighboring
        catchments, and closest streams.  Start with a generic 
        class that can then be made specific for different 
        use cases.
    '''

    def __init__(self, 
                 catchment_shp=None, 
                 streamline_shp=None, 
                 catch_idx=None,
                 stream_idx=None):
        ''' 
            Parameters
            ----------
            catchment_shp: string or Path
                path and filename for shapefile with catchments 
                used in the analysis.  String should be passed
                with right format to be cast as a pathlib.Path()

            streamline_shp: string or Path
                path and filename for shapefile with streams
                used in the analysis. String should be passed
                with right format to be cast as a pathlib.Path()

            catch_idx: string
                name of column used to set_index for 
                geodataframe of catchments. Defaults to None which
                gives sequential indexing from shapefile.

            stream_idx: string
                name of column used to set_index for
                geodataframe of streamlines. Defaults to None which
                gives sequential indexing from shapefile.
            
            Attributes
            ----------
            catchment_df: geopandas DataFrame
                dataframe of catchment polygons
            
            streamline_df: geopandas DataFrame
                dataframe of stream lines

            Methods
            -------
            get_home(self, Well or lat/long): returns geopandas dataframe
                 polygon for catchment containing a passed Well object, geopandas
                 dataframe, or (lat/long) location - usually called by get_geometries

            get_neighbors(self, Well or lat/long or gpd): returns geopandas dataframe
                 of polygons for catchments neighboring a home catchment given
                 either a geopandas dataframe with the home catchment, a Well
                 object, an (lat/long) location - usually called by get_geometries

            get_geometries(self, gpd): returns two geopandas dataframes 
            
                The first has streamlines
                contained in the home+neighbors. 
                
                The second has the location of the closest point
                for the stream in each catchment and the well, finally the 
                inverse-distance and inverse-distance-squared apportionment
                for the neighboring closest points.

            Returns
            -------
            None
        '''
        if catchment_shp is None:
            sys.exit('need to provide catchment shapefile')
        if streamline_shp is None:
            sys.exit('need to provide streamline shapefile')

        if not isinstance(catchment_shp, Path):
            catchment_shp = Path(catchment_shp)

        if not isinstance(streamline_shp, Path):
            streamline_shp = Path(streamline_shp)

        self.catchment_df = gpd.read_file(catchment_shp)
        self.streamline_df = gpd.read_file(streamline_shp)

        # make sure they are in the same projection and if not
        # reproject streamlines into catchment projection

        if not self.catchment_df.crs == self.streamline_df.crs:
            self.streamline_df = self.streamline_df.to_crs(self.catchment_df.crs)

        if catch_idx is not None:
            self.catchment_df[catch_idx] = self.catchment_df[catch_idx].astype(int)
            self.catchment_df.set_index(catch_idx, inplace=True)
        if stream_idx is not None:
            self.streamline_df[stream_idx] = self.streamline_df[stream_idx].astype(int)
            self.streamline_df.set_index(stream_idx, inplace=True)

        return
    
    def get_home(self, well):
        ''' Return the catchment that contains a well.

            Parameters
            ----------
            well: geopandas dataframe, location in lat/long as an array or list,
                  or a hicap_analysis Well() object 


            Returns
            -------
            geopandas dataframe with one entry for the home catchment.
        
        '''
        # if its a list, first cast to np.array, then 
        # check if it is a np.array, geodataframe, or Well() object
        # return an error if neither of those
        if isinstance(well, list):
            well = np.array(well)

        if isinstance(well, np.ndarray):
            welldict = [{'name':'testwell',
                         'lat': well[0],
                         'long': well[1]}]
            welltemp = pd.DataFrame(welldict)
            well_df = gpd.GeoDataFrame(welltemp,
                                       geometry=gpd.points_from_xy(welltemp.long,
                                                                   welltemp.lat),
                                        crs='EPSG:4326')
            
            well_df = well_df.to_crs(self.catchment_df.crs)

        elif isinstance(well, gpd.GeoDataFrame):
            well_df = well.to_crs(self.catchment_df.crs)

        elif isinstance(well, Well):
            if Well.x is not None:
                welldict = [{'name':'testwell',
                         'x': well.x,
                         'y': well.y}]
                welltemp = pd.DataFrame(welldict)

                #assume if x,y is given the project is the same as the catchment
                well_df = gpd.GeoDataFrame(welltemp,
                                            geometry=shapely.Point(well.x, well.y),
                                            crs=self.catchment_df.crs)
                
            elif Well.long is not None:
                welldict = [{'name':'testwell',
                         'lat': well.lat,
                         'long': well.long}]
                welltemp = pd.DataFrame(welldict)
                well_df = gpd.GeoDataFrame(welltemp,
                                       geometry=gpd.points_from_xy(welltemp.long,
                                                                   welltemp.lat),
                                        crs='EPSG:4326')
            
                well_df = well_df.to_crs(self.catchment_df.crs)

            else:
                sys.exit('Well object needs to have x,y or lat,long in get_home()')
        else:
            sys.exit('Need to pass either (lat,long), geodataframe, or Well object')

        home_catch = gpd.sjoin(self.catchment_df, well_df, predicate='contains')
        home_catch.drop(columns='index_right', inplace=True)

        return home_catch

    def get_neighbors(self, home_info):
        ''' Return the neighboring catchments to the home catchment
            as a geodataframe.

            Parameters
            ----------
            home_info: either a geopandas dataframe with the home catchment, a Well
                        object, or a lat/long location. If it is a Well object
                        or lat/long location, then get_home() is called first. 


            Returns
            -------
            geopandas dataframe with neighbors to the home catchment.
        
        '''
        # if its a list, first cast to np.array, then 
        # check if it is a np.array, geodataframe, or Well() object
        # return an error if neither of those
        if isinstance(home_info, list):
            home_catch = self.get_home(home_info)

        elif isinstance(home_info, Well):
            home_catch = self.get_home(home_info)

        elif isinstance(home_info, gpd.GeoDataFrame):
            home_catch = home_info.copy()

        else:
            sys.exit('need to pass geodataframe, lat/long, or Well to get_neighbors()')
        
        neighbors_df = gpd.sjoin(self.catchment_df, home_catch)
        # neighbors_df.drop(columns='index_right', inplace=True)

        return neighbors_df
    
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
            well_geometries.append(WellGeometry(n, 
                                                w_df,
                                                home_df,
                                                neighbors_df,
                                                target_streams,
                                                nearest_df))
        return well_geometries

        
class WellGeometry():
    ''' Simple object with information about a well- the catchment it
        is located in.  The neighboring catchments.  The streams in those
        catchments. And a table of the closest points on the streams to the 
        well, distances, and apportionments.

        If the geoprocessor is set up to pass transmissivity, storativity,
        and streambed conductance, these are added to the object.
    '''

    def __init__(self, name, well_df, home_df, neighbors_df,
                streams_df, close_points_df, transmissivity=None, storativity=None,
                streambed_cond=None, aq_type=None):
        '''
            Parameters
            ----------
            name: string
                name or ID for well
            well_df: geodataframe
                geodataframe with one point for the well of interest
            home_df: geodataframe
                catchment containing the well
            neighbors_df: geodataframe
                neighboring catchments 
            streams_df: geodataframe
                stream linework in the home and neighboring catchments
            close_points_df: geodataframe
                closest points and apportionments

            Some calls to WellGeometry may have this information for the object,
            but the default is None

            transmissivity: float,
                transmissivity of catchment from shapefile
            storativity: float,
                transmissivity of catchment using default for glacial or bedrock
            streambed_cond: float,
                streambed conductance (lambda) from Kv_w/depth of well screen
            aq_type: string
                typically bedrock or drift (or unconsolidated, etc)

            Returns
            -------
            None
        '''

        self.name = name
        self.well_df = well_df
        self.home_df = home_df
        self.neighbors_df = neighbors_df
        self.streams_df = streams_df
        self.close_points_df = close_points_df
        self.transmissivity = transmissivity
        self.storativity = storativity
        self.streambed_cond = streambed_cond
        self.aq_type = aq_type


    def summarize(self):
        ''' print out the attributes of the WellGeometry object
        
        '''
        print(f'Name = {self.name}')
        print(f'Aquifer type = {self.aq_type}')
        print(f'Transmissivity = {self.transmissivity}')
        print(f'Storativity = {self.storativity}')
        print(f'Streambed conductance = {self.streambed_cond}')
        print(f'Home Catchment')
        print(self.home_df[['BASIN', 'lat', 'long', 'rate', 'depth']].to_markdown())
        print()
        print(f'Neighboring Catchments')
        print(self.neighbors_df[['BASIN_left', 'TYPE_left']].to_markdown())
        print()
        print(f'Streams')
        print(self.streams_df[['BASIN', 'SUBBASIN', 'TYPE']].to_markdown())
        print()
        print(f'Closest Points')
        print(self.close_points_df[['distance', 'inv_distance', 'inv_dist2', 'apportionment', 'apport2']].to_markdown())




