from hicap_analysis.wells import GPM2CFD, Well, WellResponse
from hicap_analysis.utilities import Q2ts
import numpy as np
import pandas as pd
import yaml, os, shutil
from math import radians, cos, sin, asin, sqrt

# def _loc_to_dist(loc0, loc1):
#     '''
#     Euclidean distance between two 2-d points assuming expressed in consistent units
#     '''
#     dist = np.sqrt((loc0[0]-loc1[0])**2 +(loc0[1]-loc1[1])**2)
#     return dist


def _loc_to_dist(loc0, loc1):
    '''
    Distance between two points in lat/long using Haversine Formula assuming lat/long in decimal degrees, returned in feet
    '''
    #convert decimal degrees to radians
    lon1, lat1, lon2, lat2 = map(radians, [loc0[0], loc0[1], loc1[0],loc1[1]])
    #haversine
    dlon = abs(lon2 - lon1)
    dlat = abs(lat2 - lat1)
    a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
    c = 2 * asin(sqrt(a))
    r = 3956*5280 # radius of the earth in feet
    dist = c*r
    return dist

def _print_to_screen_and_file(s, ofp):
    """[summary]

    Args:
        s ([type]): [description]
        ofp ([type]): [description]
    """
    ofp.write(f'{s}\n')
    print(s)

# some helper functions for printing out results
def _print_dd_depl(ofp, cw_dd, cw_max_depl):
    ofp.write(f'             **Drawdown**\n{"Response":30s}{"Drawdown(ft)":30s}\n')
    for ck,v in cw_dd.items():
        ofp.write(f'{ck:30s}{v:<30.4f}\n')
    ofp.write(f'          **Maximum Depletion**\n{"Response":30s}{"Depletion(cfs)":30s}\n')
    for ck,v in cw_max_depl.items():
        ofp.write(f'{ck:30s}{v:<30.4f}\n')

def _print_single_well_header(ofp, wname, wstatus):
    ofp.write('#'*50 + '\n')
    ofp.write(f'Well Name: {wname}\n')
    ofp.write(f'Well status: {wstatus}\n')
def _print_combined_well_results(ofp, cw_dat):
    ofp.write('#'*50 + '\n')
    total_dd = np.sum([v for _,v in cw_dat.drawdown.items()])
    ofp.write('Total Drawdown (ft):   {total_dd:<16.4f}')
class Project():
    def __init__(self, ymlfile):
        """[summary]
        """
        self.status_categories = ['existing', 'active', 'pending', 'new_approved', 'inactive']
        self.wells = {} # dictionary to hold well objects
        self.defaults = ['dd_days','depletion_years','pumping_days'] # allowable project default parameters
        self.stream_apportionment_dict = {}
        self.proposed_well_categories =['pending']
        self.existing_well_categories = ['existing', 'active',  'new_approved']
        self.existing_wells = []
        self.proposed_wells = []
        self.depl_method = 'Walton' # default, can specify in the yml file
        self.__dd_responses = None
        self.__stream_responses = None
    # def populate_from_yaml(self, ymlfile):
    #     """[summary]

    #     Args:
    #         ymlfile ([Path or string]): Configuration file (YAML style) for a project
    #     """
        self.ymlfile = ymlfile
        with open(ymlfile) as ifp:
            d = yaml.safe_load(ifp)


        # parse project_properties block
        if 'project_properties' in d.keys():
            self._parse_project_properties(d['project_properties'])
        else:
            raise('Configuration YAML file must have a "project_properties" block')

        # get the keys for all the remaining blocks
        self.wellkeys = [i for i in d.keys() if i.lower().startswith('well')]
        self.ddkeys = [i for i in d.keys() if i.lower().startswith('dd_resp')]
        self.streamkeys = [i for i in d.keys() if i.lower().startswith('stream')]

        # look for a timeseries file in the project_properties block to determine how to 
        # handle pumping 
        if 'pumping_timeseries_file' in d['project_properties'].keys():
            self.tsfile = d['project_properties']['pumping_timeseries_file']
            self.ts = True
            self.Q_ts = pd.read_csv(self.tsfile).set_index('sequential_day')
            # first test, if there is a time series file that all well keys are columns
            if self.ts is True:
                try:
                    assert all([i  in self.Q_ts.columns for i in self.wellkeys])
                except:
                    raise AssertionError('not all well names are represented in the time series file') 
                
        else:
            self.ts = False
            self.Q_ts = None

        # parse stream responses blocks
        if len(self.streamkeys)>0:
             self._parse_responses(self.streamkeys, d)
        else:
            print('no stream responses supplied for evaluation ')

        # parse drawdown responses blocks
        if len(self.ddkeys)>0:
            self._parse_responses(self.ddkeys, d)
        else:
            print('no drawdown responses supplied for evaluation ')

        # parse well blocks
        if len(self.wellkeys)>0:
            self._parse_wells(d, self.ts, self.Q_ts)
        else:
            raise('No wells were defined in the input file. Goodbye')

        #TODO: verify that stream apportionment and stream response are same keys for a well
        #TODO: verify that all responses called out in wells exist in yaml file

        # create well objects
        self._create_well_objects()

        # report out on yaml input to screen and logfile
        self._report_yaml_input()

    def _parse_project_properties(self, pp):
        """Method to parse all the project properties from the YAML file block

        Args:
            pp ([dict]): project properties block read from YML]
        """
        if 'depl_method' in pp.keys():
            self.depl_method = pp['depl_method']
        try:
            self.name = pp['name']
            self.T = pp['T']
            self.S = pp['S']
            self.default_parameters = {}
            self.default_parameters['default_dd_days'] = pp['default_dd_days']
            self.default_parameters['default_depletion_years'] = pp['default_depletion_years']
            self.default_parameters['default_pumping_days'] = pp['default_pumping_days']
        except:
            raise('Formatting problem with "project_properties" block')

    def _parse_responses(self, keys, d):
        """populate information about the responses to pull from in calculations

        Args:
            keys ([type]): [description]
            d ([dict]): yml file data
        """

        if keys[0].lower().startswith('dd'):
            self.__dd_responses = {}
            cr = self.__dd_responses    
        elif keys[0].lower().startswith('stream'):
            self.__stream_responses = {}
            cr = self.__stream_responses
        # populate the appropriate dictionary on self with location and name
        # information of response
        for ck in keys:
            cr[d[ck]['name']] = d[ck]['loc']
            if 'streambed_conductance' in d[ck].keys():
                cr[d[ck]['name']]['streambed_conductance'] = d[ck]['streambed_conductance']

    def _parse_wells(self, d, ts, Q_ts):
        """populate information about wells assigning apportionment values and
            the lists of proposed and existing wells

        Args:
            d (dict): yml file data
            ts (bool): flag as to whether a timeseries dataframe was read in
            Q_ts (pandas DataFrame): pumping data table for all the wells
        """
        self.__well_data = {}
        
        for ck in self.wellkeys:
            # populate dictionary using well name as key with all well data
            self.__well_data[d[ck]['name']] = d[ck]
            
            # make sure if ts is supplied that Q is not supplied for each well
            if ts is True:
                if 'Q' in d[ck].keys() or 'pumping_days' in d[ck].keys():
                    raise('ERROR:\ntime series file was supplied AND Q of pumping_days was supplied for at\n' +
                        'one well. User can only supply pumping rates in one or the other\n' +
                        'Please try again....')
            
            # populate categories of wells
            if d[ck]['status'] in self.proposed_well_categories:
                self.proposed_wells.append(d[ck]['name'])
            elif d[ck]['status'] in self.existing_well_categories:
                self.existing_wells.append(d[ck]['name'])
            # also, parse out stream apportionment for existing and proposed wells
            streamappkeys = [i for i in d[ck].keys() if 'apportion' in i]
            if len(streamappkeys) > 0:
                self.stream_apportionment_dict[d[ck]['name']] = {}
                for cak in streamappkeys:
                    self.stream_apportionment_dict[d[ck]['name']][d[ck][cak]['name']] = d[ck][cak]['apportionment']

    def _create_well_objects(self):
        """Prepare to populate a Well object for each well, using the attributes of each well and response
        """

        for ck, cw in self._Project__well_data.items():

            # update defaults as appropriate
            for currdef in self.defaults:
                if currdef not in cw.keys():
                    cw[currdef] = self.default_parameters[f'default_{currdef}']
            # calculate all necessary distances
            # first streams
            stream_dist = None
            if 'stream_response' in cw.keys():   
                stream_dist = {}
                streambed_conductance = {}
                streambed_cond_calc = 0
                for c_resp in cw['stream_response']:
                    streamx = self._Project__stream_responses[c_resp]['x']
                    streamy = self._Project__stream_responses[c_resp]['y'] 
                    stream_dist[c_resp] = _loc_to_dist([cw['loc']['x'],cw['loc']['y']], [streamx, streamy])
                    if 'streambed_conductance' in self._Project__stream_responses[c_resp].keys():
                        streambed_conductance[c_resp] = self._Project__stream_responses[c_resp]['streambed_conductance']
                        streambed_cond_calc += 1
                if streambed_cond_calc < 1:
                    streambed_conductance = None


            # next, drawdowns
            dd_dist = None
            if 'dd_response' in cw.keys():   
                dd_dist = {}
                for c_resp in cw['dd_response']:
                    ddx = self._Project__dd_responses[c_resp]['x']
                    ddy = self._Project__dd_responses[c_resp]['y'] 
                    dd_dist[c_resp] = _loc_to_dist([cw['loc']['x'],cw['loc']['y']], [ddx, ddy])
            
            if ck in self.stream_apportionment_dict.keys():
                stream_app_d =self.stream_apportionment_dict[ck]
            else:
                stream_app_d = None
            
            # sort out the time series for wells and convert to CFD
            if self.ts is True:
                Q = self.Q_ts[ck] * GPM2CFD
            else:
                Q = Q2ts(cw['pumping_days'], cw['depletion_years'], cw['Q'])

            self.wells[ck] = Well(T=self.T, S=self.S, Q=Q, 
                    theis_dd_days=cw['dd_days'], depletion_years=cw['depletion_years'],
                    stream_dist=stream_dist, drawdown_dist=dd_dist,
                    stream_apportionment=stream_app_d, depl_method = self.depl_method,
                    streambed_conductance=streambed_conductance
            )


    def _report_yaml_input(self):
        """summarize broad details of the YAML file read in
        """
        logfile = str(self.ymlfile).replace('.yml','.yml.import_report')
        logfile = logfile.replace('.yaml','.yml.import_report')
        with open(logfile, 'w') as ofp:
            print(f'Writing report to {logfile}\n\n')
            _print_to_screen_and_file('',ofp)
            _print_to_screen_and_file(f'Successfully parsed {self.ymlfile} (high five!)',ofp)
            _print_to_screen_and_file('*'*25,ofp)
            _print_to_screen_and_file('Summary follows:',ofp)
            _print_to_screen_and_file('\nWELLS:',ofp)
            # print well summary
            for ck in self.status_categories:
                cw = [k for k,v in self._Project__well_data.items() if v['status']==ck]
                if len(cw) >0:
                    _print_to_screen_and_file(f'{len(cw)} {ck} wells:',ofp)
                    [_print_to_screen_and_file(f'\t{i}',ofp) for i in cw]
            # stream response summary
            if self._Project__stream_responses is not None:
                _print_to_screen_and_file('\nSTREAM RESPONSES:',ofp)
                [_print_to_screen_and_file(f'\t{ck}', ofp) for ck in self._Project__stream_responses.keys()]
            else:
                _print_to_screen_and_file('No Stream Responses in the yml file',ofp)

            # drawdown response summary
            if self._Project__dd_responses is not None:
                _print_to_screen_and_file('DRAWDOWN RESPONSES:',ofp)
                [_print_to_screen_and_file(f'\t{ck}', ofp) for ck in self._Project__dd_responses.keys()]
            else:
                _print_to_screen_and_file('No Drawdown Responses in the yml file',ofp)


    def report_responses(self):
        # make a report file - named from the YML name
        ymlbase = self.ymlfile.name
        outfile = ymlbase.replace('.yml','.report.txt')
        # make a home for the report file
        outpath = self.ymlfile.parent / 'output'
        if not os.path.exists(outpath):
            os.mkdir(outpath)
        self.report_filename = outpath  / outfile
        with open(self.report_filename, 'w') as ofp:
            ofp.write(f'HiCap well analysis report, configured from: {ymlbase}\n')

            # Report on each well
            ofp.write('\nINDIVIDUAL PROPOSED WELL REPORTS\n')
            if len(self.proposed_wells) == 0:
                ofp.write('  there were no proposed wells in the configuration file!\n')
            else:
                for cex in self.proposed_wells:
                    _print_single_well_header(ofp, cex, self._Project__well_data[cex]["status"])
                    cw_dat = self.wells[cex]
                    _print_dd_depl(ofp, cw_dat.drawdown, cw_dat.max_depletion)
                ofp.write('\n')

            ofp.write('\n\nINDIVIDUAL EXISTING WELL REPORTS\n')
            if len(self.existing_wells) == 0:
                ofp.write('  there were no existing wells in the configuration file!\n')
            else:
                for cex in self.existing_wells:
                    _print_single_well_header(ofp, cex, self._Project__well_data[cex]["status"])
                    cw_dat = self.wells[cex]
                    _print_dd_depl(ofp, cw_dat.drawdown, cw_dat.max_depletion)
                ofp.write('\n')

            self.aggregate_results()

            ofp.write('\n\nCOMBINED PROPOSED WELL REPORTS\n' + '#'*50 + '\n')
            if len(self.proposed_wells) == 0:
                ofp.write('  there were no proposed wells in the configuration file!\n')
            else:
                _print_dd_depl(ofp, self.proposed_aggregated_drawdown, self.proposed_aggregated_max_depletion)
                    
            ofp.write('\n\nCOMBINED EXISTING WELL REPORTS\n' + '#'*50 + '\n')
            if len(self.existing_wells) == 0:
                ofp.write('  there were no existing wells in the configuration file!\n')
            else:
                _print_dd_depl(ofp, self.existing_aggregated_drawdown, self.existing_aggregated_max_depletion)
                    
            ofp.write('\n\nTOTAL COMBINED WELL REPORTS\n' + '#'*50 + '\n')
            _print_dd_depl(ofp, self.total_aggregated_drawdown, self.total_aggregated_max_depletion)
            

    def aggregate_results(self):
        # make dictionaries to contain the drawdown results
        self.existing_aggregated_drawdown = {}
        self.proposed_aggregated_drawdown = {}
        self.total_aggregated_drawdown = {}
        
        # make dictionaries to contain the max_depletion results
        self.existing_aggregated_max_depletion = {}
        self.proposed_aggregated_max_depletion = {}
        self.total_aggregated_max_depletion = {}
        
        # make dictionaries to contain the sum_depletion results
        self.existing_aggregated_sum_depletion = {}
        self.proposed_aggregated_sum_depletion = {}
        self.total_aggregated_sum_depletion = {}

        # make dictionaries to contain the sum of depletion for each base stream name
        self.existing_aggregated_base_stream_sum_depletion = {}
        self.proposed_aggregated_base_stream_sum_depletion = {}
        self.total_aggregated_base_stream_sum_depletion = {}        
        self.base_streams = []
        # first existing wells
        for cwell in self.existing_wells:
            cw_dd = self.wells[cwell].drawdown
            for ck, v in cw_dd.items():
                if ck not in self.existing_aggregated_drawdown.keys():
                    self.existing_aggregated_drawdown[ck] = v
                else:
                    self.existing_aggregated_drawdown[ck] += v
            # first sum up depletion time series per well to later get max of sum by location
            cw_dep = self.wells[cwell].depletion
            for ck, v in cw_dep.items():
                if ck not in self.existing_aggregated_sum_depletion.keys():
                    self.existing_aggregated_sum_depletion[ck] = v
                else:
                    self.existing_aggregated_sum_depletion[ck] += v
                # also parse to base stream name
                base_key = ck.split(':')[0]
                if base_key not in self.base_streams:
                    self.base_streams.append(base_key)
                if base_key not in self.existing_aggregated_base_stream_sum_depletion.keys():
                    self.existing_aggregated_base_stream_sum_depletion[base_key] = v
                else:
                    self.existing_aggregated_base_stream_sum_depletion[base_key] += v
            
            # then oldskool max per location            
            cw_max_dep = self.wells[cwell].max_depletion
            for ck, v in cw_max_dep.items():
                if ck not in self.existing_aggregated_max_depletion.keys():
                    self.existing_aggregated_max_depletion[ck] = v
                else:
                    self.existing_aggregated_max_depletion[ck] += v
                


        # next proposed wells
        for cwell in self.proposed_wells:
            cw_dd = self.wells[cwell].drawdown
            for ck, v in cw_dd.items():
                if ck not in self.proposed_aggregated_drawdown.keys():
                    self.proposed_aggregated_drawdown[ck] = v
                else:
                    self.proposed_aggregated_drawdown[ck] += v
            # first sum up depletion time series per well to later get max of sum by location
            cw_dep = self.wells[cwell].depletion
            for ck, v in cw_dep.items():
                if ck not in self.proposed_aggregated_sum_depletion.keys():
                    self.proposed_aggregated_sum_depletion[ck] = v
                else:
                    self.proposed_aggregated_sum_depletion[ck] += v
            # also parse to base stream name
                base_key = ck.split(':')[0]
                if base_key not in self.proposed_aggregated_base_stream_sum_depletion.keys():
                    self.proposed_aggregated_base_stream_sum_depletion[base_key] = v
                else:
                    self.proposed_aggregated_base_stream_sum_depletion[base_key] += v
            # then oldskool max per location            
            cw_max_dep = self.wells[cwell].max_depletion
            for ck, v in cw_max_dep.items():
                if ck not in self.proposed_aggregated_max_depletion.keys():
                    self.proposed_aggregated_max_depletion[ck] = v
                else:
                    self.proposed_aggregated_max_depletion[ck] += v

        # finally calculate the totals
        for cwell in self.existing_wells+self.proposed_wells:
            cw_dd = self.wells[cwell].drawdown
            for ck, v in cw_dd.items():
                if ck not in self.total_aggregated_drawdown.keys():
                    self.total_aggregated_drawdown[ck] = v
                else:
                    self.total_aggregated_drawdown[ck] += v
            # first sum up depletion time series per well to later get max of sum by location
            cw_dep = self.wells[cwell].depletion
            for ck, v in cw_dep.items():
                if ck not in self.total_aggregated_sum_depletion.keys():
                    self.total_aggregated_sum_depletion[ck] = v
                else:
                    self.total_aggregated_sum_depletion[ck] += v
                # also parse to base stream name
                base_key = ck.split(':')[0]
                if base_key not in self.total_aggregated_base_stream_sum_depletion.keys():
                    self.total_aggregated_base_stream_sum_depletion[base_key] = v
                else:
                    self.total_aggregated_base_stream_sum_depletion[base_key] += v                    
            # then oldskool max per location            
            cw_max_dep = self.wells[cwell].max_depletion
            for ck, v in cw_max_dep.items():
                if ck not in self.total_aggregated_max_depletion.keys():
                    self.total_aggregated_max_depletion[ck] = v
                else:
                    self.total_aggregated_max_depletion[ck] += v

    def write_responses_csv(self):
        # create a dataframe to hold the aggregated results
        cols = [f'{i}:dd (ft)' for i in self.total_aggregated_drawdown.keys()] + \
                        [f'{i}:depl (cfs)' for i in self.total_aggregated_max_depletion.keys()]
        col_base = [i.replace(':dd (ft)','').replace(':depl (cfs)','') for i in cols]
 
        rows = [f'{i}: proposed' for i in self.proposed_wells] + \
                [f'{i}: existing' for i in self.existing_wells] + \
                ['total_proposed', 'total_existing','total_combined']
        row_base = [i.replace(': existing','').replace(': proposed','') for i in rows]
        agg_df = pd.DataFrame(index=row_base, columns=col_base)

        if self._Project__stream_responses is not None:
            # now make a special case dataframe for aggregated results by base stream name
            agg_base_stream_df  =  pd.DataFrame(index=row_base, columns=self.base_streams)
            all_depl_ts =pd.DataFrame(index=
                self.wells[list(self.wells.keys())[0]].Q.index)

            # fill in the dataframe
            # individual wells
            for cn, cw in self.wells.items():
                for cresp, cdd in cw.drawdown.items():
                    agg_df.loc[cn,cresp] = cdd
                for cresp, cdepl in cw.max_depletion.items():
                    agg_df.loc[cn,cresp] = cdepl
                    basekey = cresp.split(':')[0]
                    agg_base_stream_df.loc[cn,basekey] = cdepl
                all_depl_ts = pd.concat(
                    (all_depl_ts,pd.DataFrame(index=all_depl_ts.index, 
                                            data=cw.depletion)), 
                    axis=1
                    )

            # totals
            #proposed
            for cresp, cdd in self.proposed_aggregated_drawdown.items():
                agg_df.loc['total_proposed', cresp] = cdd
            for cresp, cdepl in self.proposed_aggregated_max_depletion.items():
                agg_df.loc['total_proposed', cresp] = cdepl
            for cresp, cdepl in self.proposed_aggregated_base_stream_sum_depletion.items():
                agg_base_stream_df.loc['total_proposed', cresp] = np.max(cdepl)
                
            #existing
            for cresp, cdd in self.existing_aggregated_drawdown.items():
                agg_df.loc['total_existing', cresp] = cdd
            for cresp, cdepl in self.existing_aggregated_max_depletion.items():
                agg_df.loc['total_existing', cresp] = cdepl
            for cresp, cdepl in self.existing_aggregated_base_stream_sum_depletion.items():
                agg_base_stream_df.loc['total_existing', cresp] = np.max(cdepl)
                
            #total
            for cresp, cdd in self.total_aggregated_drawdown.items():
                agg_df.loc['total_combined', cresp] = cdd
            for cresp, cdepl in self.total_aggregated_max_depletion.items():
                agg_df.loc['total_combined', cresp] = cdepl
            for cresp, cdepl in self.total_aggregated_base_stream_sum_depletion.items():
                agg_base_stream_df.loc['total_combined', cresp] = np.max(cdepl)
                

            agg_df.columns = cols
            agg_df.index = rows
            
            # make a report file - named from the YML name
            ymlbase = self.ymlfile.name
            outfile = ymlbase.replace('.yml','.table_report.csv')
            # make a home for the report file
            outpath = self.ymlfile.parent / 'output'
            if not os.path.exists(outpath):
                os.mkdir(outpath)

            self.csv_output_filename = outpath / outfile
            agg_df.to_csv(self.csv_output_filename)
            # slap the csv dataframes into self
            self.agg_df = agg_df
            self.csv_stream_output_filename = outpath / outfile.replace('.csv','.base_stream_depletion.csv')
            agg_base_stream_df.to_csv(self.csv_stream_output_filename)
            self.csv_stream_output_ts_filename = outpath / outfile.replace('.csv','.all_ts.csv')
            all_depl_ts.to_csv(self.csv_stream_output_ts_filename)

            self.agg_base_stream_df = agg_base_stream_df
        
        else:
            for cn, cw in self.wells.items():
                for cresp, cdd in cw.drawdown.items():
                    agg_df.loc[cn,cresp] = cdd
            # totals
            #proposed
            for cresp, cdd in self.proposed_aggregated_drawdown.items():
                agg_df.loc['total_proposed', cresp] = cdd
            #existing
            for cresp, cdd in self.existing_aggregated_drawdown.items():
                agg_df.loc['total_existing', cresp] = cdd
            #total
            for cresp, cdd in self.total_aggregated_drawdown.items():
                agg_df.loc['total_combined', cresp] = cdd    
                

            agg_df.columns = cols
            agg_df.index = rows
            
            # make a report file - named from the YML name
            ymlbase = self.ymlfile.name
            outfile = ymlbase.replace('.yml','.table_report.csv')
            # make a home for the report file
            outpath = self.ymlfile.parent / 'output'
            if not os.path.exists(outpath):
                os.mkdir(outpath)

            self.csv_output_filename = outpath / outfile
            agg_df.to_csv(self.csv_output_filename)
            