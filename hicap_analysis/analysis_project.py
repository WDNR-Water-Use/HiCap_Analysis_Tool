<<<<<<< HEAD
from hicap_analysis.wells import GPM2CFD, Well, WellResponse
=======
from hicap_analysis.wells import Well, WellResponse
>>>>>>> 6aa23e2 (fixed example.yml conflicts)
import numpy as np
import pandas as pd
import yaml
def _loc_to_dist(loc0, loc1):
    '''
    Euclidean distance between two 2-d points assuming expressed in consistent units
    '''
    dist = np.sqrt((loc0[0]-loc1[0])**2 +(loc0[1]-loc1[1])**2)
    return dist

def _print_to_screen_and_file(s, ofp):
    """[summary]

    Args:
        s ([type]): [description]
        ofp ([type]): [description]
    """
    ofp.write(f'{s}\n')
    print(s)
class Project():
    def __init__(self) -> None:
        """[summary]
<<<<<<< HEAD
        """
        self.status_categories = ['existing', 'active', 'pending', 'new_approved', 'inactive']
        self.wells = {} # dictionary to hold well objects
        self.defaults = ['dd_days','depletion_years','pumping_days'] # allowable project default parameters
        self.stream_apportionment_dict = {}
        

    def populate_from_yaml(self, ymlfile):
        """[summary]

        Args:
            ymlfile ([Path or string]): Configuration file (YAML style) for a project
        """
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
            self._parse_wells(d)
        else:
            raise('No wells were defined in the input file. Goodbye')

        #TODO: verify that stream approptionment and stream response are same keys for a well
        #TODO: verify that all responses called out in wells exist in yaml file

        # create well objects
        self._create_well_objects()

        # report out on yaml input to screen and logfile
        self._report_yaml_input()


        j = 2

    def _parse_project_properties(self, pp):
        """Method to parse all the project properties from the YAML file block

        Args:
            pp ([dict]): project properties block read from YML]
        """
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
        """[summary]

        Args:
            keys ([type]): [description]
            d ([type]): [description]
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

    def _parse_wells(self, d):
        """[summary]

        Args:
            keys ([type]): [description]
            d ([type]): [description]
        """
        self.__well_data = {}
        for ck in self.wellkeys:
            # populate dictionary using well name as key with all well data
            self.__well_data[d[ck]['name']] = d[ck]

            # also, parse out stream apportionment
            streamappkeys = [i for i in d[ck].keys() if 'apportion' in i]
            if len(streamappkeys) > 0:
                self.stream_apportionment_dict[d[ck]['name']] = {}
                for cak in streamappkeys:
                    self.stream_apportionment_dict[d[ck]['name']][d[ck][cak]['name']] = d[ck][cak]['apportionment']

    def _create_well_objects(self):
        """[summary]
        """

        for ck, cw in self._Project__well_data.items():

            # update defaults as appopriate
            for currdef in self.defaults:
                if currdef not in cw.keys():
                    cw[currdef] = self.default_parameters[f'default_{currdef}']
            # calculate all necessary distances
            # first streams
            stream_dist = None
            if 'stream_response' in cw.keys():   
                stream_dist = {}
                for c_resp in cw['stream_response']:
                    streamx = self._Project__stream_responses[c_resp]['x']
                    streamy = self._Project__stream_responses[c_resp]['y'] 
                    stream_dist[c_resp] = _loc_to_dist([cw['loc']['x'],cw['loc']['y']], [streamx, streamy])
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

            self.wells[ck] = Well(T=self.T, S=self.S, Q=cw['Q']*GPM2CFD, depletion_years=cw['depletion_years'],
            theis_dd_days=cw['dd_days'],depl_pump_time=cw['pumping_days'],stream_dist=stream_dist,drawdown_dist=dd_dist,
            stream_apportionment=stream_app_d
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
            if len(self._Project__stream_responses) > 0:
                _print_to_screen_and_file('\nSTREAM RESPONSES:',ofp)
                [_print_to_screen_and_file(f'\t{ck}', ofp) for ck in self._Project__stream_responses.keys()]
            else:
                _print_to_screen_and_file('No Stream Responses in the yml file',ofp)

            # drawdown response summary
            if len(self._Project__dd_responses) > 0:
                _print_to_screen_and_file('DRAWDOWN RESPONSES:',ofp)
                [_print_to_screen_and_file(f'\t{ck}', ofp) for ck in self._Project__dd_responses.keys()]
            else:
                _print_to_screen_and_file('No Drawdown Responses in the yml file',ofp)

=======
        """
        self.status_categories = ['existing', 'active', 'pending', 'new_approved', 'inactive' ]
        pass

    def populate_from_yaml(self, ymlfile):
        """[summary]

        Args:
            ymlfile ([Path or string]): Configuration file (YAML style) for a project
        """
        self.ymlfile = ymlfile
        with open(ymlfile) as ifp:
            d = yaml.safe_load(ifp)

        # parse project_properties block
        if 'project_properties' in d.keys():
            self._parse_project_properties(d['project_properties'])
        else:
            raise('Configuration YAML file must have a "project_properties" block')

        # get the keys for all the remaining blocks
        wellkeys = [i for i in d.keys() if i.lower().startswith('well')]
        ddkeys = [i for i in d.keys() if i.lower().startswith('dd_resp')]
        streamkeys = [i for i in d.keys() if i.lower().startswith('stream')]

        # parse stream responses blocks
        if len(streamkeys)>0:
            self._parse_responses(streamkeys, d)
        else:
            print('no stream responses supplied for evaluation ')

        # parse drawdown responses blocks
        if len(streamkeys)>0:
            self._parse_responses(ddkeys, d)
        else:
            print('no drawdown responses supplied for evaluation ')

        # parse well blocks
        if len(wellkeys)>0:
            self._parse_wells(wellkeys, d)
        else:
            raise('No wells were defined in the input file. Goodbye')

        '''well2 = Well(T=pars['T'], S=pars['S'], Q=pars['Q2_gpm']*GPM2CFD, depletion_years=5,
                theis_dd_time=pars['theis_p_time'],depl_pump_time=pars['depl_pump_time'],
                stream_dist = {pars['stream_name_1']:pars['w2s1_dist'], pars['stream_name_2']:pars['w2s2_dist']},
                drawdown_dist={'muni':pars['w2muni_dist']},
                stream_apportionment={pars['stream_name_1']:pars['w2s1_appor'],pars['stream_name_2']:pars['w2s2_appor']})'''
        # report out on yaml input to screen and logfile
        self._report_yaml_input()
        j = 2

    def _parse_project_properties(self, pp):
        """Method to parse all the project properties from the YAML file block

        Args:
            pp ([dict]): project properties block read from YML]
        """
        try:
            self.name = pp['name']
            self.T = pp['T']
            self.S = pp['S']
            self.default_dd_time = pp['default_dd_time']
            self.default_depletion_time = pp['default_depletion_time']
            self.default_pump_duration = pp['default_pump_duration']
        except:
            raise('Formatting problem with "project_properties" block')

    def _parse_responses(self, keys, d):
        """[summary]

        Args:
            keys ([type]): [description]
            d ([type]): [description]
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

    def _parse_wells(self, keys, d):
        """[summary]

        Args:
            keys ([type]): [description]
            d ([type]): [description]
        """
        self.__well_data = {}
        for ck in keys:
            # populate dictionary using well name as key with all well data
            self.__well_data[d[ck]['name']] = d[ck]

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
            if len(self._Project__stream_responses) > 0:
                _print_to_screen_and_file('\nSTREAM RESPONSES:',ofp)
                [_print_to_screen_and_file(f'\t{ck}', ofp) for ck in self._Project__stream_responses.keys()]
            else:
                _print_to_screen_and_file('No Stream Responses in the yml file',ofp)

            # drawdown response summary
            if len(self._Project__dd_responses) > 0:
                _print_to_screen_and_file('DRAWDOWN RESPONSES:',ofp)
                [_print_to_screen_and_file(f'\t{ck}', ofp) for ck in self._Project__dd_responses.keys()]
            else:
                _print_to_screen_and_file('No Drawdown Responses in the yml file',ofp)

            

            
            


>>>>>>> 6aa23e2 (fixed example.yml conflicts)

    def aggregate_responses(self):
        # identify which responses and which wells eg. all of a certain status, or all,
        # or a single well
        # TODO: consider creating an aggregation container while parsing the YML (e.g. use the unique set of response
        # keys to instantiate off the couch) 
        # other option is to add key if not already in palce as looping through the wells....
        pass



    def calculate_significance(self):
        """[summary]
        """
        pass


