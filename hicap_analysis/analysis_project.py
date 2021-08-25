from hicap_analysis.wells import Well, WellResponse
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

            

            
            



    def aggregate_responses(self):
        # identify which responses and which wells eg. all of a certain status, or all,
        # or a single well
        pass



    def calculate_significance(self):
        """[summary]
        """
        pass


