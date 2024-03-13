import sys
sys.path.append('../')
import hicap_analysis.wells
import hicap_analysis.analysis_project as ap
from hicap_analysis.analysis_project import Project
from pathlib import Path
import numpy as np
import pandas as pd


# ### Path to the folder and the name of the .yml file you've created.
# ##### Note use of double backslashes in the path name
datapath = Path('./MASTER')
odir = datapath / 'output'
yml_file = 'TestExample.yml'#sys.argv[1]
# ts_path = sys.argv[2]
# only read out year 5 for time series ::: hard coded
times = range(365*4,365*5+1)


bdplobs = pd.read_csv(datapath / 'basedeplobs.dat', header=None)
bdplobs.columns = ['obsname']
bdplobs.index = bdplobs.obsname
bdplobs['obs_values'] = np.nan


ap = Project(datapath/yml_file)
ap.report_responses()
ap.write_responses_csv()




base_data = pd.read_csv(
    odir/f'{yml_file.replace(".yml","")}.table_report.base_stream_depletion.csv', index_col=0)

for cob in bdplobs.obsname:
    riv,wel,_ = cob.split(':')
    print(cob)
    bdplobs.loc[cob, 'obs_values'] = base_data.loc[wel][riv]

output_ts = [i.strip() for i in open(ts_path, 'r').readlines()]

ts_obs = []
for c_ts in output_ts:
    ts_obs.extend([f'{c_ts}__{i}' for i in times])
    
ts_df = pd.DataFrame(index = ts_obs, data = {'obsname':ts_obs,'obs_values':np.nan})
ts_data = pd.read_csv(
    odir/f'{yml_file.replace(".yml","")}.table_report.all_ts.csv', index_col=0) 
for cob in ts_df.index:
    criv,ctime = cob.split('__')
    ts_df.loc[cob,'obs_values'] = ts_data.loc[int(ctime)][criv]


allout = pd.concat([bdplobs,ts_df])
# write out all the observations
allout['obs_values'].to_csv('allobs.out', sep = ' ', header=None)