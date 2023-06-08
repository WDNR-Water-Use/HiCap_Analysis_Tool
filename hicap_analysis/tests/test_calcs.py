import sys
sys.path.append(r"C:\Users\doughwd\Documents\Water Use\VSCODE\HiCap_Analysis_Tool")
from hicap_analysis.wells import GPM2CFD
from os import pardir
import numpy as np
import pandas as pd
from pathlib import Path
import pytest


datapath = Path('hicap_analysis/tests/data')

@pytest.fixture
def theis_results():
    from hicap_analysis import wells as wo
    excel_file = datapath / 'HighCap_Analysis_Worksheet_Example.xlsm'
    p = pd.read_excel(excel_file,
        sheet_name='Property_Drawdown_Analysis',
        usecols='C:D',
        skiprows=7,
        index_col=0)
    # read the two Q values and convert to CFD
    Q = [float(i)*wo.GPM2CFD for i in [p.loc['Pumping Rate Well #1 (gpm)'],
            p.loc['Pumping Rate Well #2 (gpm)']]]
    S = float(p.loc['Storage Coefficient (unitless)'])
    T = float(p.loc['Transmissivity (ft2/day)']) 
    time = float(p.loc['Theis Time of Pumping (days)'])
    params = {'Q':Q,'S':S,'T':T, 'time':time}
    theis_res = pd.read_excel(excel_file, 
        sheet_name='Theis', 
        skiprows= 9,
        usecols=('A:B,H:I'),
        names=['well1_dd','well1_r','well2_dd','well2_r'])
      
    return {'params':params, 'theis_res':theis_res}


@pytest.fixture
def walton_results():
    from hicap_analysis import wells as wo
    excel_file = datapath / 'HighCap_Analysis_Worksheet_Example.xlsm'  
    walton_res = pd.read_excel(excel_file,
        sheet_name='Stream#1_Depletion', 
        skiprows= 104,
        usecols=('C,M:N,R,AB:AC,AK'),
        names=['t_well','dep1','dep2','t_image','rch1','rch2', 'total_dep'])
    p = pd.read_excel(excel_file,
        sheet_name='Stream#1_Depletion', 
        skiprows= 70,
        nrows=30,
        usecols=('B:D'),
        names=['par','v1','v2'],
        index_col=0,)
    Q = p.loc['Q - Pumping rate (ft^3/dy)'].T.values
    S = p.loc['S - Storage (unitless)'].T.values
    dist = p.loc['a - Distance (feet)'].T.values
    T_gpd_ft = p.loc['T - Transmissivity (gpd/ft)'].T.values
    # a little trickery to get the index of the time array for start and end of each well
    Q_start_day = [pd.to_datetime(i).day_of_year - 1 
        for i in p.loc['  First Day of Annual Pumping ='].T.values]
    Q_end_day = [pd.to_datetime(i).day_of_year - 1 
        for i in p.loc['  Last Day of Annual Pumping ='].T.values]
    params = {'Q':Q, 'S':S, 'T_gpd_ft':T_gpd_ft,
            'Q_start_day':Q_start_day,
            'Q_end_day':Q_end_day, 'dist':dist}
  
  
    return {'params':params, 'walton_res':walton_res}

@pytest.fixture
def project_spreadsheet_results():
    from hicap_analysis import wells as wo
    excel_file = datapath / 'HighCap_Analysis_Worksheet_Example.xlsm'  
    # read in common parameters
    p = pd.read_excel(excel_file, sheet_name='Property_Drawdown_Analysis', 
                skiprows= 7, nrows=12, usecols=('C:D'), index_col=0)
    p1 = pd.read_excel(excel_file, sheet_name='Property_Drawdown_Analysis', 
                skiprows= 19, nrows=12, usecols=('C:D'), index_col=0)
    p2 = pd.read_excel(excel_file, sheet_name='Property_Drawdown_Analysis', 
                skiprows= 31, nrows=12, usecols=('C:D'), index_col=0)
    p3 = pd.read_excel(excel_file, sheet_name='Property_Drawdown_Analysis', 
                skiprows= 57, nrows=50, usecols=('C:F'), index_col=0)
    p4 = pd.read_excel(excel_file, sheet_name='Cumulative_Impact_Analysis', 
                skiprows= 22, nrows=5, usecols=('C:D'), index_col=0)
    p5 = pd.read_excel(excel_file, sheet_name='Cumulative_Impact_Analysis', 
                skiprows= 36, nrows=10, usecols=('H:I'), index_col=0)
                
    params= {'T': p.loc['Transmissivity (ft2/day)'].values[0],
                    'S': p.loc['Storage Coefficient (unitless)'].values[0],
                    'Q1_gpm': p.loc['Pumping Rate Well #1 (gpm)'].values[0],
                    'Q2_gpm': p.loc['Pumping Rate Well #2 (gpm)'].values[0],
                    'w1muni_dist': p3.loc['Distance from Well #1 to Municpal Well'].values[0],
                    'w2muni_dist': p3.loc['Distance from Well #2 to Municpal Well'].values[0], 
                    'w1sprng1_dist': p3.loc['Distance from Well #1 to Spring'].values[0],
                    'w2sprng1_dist': p3.loc['Distance from Well #2 to Spring'].values[0],
                    'muni_dd_combined_proposed': p3.loc['Distance from Well #1 to Municpal Well'].values[-1],
                    'sprng1_dd_combined_proposed':p3.loc['Distance from Well #1 to Spring'].values[-1],
                    'well1_5ftdd_loc': p3.loc[' Well #1 5-ft Drawdown (feet)'].values[0],
                    'well1_1ftdd_loc': p3.loc[' Well #1 1-ft Drawdown (feet)'].values[0],
                    'theis_p_time': p.loc['Theis Time of Pumping (days)'].values[0],
                    'stream_name_1': p1.loc['Stream Name'].values[0],
                    'stream_name_2': p2.loc['Stream Name'].values[0],
                    'depl_pump_time':p1.loc['Stream Depletion Duration Period (Days)'].values[0],
                    'w1s1_dist': p1.loc['Well #1 - Distance to Stream (feet)'].values[0],     
                    'w1s1_appor': p1.loc['Well #1 - Fraction Intercepting Stream (.1-1)'].values[0],
                    'w2s1_dist': p1.loc['Well #2 - Distance to Stream (feet)'].values[0],
                    'w2s1_appor': p1.loc['Well #2 - Fraction Intercepting Stream (.1-1)'].values[0],
                    'w1s2_dist': p2.loc['Well #1 - Distance to Stream (feet)'].values[0],
                    'w1s2_appor': p2.loc['Well #1 - Fraction Intercepting Stream (.1-1)'].values[0],
                    'w2s2_dist': p2.loc['Well #2 - Distance to Stream (feet)'].values[0],
                    'w2s2_appor': p2.loc['Well #2 - Fraction Intercepting Stream (.1-1)'].values[0],
                    's1_4yr_depl_cfs': p3.loc['Stream #1 depletion after year 4 (cfs)'].values[0],
                    's2_4yr_depl_cfs': p3.loc['Stream #2 depletion after year 4  (cfs)'].values[0],
                    'muni_dd_total_combined': p4.loc['Cumulative Impact Drawdown (ft)'].values[0],
                    'stream1_depl_existing':p5.iloc[0].values[0],
                    'stream1_depl_total_combined':p5.iloc[3].values[0]
                     }            
    return params

def test_project_spreadsheet(project_spreadsheet_results):
    from hicap_analysis.wells import Well, GPM2CFD
    pars = project_spreadsheet_results
    # set up the Project with multiple wells and multiple streams and make calculations
    well1 = Well(T=pars['T'], S=pars['S'], Q=pars['Q1_gpm']*GPM2CFD, depletion_years=5,
                theis_dd_days=pars['theis_p_time'],depl_pump_time=pars['depl_pump_time'],
                stream_dist = {pars['stream_name_1']:pars['w1s1_dist'], pars['stream_name_2']:pars['w1s2_dist']},
                drawdown_dist={'muni':pars['w1muni_dist']},
                stream_apportionment={pars['stream_name_1']:pars['w1s1_appor'],pars['stream_name_2']:pars['w1s2_appor']})
    well2 = Well(T=pars['T'], S=pars['S'], Q=pars['Q2_gpm']*GPM2CFD, depletion_years=5,
                theis_dd_days=pars['theis_p_time'],depl_pump_time=pars['depl_pump_time'],
                stream_dist = {pars['stream_name_1']:pars['w2s1_dist'], pars['stream_name_2']:pars['w2s2_dist']},
                drawdown_dist={'muni':pars['w2muni_dist']},
                stream_apportionment={pars['stream_name_1']:pars['w2s1_appor'],pars['stream_name_2']:pars['w2s2_appor']})
    dd1 = well1.drawdown['muni']
    dd2 = well2.drawdown['muni']
 
    assert np.allclose(dd1+dd2, pars['muni_dd_combined_proposed'], atol=0.1)
    

    depl1 = well1.depletion
    depl2 = well2.depletion
    stream1_max_depl = np.max(depl1[pars['stream_name_1']]) + np.max(depl2[pars['stream_name_1']])
    stream2_max_depl = np.max(depl1[pars['stream_name_2']]) + np.max(depl2[pars['stream_name_2']])
    assert np.allclose(stream1_max_depl, pars['s1_4yr_depl_cfs'], atol=1e-2)
    assert np.allclose(stream2_max_depl, pars['s2_4yr_depl_cfs'], atol=1e-2)
    
def test_theis(theis_results):
    """Test for the theis calculations - compared with two wells at multiple distances
        in the example spreadsheet

    Args:
        theis_results (@fixture, dict): parameters and results from example spreadsheet
    """


    from hicap_analysis import wells as wo
    pars = theis_results['params']
    dist = theis_results['theis_res'].well1_r
    
    time = pars['time']
    dd = [wo._theis(pars['T'], pars['S'], time, dist, currQ) for currQ in pars['Q']]
    assert np.allclose(dd[0],theis_results['theis_res'].well1_dd, atol=0.5)
    assert np.allclose(dd[1],theis_results['theis_res'].well2_dd, atol=0.7)

def test_distance():
    from hicap_analysis import analysis_project as ap
    assert np.isclose(ap._loc_to_dist([89.38323, 43.07476],[89.38492, 43.07479]), 450.09, atol=0.1)
  #  ([2,3],[9,32.9]), 30.70846788753877)

def test_glover():
    """Test for the glover calculations
        against the Glover & Balmer (1954) paper
    """
    from hicap_analysis import wells as wo
    dist = [1000, 5000, 10000]
    Q = 1 * 3600 * 24 # no normalization in the paper but use to convert from CFS to CFD
    time = 365 * 5 # paper evaluates at 5 years in days
    K = 0.001 # ft/sec
    D = 100 # thickness in feet
    T = K*D*24*60*60 # converting to ft/day
    S = 0.2
    Qs = wo._glover(T,S,time, dist, Q)
    assert all(np.isnan(Qs)== False)
    assert np.allclose(Qs, [0.9365, 0.6906, 0.4259], atol=1e-3)
    
def test_sdf():
    """Test for streamflow depletion factor
        using values from original Jenkins (1968) paper
        https://doi.org/10.1111/j.1745-6584.1968.tb01641.x
        note Jenkins rounded to nearest 10 (page 42)
    """
    from hicap_analysis import wells as wo
    dist = 5280./2.
    T = 5.0e4/7.48
    S = 0.5
    sdf = wo._sdf(T,S,dist)
    assert np.allclose(sdf, 520, atol=1.5)

def test_walton(walton_results):
    """Test of a single year to be sure the Walton calculations are made correctly

    Args:
        walton_results ([type]): [description]
    """
    from hicap_analysis import wells as wo

    res = walton_results['walton_res']
    pars = walton_results['params']
    
    dep={}
    rch={}
    for idx in [0,1]:
        dep[idx] = wo._walton(pars['T_gpd_ft'][idx],
                    pars['S'][idx],
                    pars['dist'][idx],
                    res.t_well,
                    pars['Q'][idx]
                    )
        rch[idx] = wo._walton(pars['T_gpd_ft'][idx],
                    pars['S'][idx],
                    pars['dist'][idx],
                    res.t_image,
                    pars['Q'][idx]
                    )
    dep_tot = dep[0]-rch[0] + dep[1]-rch[1]
    assert np.allclose(dep[0], res.dep1)
    assert np.allclose(dep[1], res.dep2)
    assert np.allclose(rch[0], -res.rch1)
    assert np.allclose(rch[1], -res.rch2)
    assert np.allclose(dep_tot, res.total_dep)

def test_yaml_parsing(project_spreadsheet_results):
    pars = project_spreadsheet_results
    from hicap_analysis.analysis_project import Project 
    from hicap_analysis import wells as wo
    ap = Project(datapath / 'example.yml')
   # ap.populate_from_yaml(datapath / 'example.yml')
    #verify that the created well objects are populated with the same values as in the YML file
    assert set(ap.wells.keys()).difference(set(['new1','new2','Existing_CAFO','Existing_Irrig'])) == set()
    assert set(ap._Project__stream_responses.keys()).difference(set(['Upp Creek', 'no paddle'])) == set()
    assert set(ap._Project__dd_responses.keys()).difference(set(['Muni1', 'Sprng1'])) == set()
    
    # spot check some numbers
    assert ap.wells['new1'].T == 35
    assert np.isclose(wo.GPM2CFD * 1000, ap.wells['new2'].Q)
    assert ap.wells['new2'].stream_apportionment['Upp Creek'] == 0.3


    ap.report_responses()
    
    ap.write_responses_csv()

    agg_results = pd.read_csv(ap.csv_output_filename, index_col=0)
    # read in the CSV file and spot check against the spreadsheet output
    assert np.isclose(pars['muni_dd_combined_proposed'], agg_results.loc['total_proposed', 'Muni1:dd (ft)'], atol=0.1)
    assert np.isclose(pars['sprng1_dd_combined_proposed'], agg_results.loc['total_proposed', 'Sprng1:dd (ft)'], atol=0.002)
    assert np.isclose(pars['stream1_depl_existing'], agg_results.loc['total_existing', 'Upp Creek:depl (cfs)'], atol=0.005)
    assert np.isclose(pars['stream1_depl_total_combined'], agg_results.loc['total_combined', 'Upp Creek:depl (cfs)'], atol=0.01)

def test_complex_yml():
    from hicap_analysis.analysis_project import Project 

    ap = Project(datapath / 'example2.yml')
    ap.report_responses()
    ap.write_responses_csv()
    
    df_ts = pd.read_csv(ap.csv_stream_output_ts_filename, index_col=0)
    df_agg = pd.read_csv(ap.csv_stream_output_filename, index_col=0)
    
    df_ts_max = df_ts.max().to_frame()
    df_ts_max.rename(columns={0:'raw'}, inplace=True)
    s_cols_exist = [i for i in df_ts.columns if ("Spring" in i)&('93444' not in i)]
    s_cols_prop = [i for i in df_ts.columns if ("Spring" in i)&('93444' in i)]

    e_cols_exist=[i for i in df_ts.columns if ("EBranch" in i)&('93444' not in i)]
    e_cols_prop=[i for i in df_ts.columns if ("EBranch" in i)&('93444' in i)]

    s_cols_tot =s_cols_exist + s_cols_prop
    e_cols_tot =e_cols_exist + e_cols_prop
    
    df_ts_max['read'] = [df_agg.loc[i.split(':')[1],i.split(':')[0]] for i in df_ts_max.index ]
  
    assert all(np.isclose(df_ts_max.raw, df_ts_max['read']) )
    
    keys = ('SpringBrook:proposed','SpringBrook:existing','SpringBrook:combined',
         'EBranchEauClaire:proposed','EBranchEauClaire:existing','EBranchEauClaire:combined')
    vals = (s_cols_prop, s_cols_exist, s_cols_tot,
           e_cols_prop, e_cols_exist, e_cols_tot)
    for k,v in zip(keys,vals):
        df_agg_val = df_agg.loc[f'total_{k.split(":")[1]}', k.split(':')[0]]
        calc_val = np.max(df_ts[v].sum(axis=1))
        assert np.isclose(df_agg_val, calc_val)
        
    return('stoked')

def test_run_yml_example():
    import hicap_analysis.wells
    import hicap_analysis.analysis_project as ap
    from hicap_analysis.analysis_project import Project

    yml_file = 'example.yml'
    ap = Project(datapath/yml_file)
    ap.report_responses()
    ap.write_responses_csv()