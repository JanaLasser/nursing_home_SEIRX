from scseirx.model_nursing_home import SEIRX_nursing_home
import scseirx.analysis_functions as af
import pandas as pd
import numpy as np
import networkx as nx
from os.path import join


def compose_agents(measures, simulation_params, e_screen_interval,
                   r_screen_interval, e_vaccination_ratio, r_vaccination_ratio):
    '''
    Utility function to compose agent dictionaries as expected by the simulation
    model as input from the dictionary of prevention measures.
    
    Parameters
    ----------
    measures : dictionary
        Dictionary of prevention measures. Needs to include the fields 
        employee_mask and resident_mask. 
    simulation_params : dictionary
        Dictionary of (epidemiological) simulation parameters. Needs to include
        the fields employee_index_probability and resident_index_probability.
    e_screen_interval : integer
        Interval (in days) of the preventive testing screens in the employee
        agent group. Can be [2, 3, 7, None].
    r_screen_interval : integer
        Interval (in days) of the preventive testing screens in the resident
        agent group. Can be [2, 3, 7, None].
    e_vaccination_ratio : float
        Ratio of vaccinated employees.
    r_vaccination_ratio : float
        Ratio of vaccinated residents.
        
    Returns
    -------
    agent_types : dictionary of dictionaries
        Dictionary containing the fields "screening_interval", 
        "index_probability" and "mask" for the agent groups "employee" and
        "resident".
    
    '''
    agent_types = {
        'employee':{
            'screening_interval':e_screen_interval,
            'vaccination_ratio': e_vaccination_ratio,
            'index_probability':simulation_params['employee_index_probability'],
            'mask':measures['employee_mask']},

        'resident':{
            'screening_interval':r_screen_interval,
            'vaccination_ratio': r_vaccination_ratio,
            'index_probability':simulation_params['resident_index_probability'],
            'mask':measures['resident_mask']},
    }
    return agent_types


def run_model(measures, simulation_params, contact_network_src, test_type,
              index_case, e_screen_interval, r_screen_interval,
              e_vaccination_ratio, r_vaccination_ratio, N_steps=500, run=True):
    '''
    Runs a simulation with an SEIRX_nursing_home model 
    (see https://pypi.org/project/scseirx/1.3.0/), given a set of parameters.
    
    Parameters:
    -----------
    measures : dictionary
        Dictionary listing all prevention measures in place for the given
        scenario. Fields that are not specifically included in this dictionary
        will revert to SEIRX_nursing_home defaults.
    simulation_params : dictionary
        Dictionary holding simulation parameters such as "verbosity" and
        "base_transmission_risk". Fields that are not included will revert back
        to SEIRX_nursing_home defaults.
    contact_network_src : string
        Absolute or relative path pointing to the location of the contact
        network used for the simulation. Networks need to be saved in networkx's
        .bz2 format.
    test_type : string
        Test technology used in the preventive screening. Available test
        technologies are listed in the module testing_strategy.py
    index_case : string
        Agent group from which the index case is drawn. Can be "employee" or
        "resident".
    e_screen_interval : integer
        Interval (in days) of the preventive testing screens in the employee
        agent group. Can be [2, 3, 7, None].
    r_screen_interval : integer
        Interval (in days) of the preventive testing screens in the resident
        agent group. Can be [2, 3, 7, None].
    e_vaccination_ratio : float
        Ratio of vaccinated employees.
    r_vaccination_ratio : float
        Ratio of vaccinated residents.
    N_steps : integer
        Number of maximum steps per run. This is a very conservatively chosen 
        value that ensures that an outbreak will always terminate within the 
        allotted time. Most runs are terminated way earlier anyways, as soon as 
        the outbreak is over.
    run : bool
        If run == False, the model will be initialised but not run. This serves
        for debugging purposes and to check the statistics of the freshly 
        instantiated model.
        
    Returns
    -------
    model : SEIRX_nursing_home model instance holding a completed simulation run
        and all associated data.
    '''
    agent_types = compose_agents(measures, simulation_params,
                                 e_screen_interval, r_screen_interval,
                                 e_vaccination_ratio, r_vaccination_ratio)

    # load the contact graph for a single living unit in a nursing home
    G = nx.readwrite.gpickle.read_gpickle(\
            join(contact_network_src, 'interactions_single_quarter.bz2'))

    # initialize the model
    model = SEIRX_nursing_home(G, 
      simulation_params['verbosity'], 
      base_transmission_risk = simulation_params['base_transmission_risk'],
      testing = measures['testing'],
      exposure_duration = simulation_params['exposure_duration'],
      time_until_symptoms = simulation_params['time_until_symptoms'],
      infection_duration = simulation_params['infection_duration'],
      quarantine_duration = measures['quarantine_duration'],
      subclinical_modifier = simulation_params['subclinical_modifier'],
      infection_risk_contact_type_weights = \
                 simulation_params['infection_risk_contact_type_weights'],
      K1_contact_types = measures['K1_contact_types'],
      diagnostic_test_type = measures['diagnostic_test_type'],
      preventive_screening_test_type = test_type,
      follow_up_testing_interval = \
                 measures['follow_up_testing_interval'],
      liberating_testing = measures['liberating_testing'],
      index_case = index_case,
      agent_types = agent_types, 
      age_transmission_risk_discount = \
                simulation_params['age_transmission_discount'],
      age_symptom_modification = simulation_params['age_symptom_modification'],
      mask_filter_efficiency = measures['mask_filter_efficiency'],
      transmission_risk_ventilation_modifier = \
                measures['transmission_risk_ventilation_modifier'],)

    if run:
        # run the model until the outbreak is over
        for i in range(N_steps):
            # break if first outbreak is over
            if len([a for a in model.schedule.agents if \
                (a.exposed == True or a.infectious == True)]) == 0:
                break
            model.step()
        
    return model

def run_ensemble(N_runs, measures, simulation_params, contact_network_src,
                 ensmbl_dst, test_type=None, index_case='employee',
                 e_screen_interval=None, r_screen_interval=None,
                 e_vaccination_ratio=0, r_vaccination_ratio=0):
    '''
    Utility function to run an ensemble of simulations for a given parameter 
    combination.
    
    Parameters:
    ----------
    N_runs : integer
        Number of individual simulation runs in the ensemble.
    measures : dictionary
        Dictionary listing all prevention measures in place for the given
        scenario. Fields that are not specifically included in this dictionary
        will revert to SEIRX_nursing_home defaults.
    simulation_params : dictionary
        Dictionary holding simulation parameters such as "verbosity" and
        "base_transmission_risk". Fields that are not included will revert back
        to SEIRX_nursing_home defaults.
    contact_network_src : string
        Absolute or relative path pointing to the location of the contact
        network used for the simulation. Networks need to be saved in networkx's
        .bz2 format.
    test_type : string
        Test technology used in the preventive screening. Available test
        technologies are listed in the module testing_strategy.py
    index_case : string
        Agent group from which the index case is drawn. Can be "employee" or
        "resident".
    e_screen_interval : integer
        Interval (in days) of the preventive testing screens in the employee
        agent group. Can be [2, 3, 7, None].
    r_screen_interval : integer
        Interval (in days) of the preventive testing screens in the resident
        agent group. Can be [2, 3, 7, None].
    e_vaccination_ratio : float
        Ratio of vaccinated employees.
    r_vaccination_ratio : float
        Ratio of vaccinated residents.
        
    Returns:
    --------
    ensemble_results : pandas DataFrame
        Data Frame holding the observable of interest of the ensemble, namely
        the number of infected students and teachers.
    '''
    
    ensemble_results = pd.DataFrame()
    ensemble_runs = pd.DataFrame()
    for run in range(1, N_runs + 1):
        model = run_model(measures, simulation_params, contact_network_src,
                   test_type, index_case, e_screen_interval, r_screen_interval,
                   e_vaccination_ratio, r_vaccination_ratio)
        
        # collect the statistics of the single run
        R0, _ = af.calculate_finite_size_R0(model)
        infected_employees = af.count_infected(model, 'employee')
        infected_residents = af.count_infected(model, 'resident')
        data = model.datacollector.get_model_vars_dataframe()
        N_resident_screens_reactive = data['screen_residents_reactive'].sum()
        N_resident_screens_preventive = data['screen_residents_preventive'].sum()
        N_employee_screens_reactive = data['screen_employees_reactive'].sum()
        N_employee_screens_preventive = data['screen_employees_preventive'].sum()
        N_diagnostic_tests = data['N_diagnostic_tests'].max()
        N_preventive_screening_tests = data['N_preventive_screening_tests'].max()
        N_tests = N_diagnostic_tests + N_preventive_screening_tests
        transmissions = sum([a.transmissions for a in model.schedule.agents])
        pending_test_infections = data['pending_test_infections'].max()
        undetected_infections = data['undetected_infections'].max()
        predetected_infections = data['predetected_infections'].max()
        duration = len(data)
        N_agents = len(model.schedule.agents)
        test_rate = N_tests / duration / N_agents
        
        data['run'] = run
        data['step'] = range(0, len(data))
        ensemble_runs = pd.concat([ensemble_runs, data])

        # add run results to the ensemble results
        ensemble_results = ensemble_results.append({ 
                  'R0':R0,
                  'infected_residents':infected_residents,
                  'infected_employees':infected_employees,
                  'N_resident_screens_reactive':N_resident_screens_reactive,
                  'N_employee_screens_reactive':N_employee_screens_reactive,
                  'N_resident_screens_preventive':N_resident_screens_preventive,
                  'N_employee_screens_preventive':N_employee_screens_preventive,
                  'N_diagnostic_tests':N_diagnostic_tests,
                  'N_preventive_tests':N_preventive_screening_tests,
                  'N_tests':N_tests,
                  'transmissions':transmissions,
                  'pending_test_infections':pending_test_infections,
                  'undetected_infections':undetected_infections,
                  'predetected_infections':predetected_infections,
                  'duration':duration,
                  'N_agents':N_agents,
                  'test_rate':test_rate},
                ignore_index=True)
        
    ensemble_runs = ensemble_runs.reset_index(drop=True)
    ensemble_runs.to_csv(join(ensmbl_dst,
            'test-{}_index-{}_esi-{}_rsi-{}_evr-{}_rvr-{}.csv'\
        .format(test_type, index_case, e_screen_interval, r_screen_interval,
                e_vaccination_ratio, r_vaccination_ratio)), index=False)
        
    return ensemble_results


def evaluate_ensemble(ensemble_results, index_case, test_type=None, 
                      e_screen_interval=None, r_screen_interval=None,
                      e_vaccination_ratio=0, r_vaccination_ratio=0):
    '''
    Utility function to calculate ensemble statistics.
    
    Parameters:
    -----------
    ensemble_results: pandas DataFrame
        Data Frame holding the observable of interest of the ensemble, namely
        the number of transmissions from the index case R0, the number of
        infected residents and employees, the number of reactive and preventive
        screens of the resident and employee agent groups, the number of 
        diagnostic and preventive tests used, the overall number of
        transmissions, the number of transmissions that occured while a test
        was pending, the number of infections that occured directly after an
        agent had been tested negative (i.e. because of lacking test
        sensitivity), the number of infections that were detected before an 
        agend became infectious, and the duration of the simulation.
    test_type : string
        Test technology used in the preventive screening. Available test
        technologies are listed in the module testing_strategy.py
    index_case : string
        Agent group from which the index case is drawn. Can be "employee" or
        "resident".
    e_screen_interval : integer
        Interval (in days) of the preventive testing screens in the employee
        agent group. Can be [2, 3, 7, None].
    r_screen_interval : integer
        Interval (in days) of the preventive testing screens in the resident
        agent group. Can be [2, 3, 7, None].
    e_vaccination_ratio : float
        Ratio of vaccinated employees.
    r_vaccination_ratio : float
        Ratio of vaccinated residents.
        
        
    Returns:
    --------
    row : dictionary
        Dictionary holding the values for the screening parameters (test_type, 
        index_case, e_screen_interval, r_screen_interval) and the values of the 
        respective observables of interest.
    '''
    # add ensemble statistics to the overall results
    row = {'test_type':test_type,
           'index_case':index_case,
           'resident_screen_interval':r_screen_interval,
           'employee_screen_interval':e_screen_interval,
           'employee_vaccination_ratio':e_vaccination_ratio,
           'resident_vaccination_ratio':r_vaccination_ratio}
    
    for col in ['R0', 'infected_residents', 'infected_employees', 
                'N_resident_screens_reactive', 'N_employee_screens_reactive',
                'N_resident_screens_preventive', 'N_employee_screens_preventive',
                'N_diagnostic_tests', 'N_preventive_tests', 'transmissions', 
                'pending_test_infections', 'undetected_infections',
                'predetected_infections', 'duration', 'N_tests', 'N_agents', 'test_rate']:

        row.update(af.get_statistics(ensemble_results, col))
    
    return row