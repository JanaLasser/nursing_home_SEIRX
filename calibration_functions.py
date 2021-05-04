from scseirx.model_nursing_home import SEIRX_nursing_home
import scseirx.analysis_functions as af
import pandas as pd
import numpy as np
import networkx as nx
from os.path import join


def compose_agents(measures, simulation_params):
    '''
    Utility function to compose agent dictionaries as expected by the simulation
    model as input from the dictionary of prevention measures.
    
    Parameters
    ----------
    measures : dictionary
        Dictionary of prevention measures. Needs to include the fields 
        employee_mask and resident_mask. 
    simulation_params: dictionary
        Dictionary of (epidemiological) simulation parameters. Needs to include
        the fields employee_index_probability and resident_index_probability.
        
    Returns
    -------
    agent_types : dictionary of dictionaries
        Dictionary containing the fields "screening_interval", 
        "index_probability" and "mask" for the agent groups "employee" and
        "resident".
    
    '''
    agent_types = {
        'employee':{
            'screening_interval':measures['employee_screening_interval'],
            'index_probability':simulation_params['employee_index_probability'],
            'mask':measures['employee_mask']},

        'resident':{
            'screening_interval':measures['resident_screening_interval'],
            'index_probability':simulation_params['resident_index_probability'],
            'mask':measures['resident_mask']},
    }
    return agent_types


def run_model(contact_weight, measures, simulation_params, contact_network_src, 
              N_steps=100):
    '''
    Runs a simulation with an SEIRX_nursing_home model 
    (see https://pypi.org/project/scseirx/1.3.0/), given a set of parameters.
    
    Parameters:
    -----------
    contact_weight : float
        Multiplicative weight for the decrease of infection risk during a 
        contact of type intermediate and far, as compared to a household (close)
        contact. Note: differentiation between "intermediate" and "far" does not
        lead to improved calibration accuracy and is therefore not done.
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
    N_steps : integer
        Number of maximum steps per run. Since the empirically observed
        outbreaks only last for a maximum of 27 days, it is more than sufficient
        for the purpose of comparison to truncate simulations at 100 steps.
        
    Returns
    -------
    model : SEIRX_nursing_home model instance holding a completed simulation run
        and all associated data.
    '''
    # in the calibration szenario, all index cases in nursing homes are 
    # assumed to be introduced through employees, since residents did not
    # have any visitors and were often not allowed to go outside
    index_case = 'employee'
        
    agent_types = compose_agents(measures, simulation_params)

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
                 {'very_far':0, 'far':contact_weight,
                  'intermediate':contact_weight, 'close':1},
      K1_contact_types = measures['K1_contact_types'],
      diagnostic_test_type = measures['diagnostic_test_type'],
      preventive_screening_test_type = \
                 measures['preventive_screening_test_type'],
      follow_up_testing_interval = \
                 measures['follow_up_testing_interval'],
      liberating_testing = measures['liberating_testing'],
      index_case = index_case,
      agent_types = agent_types, 
      age_transmission_risk_discount = \
                simulation_params['age_transmission_discount'],
      age_symptom_discount = simulation_params['age_symptom_discount'],
      mask_filter_efficiency = measures['mask_filter_efficiency'],
      transmission_risk_ventilation_modifier = \
                measures['transmission_risk_ventilation_modifier'],)

    # run the model until the outbreak is over
    for i in range(N_steps):
        model.step()
        
    return model


def run_ensemble(N_runs, contact_weight, measures, simulation_params, 
                 contact_network_src, emp_data_src, comp_period, ensmbl_dst):
    '''
    Utility function to run an ensemble of simulations during parameter 
    calibration. The parameter "contact_weight" is calibrated here.
    
    Parameters:
    ----------
    N_runs : integer
        Number of individual simulation runs in the ensemble.
    contact_weight : float
        Multiplicative weight for the decrease of infection risk during a 
        contact of type intermediate and far, as compared to a household (close)
        contact. Note: differentiation between "intermediate" and "far" does not
        lead to improved calibration accuracy and is therefore not done.
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
    emp_data_src : string
        Absolute or relative path pointing to the location of the empirically
        observed outbreaks.
    comp_period : integer
        Number of days in the empirical observations (and time steps in the
        simulation) that will be used to compare the two. 
        
    Returns:
    --------
    ensemble_results : pandas DataFrame
        Data Frame holding the observable of interest of the ensemble, namely
        the number of infected students and teachers.
    '''
    # get the test turnover time in days (rather than as string)
    turnover_times = {'same_day':0, 'one_day':1, 'two_day':2, 'three_day':3}
    turnover = measures['diagnostic_test_type'].split('_')[0] + '_' + \
               measures['diagnostic_test_type'].split('_')[1]
    turnover = turnover_times[turnover]
    
    # data of observed outbreaks, holding information about the cumulative
    # number of infected residents and employees over time
    observed_outbreaks = pd.DataFrame()
    for ob in range(1, 5):
        tmp = pd.read_csv(join(emp_data_src, 'outbreak_{}.csv'.format(ob)))
        tmp['outbreak'] = ob
        if len(tmp) < comp_period:
            for i in range(0, comp_period - len(tmp) + 1):
                tmp = tmp.append(tmp.loc[tmp.index[-1]], ignore_index=True)
                tmp['t'] = range(len(tmp))
        observed_outbreaks = pd.concat([observed_outbreaks, tmp])
    observed_outbreaks = observed_outbreaks.reset_index()
        
    # the total number of employees and residents in the empirical setting
    # is needed to normalise the obseved number of infected agents correctly
    emp_agent_numbers = pd.read_csv(join(emp_data_src,\
                            'total_agent_numbers.csv'), index_col='outbreak')
    
    ensemble_results = pd.DataFrame()
    ensemble_runs = pd.DataFrame()
    for run in range(1, N_runs + 1):
        model = run_model(contact_weight, measures, simulation_params, 
                          contact_network_src)
        
        # collect the statistics of the single run
        R0, _ = af.calculate_finite_size_R0(model)
        # subtract the index case from the employee infected count
        infected_employees = af.count_infected(model, 'employee') - 1
        infected_residents = af.count_infected(model, 'resident')
        infected_total = infected_employees + infected_residents
        
        data = model.datacollector.get_model_vars_dataframe()
        
        # If there are no transmissions from the source case, we discard the run
        # as the "outbreak" will go unobserved. 
        # There can be rare cases in which there are transmissions but all cases
        # are asymptomatic. In these cases, there are no diagnostic tests and
        # the outbreak will remain "unobserved". We discard these cases
        if infected_total == 0 or data['N_diagnostic_tests'].sum() == 0:
            continue
            
        # get the simulation data from the point in time at whicht he first test
        # result arrived, i.e. the day the first diagnostic test was performed
        # plus the turnover time of the diagnostic test technology
        data = data.loc[data[data['N_diagnostic_tests'] > 0].index[0]+turnover:]
        data = data.reset_index()
        
        # There can be rare cases, where the source case is symptomatic, with
        # a very short infection duration and does not infect another agent.
        # In these cases, the "outbreak" duration will be shorter than the
        # required period until the first diagnostic test + the test turnover
        # time and therefore the data-frame will have zero-length after the
        # above truncating operation. We discard these cases.
        if len(data) == 0:
            continue
            
        data['run'] = run
        data['step'] = range(0, len(data))
        ensemble_runs = pd.concat([ensemble_runs, data])
        
        N_resident_screens_reactive = data['screen_residents_reactive'].sum()
        N_resident_screens_preventive = data['screen_residents_preventive'].sum()
        N_employee_screens_reactive = data['screen_employees_reactive'].sum()
        N_employee_screens_preventive = data['screen_employees_preventive'].sum()
        N_diagnostic_tests = data['N_diagnostic_tests'].max()
        N_preventive_screening_tests = data['N_preventive_screening_tests'].max()
        transmissions = sum([a.transmissions for a in model.schedule.agents])
        pending_test_infections = data['pending_test_infections'].max()
        undetected_infections = data['undetected_infections'].max()
        predetected_infections = data['predetected_infections'].max()
        
        row = {
            'R0':R0,
            'infected_residents':infected_residents,
            'infected_employees':infected_employees,
            'infected_total':infected_total,
            'N_resident_screens_reactive':N_resident_screens_reactive,
            'N_employee_screens_reactive':N_employee_screens_reactive,
            'N_resident_screens_preventive':N_resident_screens_preventive,
            'N_employee_screens_preventive':N_employee_screens_preventive,
            'N_diagnostic_tests':N_diagnostic_tests,
            'N_preventive_tests':N_preventive_screening_tests,
            'transmissions':transmissions,
            'pending_test_infections':pending_test_infections,
            'undetected_infections':undetected_infections,
            'predetected_infections':predetected_infections}
        
        # calculate the difference between the characteristics of the run and
        # the four observed outbreaks
        for outbreak in range(1, 5):
            employee_distance, resident_distance = \
                calculate_infected_over_time_difference(data, outbreak, 
                        emp_agent_numbers, observed_outbreaks, comp_period)
        
            row.update({
                    'infected_employee_distance_outbreak_{}'.format(outbreak):\
                        employee_distance,
                    'infected_resident_distance_outbreak_{}'.format(outbreak):\
                        resident_distance
                    })
        
        # add run results to the ensemble results
        ensemble_results = ensemble_results.append(row, ignore_index=True)
    
    ensemble_runs = ensemble_runs.reset_index(drop=True)
    ensemble_runs.to_csv(join(ensmbl_dst, 'cw-{:1.3f}.csv'\
            .format(contact_weight)), index=False)
    return ensemble_results


def calculate_infected_over_time_difference(sim_data, outbreak, 
                        emp_agent_numbers, observed_outbreaks, comp_period):
    '''
    Calculates the difference between the expected and simulated number of 
    infected employees and residents (separately), given the empirically 
    observed outbreak data in Austrian nursing homes and the number of infected 
    in a given simulation.
    
    Parameters:
    -----------
    sim_data : pandas DataFrame
        Data frame with the output of a single simulation in terms of the number
        of susceptible, exposed, infected and recovered employees and residents
        at every time-step of the simulation.
    outbreak : integer
        Label of the empirically observed outbreak the current simulation run
        will be compared to.
    emp_agent_numbers : pandas DataFrame
        Data frame holding information about the total number of employees and
        residents present in the homes in which the empirically observed
        outbreaks took place.
    observed_outbreaks : pandas DataFrame
        Data frame holding information about the empirically observed outbreaks
        in terms of the number of employees and residents that tested positive
        over time.
    comp_period : integer
        Number of days in the empirical observations (and time steps in the
        simulation) that will be used to compare the two. 
        
    Returns:
    --------
    distances : tuple
        Tuple of floats holding the sum of squares of the difference between the
        cumulative infected employees [residents] in the simulation and in the
        empirically observed outbreak.
    '''
    # ********* data from the simulation ********* 
    N_sim_employees = \
            sim_data.loc[0]['S_employee'] + sim_data.loc[0]['E_employee'] +\
            sim_data.loc[0]['I_employee'] + sim_data.loc[0]['R_employee']
    N_sim_residents = \
            sim_data.loc[0]['S_resident'] + sim_data.loc[0]['E_resident'] +\
            sim_data.loc[0]['I_resident'] + sim_data.loc[0]['R_resident']
    
    # calculate the total number of all employees [residents] that are exposed,
    # infected or recovered over time (cumulative infected)
    sim_data['cumulative_I_employee'] = sim_data['E_employee'] + \
            sim_data['I_employee'] + sim_data['R_employee']
    sim_data['cumulative_I_resident'] = sim_data['E_resident'] + \
            sim_data['I_resident'] + sim_data['R_resident']

    # normalise the number of infected with the total number of employees 
    # [residents] in the simulation
    sim_data['cumulative_I_employee'] = sim_data['cumulative_I_employee'] / \
                                        N_sim_employees
    sim_data['cumulative_I_resident'] = sim_data['cumulative_I_resident'] / \
                                        N_sim_residents

    # ********* empirical outbreak data ********* 
    N_emp_employees = emp_agent_numbers.loc[outbreak,'total_employees']
    N_emp_residents = emp_agent_numbers.loc[outbreak,'total_residents']
    emp_data = observed_outbreaks[observed_outbreaks['outbreak'] == outbreak]\
        .copy()
    
    # normalisation
    emp_data['cumulative_I_employee'] = emp_data['cumulative_I_employee'] / \
                                        N_emp_employees
    emp_data['cumulative_I_resident'] = emp_data['cumulative_I_resident'] / \
                                        N_emp_residents

    # truncate both vectors to have the same length, do not use more than the
    # first comp_period days to compare outbreaks.
    min_T = min(len(emp_data['cumulative_I_employee']),\
                len(sim_data['cumulative_I_employee']))
    if min_T > comp_period: min_T = comp_period
    sim_data = sim_data[0:min_T]
    emp_data = emp_data[0:min_T]
    
    # the error metric is the sum of the squared distances between the 
    # cumulative number of infected employees [residents] in the simulation and
    # the observed outbreak
    diff_employee = (emp_data['cumulative_I_employee'].values - \
                     sim_data['cumulative_I_employee'].values) ** 2
    distance_employee = diff_employee.sum()
    diff_resident = (emp_data['cumulative_I_resident'].values - \
                     sim_data['cumulative_I_resident'].values) ** 2
    distance_resident = diff_resident.sum()
    
    distances = (distance_employee, distance_resident)
    return distances


def evaluate_ensemble(ensemble_results, contact_weight):
    '''
    Utility function to calculate ensemble statistics of a range of observables.
    
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
    contact_weight : float
        Multiplicative weight for the decrease of infection risk during a 
        contact of type intermediate and far, as compared to a household (close)
        contact. Note: differentiation between "intermediate" and "far" does not
        lead to improved calibration accuracy and is therefore not done.
  
    Returns:
    --------
    row : dictionary
        Dictionary holding the values for the calibration parameter 
        (contact_weight) and the values of the  respective observables of 
        interest.
    '''
    
    # add simulation ensemble statistics to the overall results
    row = {'contact_weight':contact_weight}
    
    eval_cols = ['R0', 'infected_residents','infected_employees','infected_total',
                'N_resident_screens_reactive', 'N_employee_screens_reactive',
                'N_resident_screens_preventive', 'N_employee_screens_preventive',
                'N_diagnostic_tests', 'N_preventive_tests', 'transmissions', 
                'pending_test_infections', 'undetected_infections',
                'predetected_infections']
    
    eval_cols = eval_cols + \
        ['infected_employee_distance_outbreak_{}'.format(i) for i in range(1, 5)]
    eval_cols = eval_cols + \
        ['infected_resident_distance_outbreak_{}'.format(i) for i in range(1, 5)]
    
    for col in eval_cols:
        row.update(af.get_statistics(ensemble_results, col))

    row.update({'runs':len(ensemble_results)})
        
    return row
