# Agent based simulations of the spread of COVID 19 in nursing homes
This repository includes code to run agent based simulations for the spread of a disease (in this case SARS-CoV-2) in a nursing home. The simulations make use of the [small community SEIRX (scseirx) package](https://pypi.org/project/scseirx/) that was developed for this application and the application in the context of [schools](https://www.medrxiv.org/content/10.1101/2021.04.13.21255320v1).  
**Note**: before running any code in this repository, you will have to download the corresponding "data" and "plots" foder from the [OSF-repository](https://osf.io/hyd4r/) belonging to this project and move them to this repository.

## Reference
To reference this work, please cite this [preprint](https://osf.io/5xdqe/). 

## Data
Whenever folders are referenced in the following description, they are located in the [OSF-repository](https://osf.io/hyd4r/) belonging to this project. Exceptions to that rule are indicated accordingly.

## Repository content
The content of the repository is divided into four parts:
1. The construction of contact networks between residents and employees, based on real-world data from living conditions in nursing homes.
2. The calibration of the simulation using empirical observations of SARS-CoV-2 outbreaks in nursing homes.
3. Usage of the calibrated simulation to investigate different testing strategies to keep outbreaks in check.
4. Usage of the calibrated simulation to investigate the impact of different levels of vaccination in the resident and employee population on outbreak sizes.
Simulation parameters used for parts (2), (3) and (4) are stored in the folder ```/params``` in this repository.

### Contact networks
The data used to create the contact networks can be found at ```data/living_conditions```. The tables contain the living area, room number and table number for each resident on a weekly basis over a time period of five months between February 2020 and July 2020. To create the contact networks, the resident living conditions in the first week of the data set are used. Contact networks are created in the script ```parse_static_networks.ipynb``` and stored in the folder ```data/contact_networks```. 

### Calibration
The data used to calibrate the simulation can be found at ```data/outbreaks```. The tables contain information about the time in days and the cumulative number of infected employees and residents for a total of four empirically observed outbreaks. In addition, the total number of residents and employees present in a nursing home at the time of the outbreak for all four outbreaks is stored in the file ```data/outbreaks/total_agent_numbers.csv```.  

To calibrate the simulation, the transmission risk for contacts of type "loose" (far) and "intermediate" is varied until the simulated cumulative numbers of infected employees and residents over time match the empirically observed numbers as close as possible. The calibration is performed in the script ```data_creation_calibration.ipynb```. Here, the parameter space is first scanned in coarse steps and scanning is then refined around the optimum detected in the coarse scan. Simulation parameters and measures are stored in the files ```calibration_simulation_parameters.json``` and ```calibration_measures.json```. Aggregated ensemble results are stored in the folder ```data/calibration_results``` for the coarse (```calibration_results_coarse_N1000.csv```) and fine (```calibration_results_fine_N5000.csv```) parameter grid scan. Individual ensembles are stored in the compressed archive ```ensembles.zip``` and need to be unzipped before they can be analysed. Calibration results are visualised and compared to the empirical data in the script ```analysis_calibration.ipynb```. Visualizations are saved to the folder ```plots```.

### Investigation of testing strategies
The main functionality for running simulations resides in the ```data_creation_functions.py``` module. The main functionality for data analysis resides in the ```analysis_functions.py``` module. All results are provided for the variant B.1.1.7 (UK variant). Parameter files and result files corresponding to the variant are indicated with a ```_UK_variant``` in the respective file names.

#### Baseline outbreak sizes with no interventions
To investigate, what the baseline outbreak sizes with no intervention measures are in our system, we use the calibrated system to run a set of simulations with resident and employee index cases, respectively, but no prevention measures at all. Simulations and analysis are performed in the script ```data_creation_no_intervention.ipynb``` and ```analysis_no_intervention.ipynb```. Simulation parameters and measures are stored in the files ```no_intervention_simulation_parameters.json``` and ```no_intervention_measures.json```. Aggregated simulation results are stored in the folder ```data/simulation_results```, file ```simulation_results_no_measures_5000.csv```.  Individual ensembles are stored in the compressed archive ```ensembles/no_measures.zip``` and need to be unzipped before they can be analysed.

#### Outbreak sizes with only TTI
To investigate, what the outbreak sizes with only TTI (test-trace-isolate) are in our system, we use the calibrated system to run a set of simulations with resident and employee index cases, respectively and only diagnostic testing and isolation of contact persons. Simulations and analysis are performed in the script ```data_creation_TTI.ipynb``` and ```analysis_TTI.ipynb```. Simulation parameters and measures are stored in the files ```TTI_simulation_parameters.json``` and ```TTI_measures.json```. Aggregated simulation results are stored in the folder ```data/simulation_results```, file ```simulations_TTI_5000.csv```. Individual ensembles are stored in the compressed archive ```ensembles/TTI.zip``` and need to be unzipped before they can be analysed.

#### Test technologies and turnover times
We then use calibrated system to investigate different testing strategies. Of interest are both (i) the testing technology used (antigen tests, PCR tests or LAMP tests) and (ii) the frequency of preventive resident and employee testing. Simulations are performed in the script ```data_creation_testing_strategy.ipynb```. Simulation parameters and measures are stored in the files ```testing_strategy_simulation_parameters.json``` and ```testing_strategy_measures.json```. Aggregated simulation results are stored in the folder ```data/simulation_results```, file ```simulations_testing_strategy.csv```. Individual ensembles are stored in the compressed archive ```ensembles/testing_strategy.zip``` and need to be unzipped before they can be analysed.  

Simulation results are analysed and visualised in the script ```analysis_testing_strategy.ipynb```. This script also creates the tables A2, A3, A8 and A9 of the supplementary materials. Tables are stored in the folder ```data/results_tables```, files ```table_testing_strategy_employee.csv``` and ```table_testing_strategy_employee.csv```. Visualizations are saved to the folder ```plots```. 

### Investigation of vaccinations
The main functionality for running simulations resides in the ```data_creation_functions.py``` module. The main functionality for data analysis resides in the ```analysis_functions.py``` module. All results are provided for the variant B.1.1.7 (UK variant). Parameter files and result files corresponding to the variant are indicated with a ```_UK_variant``` in the respective file names.

#### Different vaccination prevalence
To investigate the effect of vaccinations only, we use the calibrated system to perform simulations for different ratios of vaccinated employees and residents. Simulation parameters and measures are stored in the files ```vaccination_simulation_parameters.json``` and ```vaccination_measures.json```. Visualizations are saved to the folder ```plots```.  Simulations are performed in the script ```data_creation_vaccination.ipynb```. Aggregated simulation results are stored in the folder ```data/simulation_results```, file ```simulations_vaccination_rate_5000.csv```. Individual ensembles are stored in the compressed archive ```ensembles/vaccination.zip``` and need to be unzipped before they can be analysed.  

Simulation results are analysed and visualised in the script ```analysis_vaccination.ipynb```. Visualizations are saved to the folder ```plots```. 

#### Different vaccination scenarios and testing strategies
To investigate the effect of vaccinations and testing strategies together, we first define four different vaccination scenarios:
* In a scenario where vaccinations are hard to come by and residents are prioritised to get the vaccine, 50% of residents and 0% of employees are vaccinated.
* In a scenario where vaccinations are hard to come by and employees are prioritised to get the vaccine, 0% of residents and 50% of employees are vaccinated.
* In a scenario where vaccinations are ubiquitous but employees are hesitant to get vaccinated, 90% of residents and 50% of employees are vaccinated.
* In a scenario where vaccinations are ubuqitous and both a majority of residents and employees want to get vaccinated, 90% of residents and 90% of employees are vaccinated.

For each of these four vaccination scenarios, we use the calibrated system to perform simulations for different testing strategies. We focus on PCR and antigen tests with same-day turnover for preventive screening, as these are the tests that are most widely used in practice. Simulations are performed in the script ```data_creation_testing_and_vaccination.ipynb```. Simulation parameters and measures are stored in the files ```testing_and_vaccination_simulation_parameters.json``` and ```testing_and_vaccination_measures.json```. Aggregated simulation results are stored in the folder ```data/simulation_results```, file ```simulations_testing_and_vaccination_5000.csv```. Individual ensembles are stored in the compressed archive ```ensembles/testing_and_vaccination.zip``` and need to be unzipped before they can be analysed.    

Simulation results are analysed and visualised in the script ```analysis_testing_and_vaccination.ipynb```.  Visualizations are saved to the folder ```plots```.  