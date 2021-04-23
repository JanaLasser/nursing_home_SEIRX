# Agent based simulations of the spread of COVID 19 in nursing homes
This repository includes code to run agent based simulations for the spread of a disease (in this case SARS-CoV-2) in a nursing home. The simulations make use of the [small community SEIRX (scseirx) package](https://pypi.org/project/scseirx/) that was developed for this application and the application in the context of schools [link needed].

## Reference
To reference this work, please cite this [preprint](https://osf.io/5xdqe/). 

## Data
Whenever folders are referenced in the following description, they are located in the OSF-repository belonging to this project [link needed]. Exceptions to that rule are indicated accordingly.

## Repository content
The content of the repository is divided into four parts:
1. The construction of contact networks between residents and employees, based on real-world data from living conditions in nursing homes.
2. The calibration of the simulation using empirical observations of SARS-CoV-2 outbreaks in nursing homes.
3. Usage of the calibrated simulation to investigate different testing strategies to keep outbreaks in check.
4. Usage of the calibrated simulation to investigate the impact of different levels of vaccination in the resident and employee population on outbreak sizes.
Simulation parameters used for parts (2), (3) and (4) are stored in the folder ```/params``` in this repository.

### Contact networks
The data used to create the contact networks can be found at ```data/living_conditions``` [not yet uploaded, pending data protection clearance]. The tables contain the living area, room number and table number for each resident on a weekly basis over a time period of five months between February 2020 and July 2020. To create the contact networks, the resident living conditions in the first week of the data set are used. Contact networks are created in the script ```parse_static_networks.ipynb``` and stored in the folder ```data/contact_networks```. 

### Calibration
The data used to calibrate the simulation can be found at ```data/outbreaks```. The tables contain information about the time in days and the cumulative number of infected employees and residents for a total of four empirically observed outbreaks. In addition, the total number of residents and employees present in a nursing home at the time of the outbreak for all four outbreaks is stored in the file ```data/outbreaks/total_agent_numbers.csv```.  

To calibrate the simulation, the transmission risk for contacts of type "loose" (far) and "intermediate" is varied until the simulated cumulative numbers of infected employees and residents over time match the empirically observed numbers as close as possible. The calibration is performed in the script ```data_creation_calibration.ipynb```. Here, the parameter space is first scanned in coarse steps and scanning is then refined around the optimum detected in the coarse scan. Calibration results are stored in the folder ```data/calibration_results``` Calibration results are visualised and compared to the empirical data in the script ```analysis_calibration.ipynb```. Visualizations are saved to the folder ```plots```.

### Investigation of testing strategies
The calibrated simulation is used to investigate different testing strategies. Of interest are both (i) the testing technology used (antigen tests, PCR tests or LAMP tests) and (ii) the frequency of preventive resident and employee testing. Simulations are performed in the script ```data_creation_testing_strategy.ipynb```. Simulation results are stored in the folder ```data/simulation_results```. Simulation results are analysed and visualised in the script ```analysis_testing_strategy.ipynb```. Visualizations are saved to the folder ```plots```. 

### Investigation of vaccinations
TODO