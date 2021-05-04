def plot_errors(ax, results, best_weight, xmin=0.1):
    '''
    Helper function to plot the individual (employee, resident) and total error 
    terms given the results from a parameter screen in the calibration process.
    
    Parameters:
    -----------
    ax : matplotlib axis object
        Axis on which the errors will be plotted.
    results : pandas DataFrame
        Data Frame holding the results of a calibration run (parameter scan).
        Needs to include the columns "infected_employee_distance_total_mean" and
        "infected_resident_distance_total_mean", as well as an index 
        "contact_weight".
    best_weight : float
        Optimum for the contact_weight that minimizes the total error.
    xmin : float
        Optional parameter to adjust the extent of the x-axis in the plot.
    '''
    ax.plot(results.index.get_level_values('contact_weight'), 
            results['infected_employee_distance_total_mean'],
           'o-', color='FireBrick', label='$e_1$')

    ax.plot(results.index.get_level_values('contact_weight'), 
            results['infected_resident_distance_total_mean'],
           'o-', color='DarkBlue', label='$e_2$')

    ax.plot(results.index.get_level_values('contact_weight'), 
            results['distance_total'],
           'o-', color='grey', label='$E$')
    ax.fill_between(results.index.get_level_values('contact_weight'), 
            results['distance_total'],
            results['distance_total'] - \
            results['distance_total_std'], 
            color='grey', alpha=0.2)
    ax.fill_between(results.index.get_level_values('contact_weight'), 
            results['distance_total'],
            results['distance_total'] + \
            results['distance_total_std'], 
            color='grey', alpha=0.2)

    ymax = 20
    ax.plot([best_weight, best_weight], [0, ymax],'--', color='k', label='optimum')
    ax.text(0.165, 12.5, '$q_1$ = {:1.2f}'.format(1 - best_weight))

    ax.legend(loc=9, fontsize=11.3)
    ax.set_xlim(xmin, 0.5)
    ax.set_ylim(-0.5, ymax)
    ax.set_xticks([0.1, 0.2, 0.3, 0.4, 0.5])
    ax.set_xticklabels([1 - 0.1, 1 - 0.2, 1 - 0.3, 1 - 0.4, 1 - 0.5])
    ax.set_xlabel('$q_1$', fontsize=16)
    ax.set_ylabel('sum of squared differences', fontsize=16)
    

def plot_emp_sim_data(ax, emp_data, agg, comp_period):
    '''
    Helper function to plot the empirical and simulated data for comparison.
    
    Parameters:
    -----------
    ax : matplotlib axis object
        Axis on which the errors will be plotted.
    emp_data : pandas DataFrame
        Data Frame holding the data of the empirically observed oubtreaks. Needs
        to include the columns "cumulative_I_employee" and 
        "cumulative_I_resident", holding information about the sumulative sum of 
        exposed, infected and recovered employees [residents] in a nursing home 
        outbreak. 
    agg : pandas DataFrame
        Data Frame holding the results of the simulated ensemble with the 
        optimal contact weight, aggregated by day. Needs to hold the columns
        "I_total_employee" and "I_total_residents", which are the sum of
        exposed, infected and recovered employees [residents] at a given step
        (day) in the simulation.
    comp_period : int
        Number of days that were used to compare the empirical and simulated
        data. Used here to determine the extent of the x-axis in the plot.
    '''
    colors = ['FireBrick', 'DarkBlue']
    agent_types = ['employee', 'resident']
    

    # simulated data
    for agent_type, color in zip(agent_types, colors):

        ax.plot(agg.index, agg['I_total_{}'.format(agent_type)]['mean'],
                color=color, label='sim. I {}.'.format(agent_type[0:3]))
        ax.fill_between(agg.index, agg['I_total_{}'.format(agent_type)]['percentile_10'],
                                   agg['I_total_{}'.format(agent_type)]['percentile_90'],
                                   color=color, alpha=0.2)

    ax.plot(np.zeros(1), np.zeros([1, 2]), color='w', alpha=0, label=' ')
    
    # empirical data
    for ob, marker in zip(list(range(1, 5)), ['+', '*', '^', 'o']):
        ob_data = emp_data[emp_data['outbreak'] == ob]
        for agent_type, color in zip(agent_types, colors):
            ax.scatter(ob_data['t'], ob_data['cumulative_I_{}'.format(agent_type)],
                       marker=marker, color=color, alpha=0.3,
                       label='OB {}: I {}.'.format(ob, agent_type[0:3]))

    ax.set_xlim(0, comp_period + 0.5)
    ax.set_ylim(0, 27)
    
    # legend
    from matplotlib.patches import Patch
    from matplotlib.pyplot import Line2D
    alpha=0.2
    emp_handle = Patch(facecolor='FireBrick',label='employees', alpha=alpha)
    res_handle = Patch(facecolor='DarkBlue',label='residents', alpha=alpha)

    sim_handle = Line2D((0,1),(0,0), color='k', linewidth=5)
    case1_handle = Line2D((0,1),(0,0), color='k', marker='o', linestyle='',
                              markersize=10)
    case2_handle = Line2D((0,1),(0,0), color='k', marker='*', linestyle='',
                              markersize=10)
    case3_handle = Line2D((0,1),(0,0), color='k', marker='X', linestyle='',
                              markersize=10)
    case4_handle = Line2D((0,1),(0,0), color='k', marker='^', linestyle='',
                              markersize=10)

    legend = ax.legend([case1_handle, case2_handle,
               case3_handle, case4_handle, res_handle, emp_handle, sim_handle],
              ['outbreak 1', 'outbreak 2', 'outbreak 3', 'outbreak 4', 'residents',
               'employees', 'simulation'], loc=9, ncol=2, fontsize=11.3)
    
    ax.set_xlabel('days', fontsize=16)
    ax.set_ylabel('N infected', fontsize=16)