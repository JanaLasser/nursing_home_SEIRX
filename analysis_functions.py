import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from os.path import join
import matplotlib.gridspec as gridspec
import seaborn as sns
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable

# Since these plots were also used to advise Austrian decision makers, we 
# provide plot label mappings to German for easier communication of results.
label_map = {
    'english':{
        'metric_name_map':{
            'infected_residents_mean':'follow-up cases residents (mean)',
            'infected_residents_median':'follow-up cases residends (median)',
            'infected_residents_0.90':'follow-up cases residents (90th percentile)',
            'infected_residents':'follow-up cases residents',
            'R0_mean':'$R_\\mathrm{eff}$',
            'R0_mean':'$R_\\mathrm{eff}$',
            'test_rate_mean':'tests / day / agent (mean)'},

        'index_case_map':{
            'employee':'index case employee',
            'resident':'index case resident'},
        
        'frequency_name_map':{
            np.nan:'never',
            2:'3 times\na week',
            3:'twice\na week',
            7:'once\na week'},
    
        'test_name_map':{
            'same_day_antigen':'same-day antigen',
            'one_day_PCR':'one day PCR',
            'two_day_PCR':'two days PCR',
            'same_day_PCR':'same-day PCR',
            'same_day_LAMP':'same-day RT-LAMP',
            None:''},
        
        'index_case_label':'index case',
        'employee_label':'employee',
        'resident_label':'resident',
        
        'xlabels':{
            'vaccination':'vaccinated employees',
            'testing_and_vaccination':'screening frequency employees',
            'testing_strategy':'screening frequency employees'},
        
        'ylabels':{
            'vaccination':'vaccinated residents',
            'testing_and_vaccination':'screening frequency residents',
            'testing_strategy':'screening frequency residents'},
    },
    'german':{
        'metric_name_map':{
            'infected_residents_mean':'Folgefälle BewohnerInnen (Mittelwert)',
            'infected_residents_median':'Folgefälle BewohnerInnen (Median)',
            'infected_residents_0.90':'Folgefälle BewohnerInnen (90. Percentile)',
            'infected_residents':'Folgefälle BewohnerInnen',
            'R0_mean':'$R_\\mathrm{eff}$',
            'R0_mean':'$R_\\mathrm{eff}$',
            'test_rate_mean':'Tests / Tag / Agent (Mittelwert)'},

        'index_case_map':{
            'employee':'kein Besuch (Indexfall MitarbeiterIn)',
            'resident':'Besuch (Indexfall BewohnerIn)'},
        
        'frequency_name_map':{
            np.nan:'nie',
            2:'3 mal\npro Woche',
            3:'2 mal\npro Woche',
            7:'ein mal\npro Woche'},
    
        'test_name_map':{
            'same_day_antigen':'Antigen (selber Tag)',
            'one_day_PCR':'PCR (1 Tag)',
            'two_day_PCR':'PCR (2 Tage)',
            'same_day_PCR':'PCR (selber Tag)',
            'same_day_LAMP':'RT-LAMP (selber Tag)',
            None:''},
        
        'index_case_label':'Indexfall',
        'employee_label':'MitarbeiterIn',
        'resident_label':'BewohnerIn',
        
        'xlabels':{
            'vaccination':'geimpfte MitarbeiterInnen',
            'testing_and_vaccination':'Testfrequenz MitarbeiterInnen',
            'testing_strategy':'Testfrequenz MitarbeiterInnen'},
        
        'ylabels':{
            'vaccination':'geimpfte BewohnerInnen',
            'testing_and_vaccination':'Testfrequenz BewohnerInnen',
            'testing_strategy':'Testfrequenz BewohnerInnen'},
    }
}

def round_decimals_up(number:float, decimals:int=2):
    """
    Returns a value rounded up to a specific number of decimal places.
    """
    import math
    if not isinstance(decimals, int):
        raise TypeError("decimal places must be an integer")
    elif decimals < 0:
        raise ValueError("decimal places has to be 0 or more")
    elif decimals == 0:
        return math.ceil(number)

    factor = 10 ** decimals
    return math.ceil(number * factor) / factor


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
    
    
def get_image(df, subset, screening_params, metric):
    subset = df.loc[subset]
    img = np.zeros((len(screening_params), len(screening_params)))
    for i, p_index in enumerate(screening_params):
        for j, e_index in enumerate(screening_params):
            img[i, j] = subset.loc[e_index, p_index][metric]
    return img

def plot_heatmap(ax, img, screening_params, vmin, vmax, xticks, yticks,
                  xlabel, ylabel, ticklabel_fontsize=9, 
                 cmap=plt.get_cmap('coolwarm')):

    im = ax.imshow(img, origin='lower', vmin=vmin, vmax=vmax, cmap=cmap)
    
    if xticks:
        ax.set_xticks(range(len(screening_params)))
        ax.set_xticklabels(screening_params, fontsize=ticklabel_fontsize)
        ax.set_xlabel(xlabel, fontsize=12)
    else:
        ax.set_xticks([])
    if yticks:    
        ax.set_yticks(range(len(screening_params)))
        ax.set_yticklabels(screening_params, fontsize=ticklabel_fontsize)
        ax.set_ylabel(ylabel, fontsize=12)
    else:
        ax.set_yticks([])
    
    return im

def annotate_heatmap(ax, img):
    for i in range(img.shape[0]):
        for j in range(img.shape[1]):
            y_pos = i - 0.05
            x_pos = j - 0.15
            ax.text(x_pos, y_pos, '{:1.2f}'.format(img[i, j]))
            
            
def load_ensmbl(p, test_name_map, ensmbl_dst, variant='wild_type'):
    test_type, index_case,\
    e_screen_interval, r_screen_interval, testing_strat_label, \
    e_vacc_ratio, r_vacc_ratio, vacc_scen_label = p
    
    ensmbl_dst = join(ensmbl_dst, variant)
    ensmbl_name = 'test-{}_index-{}_esi-{}_rsi-{}_evr-{}_rvr-{}.csv'

    ensmbl_file = ensmbl_name.format(test_type, index_case, e_screen_interval,
                r_screen_interval, e_vacc_ratio, r_vacc_ratio)
    
    ensmbl = pd.read_csv(join(ensmbl_dst, ensmbl_file))
    ensmbl['vaccination_scenario'] = vacc_scen_label
    ensmbl['employee_vaccination_ratio'] = e_vacc_ratio
    ensmbl['resident_vaccination_ratio'] = r_vacc_ratio
    ensmbl['testing_scenario'] = testing_strat_label + '\n' + \
        test_name_map[test_type]
    ensmbl['index_case'] = index_case
    ensmbl['infected_residents'] = ensmbl['E_resident'] + \
                                   ensmbl['I_resident'] + ensmbl['R_resident']
    
    return ensmbl


def load_baseline_ensmbl(p, test_name_map, baseline_ensmbl_dst, variant='wild_type'):
    test_type, index_case, e_screen_interval, r_screen_interval,\
        testing_strat_label, vacc_strat_label = p
    
    e_vacc_ratio = 0
    r_vacc_ratio = 0
    
    baseline_ensmbl_dst = join(baseline_ensmbl_dst, variant)

    try:
        ensmbl_name = 'test-{}_index-{}_esi-{}_rsi-{}.csv'
        ensmbl_file = ensmbl_name.format(test_type, index_case, e_screen_interval,
                    r_screen_interval)
        ensmbl = pd.read_csv(join(baseline_ensmbl_dst, ensmbl_file))
        
    except FileNotFoundError:
        ensmbl_name = 'test-{}_index-{}_esi-{}_rsi-{}_evr-{}_rvr-{}.csv'
        ensmbl_file = ensmbl_name.format(test_type, index_case, e_screen_interval,
                    r_screen_interval, e_vacc_ratio, r_vacc_ratio)
        ensmbl = pd.read_csv(join(baseline_ensmbl_dst, ensmbl_file))
    

    
    ensmbl['vaccination_scenario'] = vacc_strat_label
    ensmbl['testing_scenario'] = testing_strat_label + '\n' + \
        test_name_map[test_type]
    ensmbl['index_case'] = index_case
    ensmbl['infected_residents'] = ensmbl['E_resident'] + \
                                   ensmbl['I_resident'] + ensmbl['R_resident']
    return ensmbl


def load_no_test_ensmbl(p, test_name_map, no_test_ensmbl_dst, variant='wild_type'):
    test_type, index_case, e_screen_interval, r_screen_interval,\
        testing_strat_label, e_vacc_ratio, r_vacc_ratio, vacc_strat_label = p
    
    no_test_ensmbl_dst = join(no_test_ensmbl_dst, variant)

    ensmbl_name = 'test-{}_index-{}_esi-{}_rsi-{}_evr-{}_rvr-{}.csv'
    ensmbl_file = ensmbl_name.format(test_type, index_case, e_screen_interval,
                r_screen_interval, e_vacc_ratio, r_vacc_ratio)
    ensmbl = pd.read_csv(join(no_test_ensmbl_dst, ensmbl_file))
    
    ensmbl['vaccination_scenario'] = vacc_strat_label
    ensmbl['testing_scenario'] = testing_strat_label + '\n' + \
        test_name_map[test_type]
    ensmbl['index_case'] = index_case
    ensmbl['infected_residents'] = ensmbl['E_resident'] + \
                                   ensmbl['I_resident'] + ensmbl['R_resident']
    return ensmbl


def load_TTI_ensmbl(p, TTI_ensmbl_dst, variant='wild_type'):
    test_type, index_case, e_screen_interval, r_screen_interval,\
        testing_strat_label, e_vacc_ratio, r_vacc_ratio, vacc_strat_label = p
    
    TTI_ensmbl_dst = join(TTI_ensmbl_dst, variant)
    
    ensmbl_name = 'test-{}_index-{}_esi-{}_rsi-{}_evr-{}_rvr-{}.csv'
    ensmbl_file = ensmbl_name.format(test_type, index_case, e_screen_interval,
                r_screen_interval, e_vacc_ratio, r_vacc_ratio)
    ensmbl = pd.read_csv(join(TTI_ensmbl_dst, ensmbl_file))
    
    ensmbl['vaccination_scenario'] = vacc_strat_label
    ensmbl['testing_scenario'] = testing_strat_label + '\n'
    ensmbl['index_case'] = index_case
    ensmbl['infected_residents'] = ensmbl['E_resident'] + \
                                   ensmbl['I_resident'] + ensmbl['R_resident']
    return ensmbl


def get_testing_scenario_data(df, testing_scenario):
    tmp = df[df['testing_scenario'] == testing_scenario]
    agg = tmp[['run', 'vaccination_scenario', 'index_case', 'infected_residents']]\
        .groupby(['run', 'vaccination_scenario', 'index_case'])\
        .agg('max')\
        .reset_index()\
        .drop(columns=['run'])
    agg.loc[agg['index_case'] == 'resident', 'infected_residents'] -= 1

    return agg

def get_vaccination_ratio_data(df, vacc_ratio, agent_group='employee'):
    df = df.copy()
    vacc_ratio_col = '{}_vaccination_ratio'.format(agent_group)
    tmp = df[df[vacc_ratio_col] == vacc_ratio]
    agg = tmp[['run', vacc_ratio_col, 'index_case', 'infected_residents']]\
        .groupby(['run', vacc_ratio_col, 'index_case'])\
        .agg('max')\
        .reset_index()\
        .drop(columns=['run'])

    agg.loc[agg['index_case'] == 'resident', 'infected_residents'] -= 1
    return agg


def plot_testing_strategy_vaccination_scenario_grid(data, metric,
        screening_params, vacc_scenarios, vacc_scenario_labels,
        language, variant='', vmax=10, vmin=0, vstep=2):
    
    test_name_map = label_map[language]['test_name_map']
    index_case_map = label_map[language]['index_case_map']
    metric_name_map = label_map[language]['metric_name_map']
    xlabel = label_map[language]['xlabels']['testing_and_vaccination']
    ylabel = label_map[language]['ylabels']['testing_and_vaccination']

    cmap = plt.get_cmap('YlOrRd')
    
    # figure layout & axis setup
    fig = plt.figure(figsize=(15, 20))
    gs = fig.add_gridspec(nrows=8, ncols=5, width_ratios=[1,1,1,1, 0.1],\
                    height_ratios=[0.1,1,0.1,1,0.1,1,0.1,1], wspace=0.05, hspace=0)

    title_ax_1 = fig.add_subplot(gs[0, 0:])
    hmap_ax_1 = fig.add_subplot(gs[1, 0])
    hmap_ax_2 = fig.add_subplot(gs[1, 1])
    hmap_ax_3 = fig.add_subplot(gs[1, 2])
    hmap_ax_4 = fig.add_subplot(gs[1, 3])

    title_ax_2 = fig.add_subplot(gs[2, 0:])
    hmap_ax_5 = fig.add_subplot(gs[3, 0])
    hmap_ax_6 = fig.add_subplot(gs[3, 1])
    hmap_ax_7 = fig.add_subplot(gs[3, 2])
    hmap_ax_8 = fig.add_subplot(gs[3, 3])

    title_ax_3 = fig.add_subplot(gs[4, 0:])
    hmap_ax_9 = fig.add_subplot(gs[5, 0])
    hmap_ax_10 = fig.add_subplot(gs[5, 1])
    hmap_ax_11 = fig.add_subplot(gs[5, 2])
    hmap_ax_12 = fig.add_subplot(gs[5, 3])

    title_ax_4 = fig.add_subplot(gs[6, 0:])
    hmap_ax_13 = fig.add_subplot(gs[7, 0])
    hmap_ax_14 = fig.add_subplot(gs[7, 1])
    hmap_ax_15 = fig.add_subplot(gs[7, 2])
    hmap_ax_16 = fig.add_subplot(gs[7, 3])

    cbar_ax = fig.add_subplot(gs[3:7, 4])

    hmap_axes = [[hmap_ax_1, hmap_ax_2, hmap_ax_3, hmap_ax_4],
                 [hmap_ax_5, hmap_ax_6, hmap_ax_7, hmap_ax_8],
                 [hmap_ax_9, hmap_ax_10, hmap_ax_11, hmap_ax_12],
                 [hmap_ax_13, hmap_ax_14, hmap_ax_15, hmap_ax_16]]

    title_axes = [title_ax_1, title_ax_2, title_ax_3, title_ax_4]

    # compare scenarios in which either employees or residents are the index case
    for i, index_case_mode, test_type in zip([0, 1, 2, 3],
            ['employee', 'resident', 'employee', 'resident'],
            ['same_day_antigen','same_day_antigen','same_day_PCR','same_day_PCR']):

        df = data.loc[test_type, index_case_mode, :, :]

        # remove all axis labels and ticks for the heatmaps
        t_ax = title_axes[i]
        t_ax.set_xticks([])
        t_ax.set_yticks([])
        t_ax.set_frame_on(False)
        t_ax.set_xlim(0, 1)
        t_ax.set_ylim(0, 3)
        t_ax.text(0.23, 1, 'Test: {}, '.format(test_name_map[test_type]) +\
                  index_case_map[index_case_mode], fontsize=20)

        j = 0
        for j,ax,vacc_scenario,vacc_scenario_label in zip(range(4),hmap_axes[i],
                                    vacc_scenarios, vacc_scenario_labels):

            # set flag to set axis ticks only for heatmaps at the boundaries of 
            # the figure
            xticks = False
            yticks = False
            if i > 2:
                xticks = True
            if j in [0, 4]:
                yticks = True

            # put the testing technology in the heatmap title
            ax.set_title(vacc_scenario_label.replace('\n', ' '), fontsize=14)

            # plot heatmap of the scenario
            img = get_image(df.loc[vacc_scenario[0]], vacc_scenario[1],
                            screening_params, metric)
            if index_case_mode == 'resident':
                # if a resident is the index case, we need to subtract 1 from the
                # number of infected residents, to calculate the "outbreak size",
                # which is defined as the number of FOLLOW-UP cases, given an index
                # case
                img = img - 1
            img_plot = plot_heatmap(ax, img, screening_params, vmin, vmax, xticks,
                                    yticks, xlabel, ylabel, cmap=cmap)

            # annotate heatmap with tests / days / agent
            test_rate = get_image(df.loc[vacc_scenario[0]],
                    vacc_scenario[1], screening_params, 'test_rate_mean')

            annotate_heatmap(ax, test_rate)

    # colorbar
    norm = mpl.colors.Normalize(vmin=vmin,vmax=vmax)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    cbar = fig.colorbar(sm, cax=cbar_ax, orientation='vertical',\
                            ticks=np.arange(vmin, vmax + 1, vstep))
    yticklabels = list(range(vmin, vmax, vstep)) + ['$\geq {}$'.format(vmax)]
    cbar.ax.set_yticklabels(yticklabels)
    cbar.set_label('{}'.format(metric_name_map[metric]), fontsize=12)     

    fig.text(0.13, 0.848, 'A', fontweight='bold', fontsize=20)
    fig.text(0.32, 0.848, 'B', fontweight='bold', fontsize=20)
    fig.text(0.51, 0.848, 'C', fontweight='bold', fontsize=20)
    fig.text(0.70, 0.848, 'D', fontweight='bold', fontsize=20)
    
    fig.text(0.13, 0.659, 'E', fontweight='bold', fontsize=20)
    fig.text(0.32, 0.659, 'F', fontweight='bold', fontsize=20)
    fig.text(0.51, 0.659, 'G', fontweight='bold', fontsize=20)
    fig.text(0.70, 0.659, 'H', fontweight='bold', fontsize=20)
    
    fig.text(0.13, 0.471, 'I', fontweight='bold', fontsize=20)
    fig.text(0.32, 0.471, 'J', fontweight='bold', fontsize=20)
    fig.text(0.51, 0.471, 'K', fontweight='bold', fontsize=20)
    fig.text(0.70, 0.471, 'L', fontweight='bold', fontsize=20)
    
    fig.text(0.13, 0.281, 'M', fontweight='bold', fontsize=20)
    fig.text(0.32, 0.281, 'N', fontweight='bold', fontsize=20)
    fig.text(0.51, 0.281, 'O', fontweight='bold', fontsize=20)
    fig.text(0.70, 0.281, 'P', fontweight='bold', fontsize=20)

    plt.savefig('../plots/testing_strategy_and_vaccinations{}_{}.png'\
                     .format(variant, language[0:3]), dpi=300, transparent=True)
    plt.savefig('../plots/testing_strategy_and_vaccinations{}_{}.pdf'\
                     .format(variant, language[0:3]), transparent=True)
    
    
def plot_violins(data, metric, testing_scenarios, vaccination_scenarios,
                 language, variant='', ymin=-3, ymax=20):
    
    metric_name_map = label_map[language]['metric_name_map']
    index_case_map = label_map[language]['index_case_map']
    index_case_label = label_map[language]['index_case_label']
    employee_label = label_map[language]['employee_label']
    resident_label = label_map[language]['resident_label']
    
    fig = plt.figure(figsize=(10, 7.5))
    gs = gridspec.GridSpec(2, 2)
    gs.update(wspace=0.1, hspace=0.4)

    for i, testing_scenario in enumerate(testing_scenarios):
        ax = plt.subplot(gs[int(i/2), i%2])
        agg = get_testing_scenario_data(data, testing_scenario)
        agg = agg[agg['vaccination_scenario'].isin(vaccination_scenarios)]
        sns.violinplot(x='vaccination_scenario', y=metric, 
                  hue='index_case', data=agg, split=True, ax=ax,
                  palette=['FireBrick', 'DarkBlue'], order=vaccination_scenarios)

        ax.set_xlabel('')
        ax.text(2.4, ymax*0.8, testing_scenario, ha='center',
               fontsize=10)
        ax.set_ylim(ymin, ymax)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)

        if i in [1, 3]:
            ax.set_ylabel('')
            ax.set_yticks([])
            ax.spines['left'].set_visible(False)
            if i == 1:
                l = ax.legend(title=index_case_label, loc=2, 
                        bbox_to_anchor=[-0.3, -0.3, 1, 1], fontsize=10)
                l.get_title().set_fontsize(10)
                l.get_texts()[0].set_text(employee_label)
                l.get_texts()[1].set_text(resident_label)
            else:
                ax.legend([], [], frameon=False)
        else:
            ax.set_ylabel(metric_name_map[metric])
            ax.legend([], [], frameon=False)
            
        for violin in ax.collections[::]:
            violin.set_alpha(0.7)

    fig.text(0.135, 0.87, 'A', fontweight='bold', fontsize=14)
    fig.text(0.535, 0.87, 'B', fontweight='bold', fontsize=14)
    fig.text(0.135, 0.425, 'C', fontweight='bold', fontsize=14)
    fig.text(0.535, 0.425, 'D', fontweight='bold', fontsize=14)
    plt.savefig('../plots/vaccination_scenario_violins{}_{}.pdf'\
                .format(variant, language[0:3]), transparent=True)
    plt.savefig('../plots/vaccination_scenario_violins{}_{}.png'\
                .format(variant, language[0:3]), dpi=300)
    
    
def plot_violins2(data, metric, testing_scenarios, vaccination_scenarios,
                  language, variant='', ymin=-3, ymax=20):
    
    metric_name_map = label_map[language]['metric_name_map']
    index_case_map = label_map[language]['index_case_map']
    
    fig = plt.figure(figsize=(15, 4))
    gs = gridspec.GridSpec(1, 3)
    gs.update(wspace=0.1, hspace=0.4)

    for i, testing_scenario in enumerate(testing_scenarios):
        ax = plt.subplot(gs[i])
        agg = get_testing_scenario_data(data, testing_scenario)
        agg = agg[agg['vaccination_scenario'].isin(vaccination_scenarios)]
        sns.violinplot(x='vaccination_scenario', y=metric, 
                  hue='index_case', data=agg, split=True, ax=ax,
                  palette=['FireBrick', 'DarkBlue'], order=vaccination_scenarios)

        ax.set_xlabel('')
        ax.text(1.5, ymax*0.8, testing_scenario, ha='center',
               fontsize=11)
        ax.set_ylim(ymin, ymax)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)

        if i in [1, 2]:
            ax.set_ylabel('')
            #ax.set_yticks([])
            #ax.spines['left'].set_visible(False)
            if i == 1:
                l = ax.legend(title='', loc=5,
                    bbox_to_anchor=[0, 0.1, 1, 1], fontsize=10)
                l.get_title().set_fontsize(10)
                l.get_texts()[0].set_text(index_case_map['employee'])
                l.get_texts()[1].set_text(index_case_map['resident'])
            else:
                ax.legend([], [], frameon=False)
        else:
            ax.set_ylabel(metric_name_map[metric])
            ax.legend([], [], frameon=False)
            
        for violin in ax.collections[::]:
            violin.set_alpha(0.7)

    fig.text(0.135, 0.82, 'A', fontweight='bold', fontsize=14)
    fig.text(0.4, 0.82, 'B', fontweight='bold', fontsize=14)
    fig.text(0.67, 0.82, 'C', fontweight='bold', fontsize=14)
    
    plt.savefig('../plots/vaccination_scenario_violins2{}_{}.pdf'\
                .format(variant, language[0:3]), transparent=True)
    plt.savefig('../plots/vaccination_scenario_violins2{}_{}.png'\
                .format(variant, language[0:3]), dpi=300)
    
    
def plot_vaccination_heatmap(data, metric, vaccination_ratios, language, 
                             variant=''):

    metric_name_map = label_map[language]['metric_name_map']
    index_case_map = label_map[language]['index_case_map']
    xlabel = label_map[language]['xlabels']['vaccination']
    ylabel = label_map[language]['ylabels']['vaccination']

    cmap = plt.get_cmap('YlOrRd')
    
    # figure layout & axis setup
    fig, axes = plt.subplots(1, 2, figsize=(15, 9))

    vmin=0
    vmax=10
    vstep=1

    # compare scenarios in which either employees or residents are the
    # index case
    for i, index_case_mode, ax in zip([0, 1], ['employee', 'resident'], axes):

        ax.set_title(index_case_map[index_case_mode], fontsize=20)

        # set flag to set axis ticks only for heatmaps at the boundaries of 
        # the figure
        xticks = True
        yticks = False
        if i == 0:
            yticks = True

        # plot heatmap of the scenario
        img = get_image(data, index_case_mode, vaccination_ratios, metric)
        if index_case_mode == 'resident':
            # if a resident is the index case, we need to subtract 1 from the
            # number of infected residents, to calculate the "outbreak size",
            # which is defined as the number of FOLLOW-UP cases, given an index
            # case
            img = img - 1
        img_plot = plot_heatmap(ax, img,
                ['{:1.0f}%'.format(100*i) for i in vaccination_ratios],
                vmin, vmax, xticks, yticks, xlabel, ylabel, cmap=cmap)

        mask = np.where(img < 1, 1, 0)
        for i in range(img.shape[0]):
            for j in range(img.shape[1]):
                if mask[i, j]:
                    rect = plt.Rectangle((j-0.5, i-0.5), 1, 1, color='g',
                                    alpha=1, fill=False, ec='#D6DBDF', lw=4)
                    ax.add_patch(rect)

        ax.xaxis.label.set_size(18)
        ax.yaxis.label.set_size(20)

    # colorbar
    divider = make_axes_locatable(axes[1])
    cbar_ax = divider.append_axes('right', size='4%', pad=0.1)

    norm = mpl.colors.Normalize(vmin=vmin,vmax=vmax)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    cbar = fig.colorbar(sm, cax=cbar_ax, orientation='vertical',\
                            ticks=np.arange(vmin, vmax + 1, vstep))
    yticklabels = list(range(vmin, vmax, vstep)) + ['$\geq {}$'.format(vmax)]
    cbar.ax.set_yticklabels(yticklabels)
    cbar.set_label('{}'.format(metric_name_map[metric]), fontsize=18)      

    # dummy axis to preserve spacing
    divider = make_axes_locatable(axes[0])
    cbar_ax = divider.append_axes('right', size='4%', pad=0.0)
    cbar_ax.set_axis_off()

    fig.text(0.07, 0.8525, 'A', fontweight='bold', fontsize=22)
    fig.text(0.512, 0.8525, 'B', fontweight='bold', fontsize=22)

    plt.tight_layout()

    plt.savefig('../plots/vaccinations{}_{}.pdf'.format(variant, language[0:3]))
    plt.savefig('../plots/vaccinations{}_{}.png'.format(variant, language[0:3]))


def plot_testing_strategy_heatmaps(data, metric, test_types, sim_name, 
    screening_params, cmap, vmin, vmax, vstep, language, variant,
    plt_dst, fname_addition=''):
    
    index_case_map = label_map[language]['index_case_map']
    test_name_map = label_map[language]['test_name_map']
    xlabel = label_map[language]['xlabels']['testing_strategy']
    ylabel = label_map[language]['ylabels']['testing_strategy']
    metric_name_map = label_map[language]['metric_name_map']
    frequency_name_map = label_map[language]['frequency_name_map']

    # figure layout & axis setup
    fig = plt.figure(figsize=(15, 12))
    gs = fig.add_gridspec(nrows=4, ncols=4, width_ratios=[1,1,1, 0.05],\
                           height_ratios=[0.1,1,0.1,1], wspace=0.05, hspace=0)

    title_ax_1 = fig.add_subplot(gs[0, 0:])
    hmap_ax_1 = fig.add_subplot(gs[1, 0])
    hmap_ax_2 = fig.add_subplot(gs[1, 1])
    hmap_ax_3 = fig.add_subplot(gs[1, 2])

    title_ax_2 = fig.add_subplot(gs[2, 0:])
    hmap_ax_4 = fig.add_subplot(gs[3, 0])
    hmap_ax_5 = fig.add_subplot(gs[3, 1])
    hmap_ax_6 = fig.add_subplot(gs[3, 2])

    cbar_ax = fig.add_subplot(gs[1:, 3])

    hmap_axes = [[hmap_ax_1, hmap_ax_2, hmap_ax_3], [hmap_ax_4, hmap_ax_5, hmap_ax_6]]
    title_axes = [title_ax_1, title_ax_2]

    # compare scenarios in which either employees or residents are the index case
    for i, index_case_mode in enumerate(['employee', 'resident']):
        df = data.loc[:,:,:, index_case_mode]
        
        # set flag to set axis ticks only for heatmaps at the boundaries of 
        # the figure
        t_ax = title_axes[i]
        t_ax.set_xticks([])
        t_ax.set_yticks([])
        t_ax.set_frame_on(False)
        t_ax.set_xlim(0, 1)
        t_ax.set_ylim(0, 3)
        t_ax.text(0.38, 1, index_case_map[index_case_mode], fontsize=20)
        
        # compare different test result turnover times for PCR tests
        for j, ax, test_type in zip(range(3), hmap_axes[i], test_types):
            xticks = False
            yticks = False
            if i > 0:
                xticks = True
            if j in [0, 3]:
                yticks = True
                
            # put the turnover time in the heatmap title
            ax.set_title('Test: {}'\
                    .format(test_name_map[test_type]), fontsize=14)
            
            # plot heatmap of the scenario
            img = get_image(df, test_type, screening_params, metric)
            if index_case_mode == 'resident' and metric in \
                ['infected_residents_mean', 'infected_residents_median',
                 'infected_residents_0.90', 'infected_residents']:
                img = img - 1
            img_plot = plot_heatmap(ax, img, screening_params, vmin, vmax,
                    xticks, yticks, xlabel, ylabel, cmap=cmap)
            
            # annotate heatmap with tests / days / agent
            #test_rate = get_image(df, test_type, screening_params, 
            #                            'test_rate_mean')
            # annotate_heatmap(ax, test_rate)

    # colorbar
    norm = mpl.colors.Normalize(vmin=vmin,vmax=vmax)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    cbar = fig.colorbar(sm, cax=cbar_ax, orientation='vertical',\
                    ticks=np.arange(vmin, vmax+vstep, vstep))
    yticklabels = ['{:1.2f}'.format(i) for i in np.arange(vmin, vmax, vstep)] + \
            ['$\geq {}$'.format(vmax)]
    cbar.ax.set_yticklabels(yticklabels)
    cbar.set_label('{}'.format(metric_name_map[metric]), fontsize=12)    

    # saving of plots
    plt.savefig(join(plt_dst, '{}{}{}_{}.png'\
        .format(sim_name, variant, fname_addition, language[0:3])),
                dpi=300, transparent=True)
    plt.savefig(join(plt_dst,'{}{}{}_{}.pdf'\
        .format(sim_name, variant, fname_addition, language[0:3])))