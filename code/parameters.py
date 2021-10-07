import numpy as np


def get_params(country='Germany', modeltype='base', scenario='scenario2'):
    mapping_raw = {'Germany':raw_deu, 'Czech':raw_cze, 'Denmark':raw_dsk, 'Portugal':raw_por}
    mapping_params = {'Germany':{}, 'Czech':params_cze, 'Denmark':params_dsk, 'Portugal':params_por}
    params = params_base.copy()
    params.update(mapping_params[country])

    mapping_modelparams = {'base':{}, 'ramp':params_ramp, 'logistic':params_logistic}
    params.update(mapping_modelparams[modeltype])

    params['Rt_base'] = Rtbase[scenario]

    raw = mapping_raw[country]
    return calc_y0(params, **raw)

def calc_y0(params, V_raw, R_raw, cases, ICU, W_V, W_R, serop):
    darkfigure = serop/(R_raw/1e6)

    y0 = {
        'R': darkfigure * (R_raw-W_R),
        'W': W_V + darkfigure*W_R,
        'ICU': ICU,
    }
    y0['V']= V_raw - W_V - y0['R']*(1-(V_raw-W_V)/1e6)

    E_stay = 1./params['rho']
    I_stay = 1./(params['gamma']+params['delta']+params['Theta'])

    immune_infected = (y0['R']+y0['V'])*(1-params['eta'])
    S_approx = 1e6 - sum(list(y0.values())) - darkfigure*(E_stay+I_stay)*cases
    breakthrough = immune_infected / (immune_infected + S_approx)

    y0.update({
        'E': darkfigure*cases*E_stay*(1-breakthrough),
        'EB': darkfigure*cases*E_stay*breakthrough,
        'I': darkfigure*cases*I_stay*(1-breakthrough),
        'IB': darkfigure*cases*I_stay*breakthrough,
    })

    y0['S'] = 1e6 - sum(list(y0.values()))

    y0.update({
        'UC': V_raw,
        'WC': 0.,
        'D': 0.,
        'C': 0.,
    })

    y0_array = [y0['S'],y0['V'],y0['W'],y0['E'],y0['EB'],y0['I'],y0['IB'],y0['ICU'],y0['R'],y0['UC'],y0['WC'],y0['D'],y0['C']]
    params['y0'] = y0_array

    return params


# Values from OWD
raw_deu = {
    'V_raw': 603900.0,
    'R_raw':  47400.0,
    'cases':    116.25,
    'ICU':       13.46,
    'W_V':    18000.0,
    'W_R':     9030.0,
    'serop':      0.112,
}

raw_cze = {
    'V_raw': 536500.0,
    'R_raw': 156600.0,
    'cases':     18.56,
    'ICU':        0.84,
    'W_V':    14650.0,
    'W_R':    29100.0,
    'serop':      0.51, #Feb and March 2021; 15. March reported = 0.13
}

raw_dsk = {
    'V_raw': 725400.0,
    'R_raw':  59800.0,
    'cases':    147.03,
    'ICU':        4.3,
    'W_V':    22450.0,
    'W_R':    11610.0,
    'serop':      0.173, #Dec 2020 serop=0.067, 30. Dec reported=0.0278
}

raw_por = {
    'V_raw': 755500.0,
    'R_raw': 102200.0,
    'cases':    191.4,
    'ICU':       12.88,
    'W_V':    18550.0,
    'W_R':    18400.0,
    'serop':      0.173,
}


params_base = {
    'Rt_base': 3.25,
    'Rt_free': 5,
    'eta': 0.75,
    'kappa': 0.8,
    'sigma': 0.5,
    'gamma': 0.1,
    'gamma_ICU': 0.13,
    'delta': 0.0019,
    'rho': 0.25,
    'omega_v_b': 1./(8*30),
    'omega_n_b': 1./(12*30),
    'chi_0': 0.1,
    'chi_1': 0.2,
    'alpha_w': 0.03,
    'alpha_u': 0.02,
    'alpha_R': 0.01,
    'e_R': 0.,
    'e_u': 0.,
    'e_w': 0.,
    'Phi_0': 0.0025,
    'phi_0': 0.0025,
    'u_base': 0.5,
    'mu': 0.267,
    'd_0': 8*30.,
    'd_mu': 0.,
    'a_Rt': 4.,
    'a_vac': 14.,
    'gamma_cutoff': 30.,
    'tau_vac1': 6*7.,
    'tau_vac2': 2*7.,
    'Theta':0.0005366,
    'Theta_ICU':0.09755,
    't_max': 360.,
    'step_size': 0.1,
    'influx': 1,
    'epsilon_free':14.,
}

params_cze = {
    'delta': 0.001717,
    'Theta':0.00040588,
    'Theta_ICU':0.0978747,
}

params_dsk = {
    'delta': 0.0017,
    'Theta':0.000429676,
    'Theta_ICU':0.09786686,
}

params_por = {
    'delta': 0.0019083,
    'Theta':0.000523277,
    'Theta_ICU':0.097566211,
}

params_ramp = {
    'epsilon_u':0.1,
    'epsilon_w':0.1,
}

params_logistic = {
    'r_base':0.015,
}


Rtbase = {
    'scenario1':5,
    'scenario2':3.75,
    'scenario3':2.5,
}


#Sweep range for the alphas:
alpharange=np.round(np.linspace(1/50, 1/400, 100), decimals=4)