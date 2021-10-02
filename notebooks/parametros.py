import numpy as np

# Values from OWD
V_raw = 603900.0
R_raw =  47400.0
cases =    116.25
ICU =       13.46
W_V =    18000.0
W_R =     9030.0
darkfigure = 2.5

params = {
#    'y0': y0_array,
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

y0 = {
    'V': V_raw - W_V,
    'R': darkfigure * (R_raw-W_R),
    'W': W_V + darkfigure*W_R,
    'ICU': ICU,
}

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


Rtbase = {
    'scenario1':5,
    'scenario2':3.75,
    'scenario3':2.5,
}



#Sweep range for the alphas:
alpharange=np.round(np.linspace(1/50, 1/400, 100), decimals=4)