import numpy as np

y0 = {
    'S': 276618.,
    'V': 599864.,
    'W':  97500.,
    'E':    390.,
    'EB':    39.,
    'I':    974.,
    'IB':    97.,
    'ICU':   18.,
    'R':  24500.,
    'UC':650000.,
    'WC':  5000.,
    'D': 0.,
    'C': 0.,
}

params = {
    'y0': list(y0.values()),
    'Rt_base': 3.25,
    'Rt_free': 5,
    'eta': 0.8,
    'kappa': 0.95,
    'sigma': 0.5,
    'gamma': 0.1,
    'gamma_ICU': 0.13,
    'delta': 0.0019,
    'rho': 0.25,
    'omega_v_b': 1./(6*30),
    'omega_n_b': 1./(12*30),
    'chi_0': 0.1,
    'chi_1': 0.2,
    'alpha_w': 0.021,
    'alpha_u': 0.014,
    'alpha_R': 0.007,
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
}


Rtbase = {
    'scenario1':5,
    'scenario2':3.5,
    'scenario3':2,
}


#Sweep range for the alphas:
alpharange=np.round(np.linspace(1/50, 1/400, 100), decimals=4)