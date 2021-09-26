import numpy as np

y0 = {
    "S": 276618.0,
    "V": 599864.0,
    "W": 97500.0,
    "E": 390.0,
    "EB": 39.0,
    "I": 974.0,
    "IB": 97.0,
    "ICU": 18.0,
    "R": 24500.0,
    "UC": 650000.0,
    "WC": 5000.0,
}

params = {
    "y0": list(y0.values()),
    "Rt_base": 3.25,
    "Rt_free": 4.5,
    "eta": 0.8,
    "kappa": 0.95,
    "sigma": 0.5,
    "gamma": 0.1,
    "gamma_ICU": 0.13,
    "delta": 0.0019,
    "rho": 0.25,
    "omega_v_b": 1.0 / (6 * 30),
    "omega_n_b": 1.0 / (12 * 30),
    "chi_0": 0.1,
    "chi_1": 0.2,
    "alpha_w": 0.018,
    "alpha_u": 0.008,
    "alpha_R": 0.01,
    "e_R": 0.0,
    "e_u": 0.0,
    "e_w": 0.0,
    "Phi_0": 0.0025,
    "phi_0": 0.0025,
    "u_base": 0.5,
    "mu": 0.267,
    "d_0": 8 * 30.0,
    "d_mu": 0.0,
    "a_Rt": 4.0,
    "a_vac": 14.0,
    "gamma_cutoff": 30.0,
    "tau_vac1": 6 * 7.0,
    "tau_vac2": 2 * 7.0,
    "t_max": 360.0,
    "step_size": 0.1,
}


Rtbase = {
    "scenario1": 4.5,
    "scenario2": 3.25,
    "scenario3": 2,
}


# Sweep range for the alphas:
alpharange = np.round(np.linspace(1 / 50, 1 / 400, 30), decimals=4)
