import numpy as np

import sys
sys.path.append("../code")
import model

import argparse
import pickle

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
}

params = {
    'y0': list(y0.values()),
    'Rt_base': 3.5,
    'Rt_free': 5.0,
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
    'alpha_w': 0.008,
    'alpha_u': 0.01,
    'alpha_R': 0.008,
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
    't_max': 360.,
    'step_size': 0.1,
}



# Input Parameter: alphas

ICU = np.round(np.linspace(50, 200, 10), decimals=4)
alphas = 1/ICU


parser = argparse.ArgumentParser(description='Sweeps')
parser.add_argument(
    "-i", "--id", type=int, help="ID", required=True,
)


args = parser.parse_args()
print(args.id)

m = model.Model(**params)

m.alpha_R = alphas[args.id]
times, data = m.run()

with open(f"../datamodelruns/alphaR={m.alpha_R}.pickle", "wb") as f:
        pickle.dump(m, f)
