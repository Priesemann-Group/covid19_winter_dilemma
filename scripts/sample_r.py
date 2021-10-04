import numpy as np
import numpy.random as npr

import sys
sys.path.append("../code")
sys.path.append("../notebooks")
import parametros
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


params = parametros.params


# Input Parameter: alphas
normal = npr.normal(params['alpha_R'], params['alpha_R']/3, 100)
alphas = normal[normal >=0]


parser = argparse.ArgumentParser(description='Samples')
parser.add_argument(
    "-i", "--id", type=int, help="ID", required=True,
)


args = parser.parse_args()
print(args.id)

m = model.Model(**params)

m.alpha_R = alphas[args.id]
times, data = m.run()

with open(f"../datamodelruns/alpha_r/sample/alphaR={m.alpha_R}.pickle", "wb") as f:
        pickle.dump(m, f)
