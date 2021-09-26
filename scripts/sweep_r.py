import numpy as np

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


params = parametros.params.copy()


#Input Parameter: scenarios
scenarios=['scenario1', 'scenario2', 'scenario3']

# Input Parameter: alphas
alphas = parametros.alpharange

mapping=[]
for scen in scenarios:
    for alphR in alphas:
        for alphu in alphas:
            for alphw in alphas:
                mapping.append([scen, alphR, alphu, alphw])



parser = argparse.ArgumentParser(description='Sweeps')
parser.add_argument(
    "-i", "--id", type=int, help="ID", required=True,
)


args = parser.parse_args()
scenario, alphaR, alphau, alphaw = mapping[args.id-1]
print(args.id)



params["Rt_base"]=parametros.Rtbase[scenario]

m = model.Model(**params)


m.alpha_R = alphaR
m.alpha_u = alphau
m.alpha_w = alphaw

times1, data1 = m.run()

with open(f"../datamodelruns/sweep/scen={scenario}-aR={m.alpha_R}-au={m.alpha_u}-aw={m.alpha_w}.pickle", "wb") as f:
        pickle.dump(m.data, f)
