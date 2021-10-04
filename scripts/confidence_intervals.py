import numpy as np
import sys
from multiprocessing import Pool
import argparse
import pickle

sys.path.append("../code")
sys.path.append("../notebooks")
import parametros
import model


params = parametros.params.copy()

# Input Parameter: scenarios
scenarios = ["scenario1", "scenario2", "scenario3"]


def run(mapping):
    """
    Runs model and saves the data
    """
    scenario, alphaR, alphau, alphaw, i = mapping

    params["Rt_base"] = parametros.Rtbase[scenario]

    m = model.Model(**params)
    # comment

    m.alpha_R = alphaR
    m.alpha_u = alphau
    m.alpha_w = alphaw

    times1, data1 = m.run()

    fstring = f"/scratch03.local/smohr/covid19_wd_sweeps/run={3*(args.id-1)+i}-scen={scenario}.npz"
    np.savez_compressed(fstring,np.array(m.chopped_data()[:,7]))


if __name__ == "__main__":
    # Parse arguments
    parser = argparse.ArgumentParser(description="Sweeps")
    parser.add_argument(
        "-i",
        "--id",
        type=int,
        help="ID",
        required=True,
    )
    args = parser.parse_args()

    # Create mapping
    mapping = []

    alpha_w = np.random.default_rng().normal(loc=0.018, scale=0.018/3., size=10)
    alpha_u = np.random.default_rng().normal(loc=0.014, scale=0.014/3., size=10)
    alpha_R = np.random.default_rng().normal(loc=0.007, scale=0.007/3., size=10)
    alpha_w *= (alpha_w >= 0)
    alpha_u *= (alpha_w >= 0)
    alpha_R *= (alpha_w >= 0)

    for scen in scenarios:
        for i in range(10):
            mapping.append([scen, alpha_w[i], alpha_u[i], alpha_R[i], i])

    pool = Pool()
    pool.map(run, mapping)
    pool.close()
