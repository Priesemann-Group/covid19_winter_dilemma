import numpy as np
import sys
from multiprocessing import Pool
import argparse
import pickle

sys.path.append("../code")
sys.path.append("../notebooks")
import parametros
import model_ramp as model

import numpy.random as npr


params = parametros.params.copy()
params.update({"y0": parametros.y0_array})


# Input Parameter: scenarios
scenarios = ["scenario1", "scenario2", "scenario3"]

nrValues = 1000
# Input Parameter: alphas drawn from their distributions
alphasR = np.round(np.linspace(0, 2 * params["alpha_R"], nrValues), decimals=6)

alphasu = np.round(
    np.linspace(
        params["alpha_u"] - 1 / 2 * params["alpha_u"],
        params["alpha_u"] + 1 / 2 * params["alpha_u"],
        nrValues,
    ),
    decimals=6,
)
alphasw = np.round(
    np.linspace(
        params["alpha_w"] - 1 / 3 * params["alpha_w"],
        params["alpha_w"] + 1 / 3 * params["alpha_w"],
        nrValues,
    ),
    decimals=6,
)


def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i : i + n]


def run(mapping):
    """
    Runs model and saves the data
    """
    scen, aR, au, aw = mapping

    params["Rt_base"] = parametros.Rtbase[scen]

    m = model.Model(**params)
    # comment

    m.alpha_R = aR
    m.alpha_u = au
    m.alpha_w = aw

    times1, data1 = m.run()

    fstring = f"/scratch03.local/smohr/covid19_wd_sweeps/sens_sweep-scen={scen}-aR={m.alpha_R}-au={m.alpha_u}-aw={m.alpha_w}.npz"
    np.savez_compressed(fstring, np.array(m.chopped_data()))


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

    u = params["alpha_u"]
    w = params["alpha_w"]
    R = params["alpha_R"]
    # Create mapping
    # TO DO: Put the things into mapping
    mapping = []
    for sc in scenarios:
        for i in range(nrValues):
            mapping.append([sc, alphasR[i], u, w])
        for i in range(nrValues):
            mapping.append([sc, R, alphasu[i], w])
        for i in range(nrValues):
            mapping.append([sc, R, u, alphasw[i]])

    # Since we have 32 cores available let's start multiple models multithreaded
    mapping_chunks = list(chunks(mapping, 32))  # Create 32 mappings

    pool = Pool()
    pool.map(run, mapping_chunks[args.id - 1])
    pool.close()
