import numpy as np
import sys
from multiprocessing import Pool
import argparse
import pickle

sys.path.append("../code")
sys.path.append("../notebooks")
import parametros
import model_ramp as model
from scipy.stats import norm

params = parametros.params.copy()
params.update({"y0": parametros.y0_array})


# Input Parameter: scenarios
scenarios = ["scenario1", "scenario2", "scenario3"]

nrValues = 33000
sigma = 0.01 / 3

# Input Parameter: alphas drawn from their distributions
alphas = {}
for alpha in ["alpha_R", "alpha_u", "alpha_w"]:
    a = params[alpha]
    l = a - 4 * sigma
    u = a + 4 * sigma
    x = np.linspace(l, u, int(np.floor(nrValues ** (1 / 3))))
    alphas[alpha] = np.round(norm.pdf(x, loc=a, scale=sigma), decimals=6)


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

    fstring = f"/scratch03.local/smohr/covid19_wd_sweeps/sens_dist-scen={scen}-aR={m.alpha_R}-au={m.alpha_u}-aw={m.alpha_w}.npz"
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

    # Create mapping
    mapping = []
    for sc in scenarios:
        for aR in alphas["alpha_R"]:
            for au in alphas["alpha_u"]:
                for aw in alphas["alpha_w"]:
                    mapping.append([sc, aR, au, aw])

    # Since we have 32 cores available let's start multiple models multithreaded
    mapping_chunks = list(chunks(mapping, 32))  # Create 32 mappings

    pool = Pool()
    pool.map(run, mapping_chunks[args.id - 1])
    pool.close()
