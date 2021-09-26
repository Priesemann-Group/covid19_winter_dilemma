import numpy as np
import sys
from multiprocessing import Pool
import argparse
import pickle

sys.path.append("../code")
sys.path.append("../notebooks")
import parametros
import model


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


params = parametros.params.copy()


# Input Parameter: scenarios
scenarios = ["scenario1", "scenario2", "scenario3"]

# Input Parameter: alphas
alphas = parametros.alpharange


def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i : i + n]


def run(mapping):
    """
    Runs model and saves the data
    """
    scenario, alphaR, alphau, alphaw = mapping

    params["Rt_base"] = parametros.Rtbase[scenario]

    m = model.Model(**params)
    # comment

    m.alpha_R = alphaR
    m.alpha_u = alphau
    m.alpha_w = alphaw

    times1, data1 = m.run()

    with open(
        f"../datamodelruns/sweep/scen={scenario}-aR={m.alpha_R}-au={m.alpha_u}-aw={m.alpha_w}.pickle",
        "wb",
    ) as f:
        pickle.dump(m.data, f)


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
    for scen in scenarios:
        for alphR in alphas:
            for alphu in alphas:
                for alphw in alphas:
                    mapping.append([scen, alphR, alphu, alphw])

    # Since we have 32 cores available let's start multiple models multithreaded
    mapping_chunks = list(chunks(mapping, 32))  # Create 32 mappings

    pool = Pool()
    pool.map(run, mapping_chunks[i])
    pool.close()
