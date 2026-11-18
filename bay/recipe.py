import datetime as dt
import json
import os
import shutil
import time
from concurrent.futures import ProcessPoolExecutor, wait
from multiprocessing import Manager
from pathlib import Path

import numpy as np
import yaml

from .mcmc import run_chain
from .stats import NAMES


def run_one(spec, data_path, out_dir, queue=None):
    """One chain. Status lines go to `queue` if given: a worker's own prints
    never reach the notebook."""
    def report(message):
        (queue.put if queue is not None else print)(f"{spec['name']:24s} {message}")

    report("started")
    start = time.time()
    sampler = run_chain(spec, data_path, report)
    path = Path(out_dir) / f"chain_{spec['name']}.npz"
    np.savez(path,
             chain=sampler.get_chain(thin=spec["thin"]).astype(np.float32),
             log_prob=sampler.get_log_prob(thin=spec["thin"]),
             acc=sampler.acceptance_fraction,
             seconds=time.time() - start,
             names=np.array(NAMES),
             spec=json.dumps(spec),
             data=data_path)
    return str(path)


def run(recipe_path, max_workers=None):
    """Run every chain in the recipe. Returns the chain paths, sorted."""
    recipe_path = Path(recipe_path)
    recipe = yaml.safe_load(recipe_path.read_text(encoding="utf-8"))
    names = [spec["name"] for spec in recipe["chains"]]
    if len(set(names)) != len(names):
        raise ValueError("chain names must be unique")

    stamp = dt.datetime.now().strftime("%Y-%m-%d_%H%M")
    out_dir = recipe_path.parent / "results" / f"{stamp}_{recipe['notes']}"
    out_dir.mkdir(parents=True)
    shutil.copy(recipe_path, out_dir / recipe_path.name)
    data_path = str((recipe_path.parent / recipe["data"]).resolve())

    workers = min(len(names), max_workers or os.cpu_count())
    print(f"{len(names)} chains on {workers} processes -> {out_dir}")
    start = time.time()

    def log(message):
        print(f"[{(time.time() - start) / 60:6.1f} min] {message}", flush=True)

    with Manager() as manager, ProcessPoolExecutor(workers) as pool:
        queue = manager.Queue()
        futures = {pool.submit(run_one, spec, data_path, out_dir, queue): spec["name"]
                   for spec in recipe["chains"]}
        pending = set(futures)
        while pending:
            done, pending = wait(pending, timeout=1)
            while not queue.empty():
                log(queue.get())
            for future in done:                  # .result() re-raises a failed chain
                log(f"{futures[future]:24s} saved {Path(future.result()).name}")
    log("run finished")
    return sorted(future.result() for future in futures)


def load(path):
    """One saved chain: chain, log_prob, acc, names, spec (dict), data (path)."""
    with np.load(path) as f:
        out = {key: f[key] for key in f.files}
    out["spec"] = json.loads(str(out["spec"]))
    out["data"] = str(out["data"])
    return out
