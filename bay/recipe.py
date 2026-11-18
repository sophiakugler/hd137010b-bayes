import datetime as dt
import json
import os
import shutil
import sys
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
        if queue is None:
            print(f"{spec['name']:24s} {message}")
        else:
            queue.put((spec["name"], message))

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


BAR_WIDTH = 24


def _bar(done, total):
    filled = int(BAR_WIDTH * done / total) if total else 0
    return "[" + "#" * filled + "-" * (BAR_WIDTH - filled) + "]"


def _line(name, message, width):
    """One chain, one line. `message` ending in "i/n" becomes a bar."""
    phase, _, counts = message.rpartition(" ")
    done, _, total = counts.partition("/")
    if phase and done.isdigit() and total.isdigit():
        done, total = int(done), int(total)
        pct = 100 * done / total if total else 0
        return f"{name:<{width}}  {phase:<8}  {_bar(done, total)} {pct:3.0f}%  {done}/{total}"
    return f"{name:<{width}}  {message}"


class _Screen:
    """Exactly `n` lines, redrawn in place: a notebook display handle when
    there is a kernel, ANSI cursor-up otherwise."""

    def __init__(self, n):
        self.n, self.handle, self.drawn = n, None, False
        try:
            from IPython import get_ipython
            from IPython.display import display
            if get_ipython() is not None:
                self.handle = display({"text/plain": ""}, raw=True, display_id=True)
        except Exception:
            self.handle = None

    def draw(self, lines):
        if self.handle is not None:
            self.handle.update({"text/plain": "\n".join(lines)}, raw=True)
            return
        if self.drawn:
            sys.stdout.write(f"\x1b[{self.n}A")
        sys.stdout.write("".join(f"\x1b[2K{line}\n" for line in lines))
        sys.stdout.flush()
        self.drawn = True


def run(recipe_path, data_path, output_path, max_workers=None):
    """Run every chain in the recipe. Returns the chain paths, sorted.

    `data_path` is the photometry every chain reads; `output_path` is the
    directory the timestamped run directory is created in. Both are taken
    from the caller, not from the recipe's location.
    """
    recipe_path = Path(recipe_path)
    data_path = Path(data_path).resolve()
    output_path = Path(output_path)
    recipe = yaml.safe_load(recipe_path.read_text(encoding="utf-8"))
    names = [spec["name"] for spec in recipe["chains"]]
    if len(set(names)) != len(names):
        raise ValueError("chain names must be unique")

    stamp = dt.datetime.now().strftime("%Y-%m-%d_%H%M")
    out_dir = output_path / f"{stamp}_{recipe['notes']}"
    out_dir.mkdir(parents=True)
    shutil.copy(recipe_path, out_dir / recipe_path.name)

    workers = min(len(names), max_workers or os.cpu_count())
    print(f"{len(names)} chains on {workers} processes -> {out_dir}")
    start = time.time()

    width = max(len(name) for name in names)
    state = dict.fromkeys(names, "queued")
    screen = _Screen(len(names))

    def draw():
        screen.draw([_line(name, state[name], width) for name in names])

    draw()
    with Manager() as manager, ProcessPoolExecutor(workers) as pool:
        queue = manager.Queue()
        futures = {pool.submit(run_one, spec, str(data_path), out_dir, queue): spec["name"]
                   for spec in recipe["chains"]}
        pending = set(futures)
        while pending:
            done, pending = wait(pending, timeout=0.2)
            dirty = bool(done)
            while not queue.empty():
                name, message = queue.get()
                state[name], dirty = message, True
            for future in done:                  # .result() re-raises a failed chain
                state[futures[future]] = f"saved {Path(future.result()).name}"
            if dirty:
                draw()
    draw()
    print(f"run finished in {(time.time() - start) / 60:.1f} min")
    return sorted(future.result() for future in futures)


def load(path):
    """One saved chain: chain, log_prob, acc, names, spec (dict), data (path)."""
    with np.load(path) as f:
        out = {key: f[key] for key in f.files}
    out["spec"] = json.loads(str(out["spec"]))
    out["data"] = str(out["data"])
    return out
