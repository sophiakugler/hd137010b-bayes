from pathlib import Path
import sys,hashlib
import nbformat
from nbclient import NotebookClient
from jupyter_client import KernelManager

ROOT=Path(__file__).resolve().parents[1]
for name in ['04_eccentric_mcmc.ipynb','05_prior_sensitivity.ipynb']:
    path=ROOT/name
    initial_hash=hashlib.sha256(path.read_bytes()).hexdigest()
    nb=nbformat.read(path,as_version=4)
    km=KernelManager(kernel_name='python3')
    km.kernel_spec.argv[0]=sys.executable
    client=NotebookClient(nb,km=km,timeout=300,resources={'metadata':{'path':str(ROOT)}})
    try:
        client.execute()
    finally:
        if km.has_kernel:km.shutdown_kernel(now=True)
    nbformat.validate(nb)
    if hashlib.sha256(path.read_bytes()).hexdigest()!=initial_hash:
        raise RuntimeError(f'Notebook changed during execution; refusing to overwrite: {path}')
    nbformat.write(nb,path)
    print('Executed and validated',name,flush=True)
