"""Colab setup with explicit interpreter diagnostics and safe restart detection."""
from pathlib import Path
import importlib
import importlib.metadata as metadata
import json
import subprocess
import sys

SUPPORTED_PYTHON = {(3, 12), (3, 13)}
MODULES = {'numpy': 'numpy', 'scipy': 'scipy', 'pandas': 'pandas', 'h5py': 'h5py',
           'pyarrow': 'pyarrow', 'shapely': 'shapely', 'matplotlib': 'matplotlib',
           'anndata': 'anndata', 'numba': 'numba', 'llvmlite': 'llvmlite',
           'zarr': 'zarr', 'numcodecs': 'numcodecs', 'decoupler': 'decoupler',
           'typing_extensions': 'typing_extensions', 'scverse-misc': 'scverse_misc',
           'ipywidgets': 'ipywidgets'}


def check_python(version=None):
    version = tuple(version or sys.version_info[:2])
    if version not in SUPPORTED_PYTHON:
        raise RuntimeError(
            f'This session is running Python {version[0]}.{version[1]}. '
            'This notebook supports Python 3.12 and 3.13. '
            'Use Runtime > Disconnect and delete runtime, then connect to a Python 3 CPU runtime. '
            'The displayed runtime selection does not change an already running interpreter.')
    return version


def changed_loaded_packages(loaded_versions):
    """Installed files cannot replace compiled modules already loaded in this kernel."""
    changed = []
    for package, version in loaded_versions.items():
        try:
            installed = metadata.version(package)
        except metadata.PackageNotFoundError:
            installed = None
        if version and version != installed:
            changed.append(package)
    return sorted(changed)


def setup_environment(code_dir, install_packages=True):
    code_dir = Path(code_dir)
    check_python()
    loaded = {}
    for package, module in MODULES.items():
        if module not in sys.modules:
            continue
        version = getattr(sys.modules[module], '__version__', None)
        if version is None:
            # Some imported modules, notably typing_extensions, lack __version__.
            try:
                version = metadata.version(package)
            except metadata.PackageNotFoundError:
                version = None
        loaded[package] = version
    if install_packages:
        print('Installing the pinned Python 3.12/3.13 dependencies...', flush=True)
        subprocess.check_call([sys.executable, '-m', 'pip', 'install', '--quiet',
                               '--disable-pip-version-check', '--prefer-binary',
                               '-r', str(code_dir / 'requirements-colab.txt')])
    importlib.invalidate_caches()
    changed = changed_loaded_packages(loaded)
    # Prevent an old notebook's imported source from surviving a source update.
    previous = sys.modules.get('xenium_prepare')
    if previous is not None and getattr(previous, 'VERSION', None) != 'lung5k-preparation-1.1.0':
        changed.append('previous preparation code')
    if changed:
        print('\nINSTALLATION COMPLETE — ONE SESSION RESTART IS REQUIRED.')
        print('Select Runtime > Restart session, then run this setup cell again.')
        print('Already loaded modules changed: ' + ', '.join(changed))
        print('Your Drive files are unaffected. Do not continue to the processing cells yet.')
        return {'ready': False, 'restart_required': True, 'changed_packages': changed}
    # Test compiled imports in a fresh process before importing them into this kernel.
    probe = subprocess.run([sys.executable, '-c',
        'import json,sys,numpy,scipy,pandas,h5py,pyarrow,shapely,anndata,numba,decoupler; '
        'print(json.dumps({"python":sys.version.split()[0],"numpy":numpy.__version__,"decoupler":decoupler.__version__}))'],
        capture_output=True, text=True)
    if probe.returncode:
        raise RuntimeError('Dependency verification failed. The complete import error follows:\n' + probe.stderr[-12000:])
    versions = json.loads(probe.stdout.strip().splitlines()[-1])
    print('Dependency verification passed:', versions)
    return {'ready': True, 'restart_required': False, 'versions': versions}
