"""Format-specific readers returning cells, genes and explicit measurement metadata."""
from importlib import import_module


def load(platform, path, annotations, roi=None, **kwargs):
    if platform not in ('atera','visium_hd','stereo_seq'):
        raise ValueError('Native adapters available here: atera, visium_hd, stereo_seq.')
    return import_module(f'{__name__}.{platform}').load(path,annotations,roi,**kwargs)


def write_template(platform, path, output, roi=None, measurement=None):
    """Write blank annotations from measured geometry; never read count matrices."""
    if platform not in ('atera', 'visium_hd', 'stereo_seq'):
        raise ValueError('Native templates available here: atera, visium_hd, stereo_seq.')
    if measurement is not None and not isinstance(measurement, dict):
        raise ValueError('measurement must be an object.')
    return import_module(f'{__name__}.{platform}').template(path, output, roi, measurement=measurement)
