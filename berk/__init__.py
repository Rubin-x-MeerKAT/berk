"""

Berk - a package for producing and managing MeerKAT continuum survey data.

"""

#from . import catalogs 
# NOTE: startup is not imported here, so that e.g. crossmatch can be used without setting BERK_* env vars
# startup.initialise() is called when tasks is imported (i.e., by bin/berk)

from . import _version
__version__ = _version.get_versions()['version']

#__all__ = ['startUp']
