"""

Berk - a package for producing and managing MeerKAT continuum survey data.

"""

#from . import catalogs 
from . import startup

from ._version import get_versions
__version__ = get_versions()['version']
del get_versions

#__all__ = ['startUp']
