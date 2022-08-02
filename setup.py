# -*- coding: iso-8859-1 -*-
#
# meerkatapult install script

import os
import glob
from setuptools import setup
from setuptools import Extension
import versioneer

setup(name='hippoxkatapult',
      version=versioneer.get_version(),
      author='Matt Hilton',
      author_email='hiltonm@ukzn.ac.za',
      classifiers=['Development Status :: 2 - Pre-Alpha',
                   'Environment :: Console',
                   'Intended Audience :: Science/Research',
                   'License :: OSI Approved :: GNU General Public License v2 (GPLv2)',
                   'Natural Language :: English',
                   'Operating System :: POSIX',
                   'Programming Language :: Python',
                   'Topic :: Scientific/Engineering :: Astronomy'],
      description="Tools for running Oxkat on archival MeerKAT data on Hippo, UKZN's HPC facility.",
      packages=['hippoxkatapult'],
      #package_data={'meerkatapult': ['data/*']},
      scripts=['bin/hippoxkatapult', 'bin/hippoxkatapult_chain'],
      #install_requires=["astropy >= 3.2",
                        #"numpy >= 1.10",
                        #"matplotlib >= 2.0"]
)
