import os
import glob
from setuptools import setup
from setuptools import Extension
import versioneer

setup(name='berk',
      version=versioneer.get_version(),
      author='Matt Hilton',
      author_email='matt.hilton@wits.ac.za',
      classifiers=['Development Status :: 2 - Pre-Alpha',
                   'Environment :: Console',
                   'Intended Audience :: Science/Research',
                   'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
                   'Natural Language :: English',
                   'Operating System :: POSIX',
                   'Programming Language :: Python',
                   'Topic :: Scientific/Engineering :: Astronomy'],
      description="A package for producing and managing MeerKAT continuum survey data.",
      packages=['berk'],
      scripts=['bin/berk', 'bin/berk_chain', 'bin/mkat_primary_beam_correct'],
      #install_requires=["astropy >= 3.2",
                        #"numpy >= 1.10",
                        #"matplotlib >= 2.0"]
)
