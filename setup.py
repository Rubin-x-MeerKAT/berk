import os
import glob
from setuptools import setup
from setuptools import Extension
import versioneer

setup(name='berk',
      version=versioneer.get_version(),
      author='Matt Hilton',
      author_email='matt.hilton@wits.ac.za',
      packages=['berk'],
      scripts=['bin/berk', 'bin/berk_chain', 'bin/mkat_primary_beam_correct'],
)
