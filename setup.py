import os
import glob
from setuptools import setup
from setuptools import Extension
import versioneer

setup(name='berk',
      version=versioneer.get_version(),
      cmdclass=versioneer.get_cmdclass(),
      author='Matt Hilton, Unnikrishnan Sureshkumar',
      author_email='matt.hilton@wits.ac.za',
      packages=['berk'],
      scripts=['bin/berk', 'bin/berk_chain', 'bin/mkat_primary_beam_correct', 'bin/xmatch'],
)
