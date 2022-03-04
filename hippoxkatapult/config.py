"""

Settings for hippoxkatapult

These are hard-coded for now, but could be put into a config file later

"""

import os

# Point this at the version of Oxkat to fetch and use
oxkatURL="https://github.com/IanHeywood/oxkat/archive/refs/tags/v0.3.tar.gz"

# Processing and data products will be written in sub-dirs here
rootDir="/data/mjh/MeerKAT/MSS/"

# This should eventually be the scratch disk
stagingDir=rootDir+os.path.sep+"staging"
#stagingDir="/scratch/month/mjh/"

# Directory where we'll cache e.g. oxkat
cacheDir="/data/mjh/hippoxkatapult-cache/"
