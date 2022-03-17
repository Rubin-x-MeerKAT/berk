"""

Startup routines and config for hippoxkatapult

"""

import os, sys

# Settings are hard-coded for now, but could be put into a YAML config file later ---------------------------
config={}

# Processing and data products will be written in sub-dirs here
config['rootDir']="/data/mjh/MeerKAT/MSS"

# This should eventually be the scratch disk
config['stagingDir']=config['rootDir']+os.path.sep+"staging"
# stagingDir="/scratch/month/mjh/"

# Directory where we'll cache e.g. oxkat
config['cacheDir']="/data/mjh/hippoxkatapult-cache"

# Oxkat version to use
config['oxkatVersion']="0.3"
config['oxkatDir']=config['cacheDir']+os.path.sep+"oxkat-%s" % (config['oxkatVersion'])

# Set-up ----------------------------------------------------------------------------------------------------
dirsToMake=[config['stagingDir'], config['cacheDir'], config['oxkatDir']]
for d in dirsToMake:
    os.makedirs(d, exist_ok = True)

if os.path.exists(config['oxkatDir']) == False:
    topDir=os.getcwd()
    os.chdir(config['cacheDir'])
    os.system("wget %s" % (config['oxkatURL']))
    os.system("tar -zxvf v%.s.tar.gz" % (config['oxkatVersion']))
    os.chdir(topDir)
