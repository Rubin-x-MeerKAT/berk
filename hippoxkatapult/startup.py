"""

Startup routines and config for hippoxkatapult

"""

import os, sys

# Settings are hard-coded for now, but could be put into a YAML config file later ---------------------------
config={}

# Processing and data products will be written in sub-dirs here
config['rootDir']=os.environ["HIPPOXKATAPULT_ROOT"]

# This should eventually be the scratch disk
config['stagingDir']=config['rootDir']+os.path.sep+"staging"
# stagingDir="/scratch/month/mjh/"

# Directory where we do processing
config['processingDir']=config['rootDir']+os.path.sep+"processing"

# Directory where we archive data products
config['productsDir']=config['rootDir']+os.path.sep+"products"

# Directory where we'll cache e.g. oxkat
config['cacheDir']=config['rootDir']+os.path.sep+"cache"

# Oxkat version to use
config['oxkatVersion']="0.3"
config['oxkatDir']=config['cacheDir']+os.path.sep+"oxkat-%s" % (config['oxkatVersion'])
config['oxkatURL']="https://github.com/IanHeywood/oxkat/archive/refs/tags/v%s.tar.gz" % (config['oxkatVersion'])

# Image-processing (source finding scripts from Jonah)
config['catalogScriptsDir']=config['cacheDir']+os.path.sep+"catalog-scripts"

# Set-up ----------------------------------------------------------------------------------------------------
dirsToMake=[config['stagingDir'], config['processingDir'], config['productsDir'], config['cacheDir']]
for d in dirsToMake:
    os.makedirs(d, exist_ok = True)

if os.path.exists(config['oxkatDir']) == False:
    topDir=os.getcwd()
    os.chdir(config['cacheDir'])
    os.system("wget %s" % (config['oxkatURL']))
    os.system("tar -zxvf v%s.tar.gz" % (config['oxkatVersion']))
    os.chdir(topDir)

if os.path.exists(config['catalogScriptsDir']) == False:
    os.system("git clone https://github.com/mattyowl/Image-processing %s" % (config["catalogScriptsDir"]))

