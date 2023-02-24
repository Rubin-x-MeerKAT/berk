"""

Startup routines and config for hippoxkatapult

"""

import os, sys

# Settings are hard-coded for now, but could be put into a YAML config file later ---------------------------
config={}

if "HIPPOXKATAPULT_ROOT" not in os.environ.keys():
    raise Exception("You need to set the HIPPOXKATAPULT_ROOT environment variable - this defines the directory where all work will be done.")
if "HIPPOXKATAPULT_MSCACHE" not in os.environ.keys():
    raise Exception("You need to set the HIPPOXKATAPULT_MSCACHE environment variable - this defines the directory where retrieved measurement sets will be stored.")

os.makedirs(os.environ["HIPPOXKATAPULT_MSCACHE"], exist_ok = True)
os.makedirs(os.environ["HIPPOXKATAPULT_ROOT"], exist_ok = True)

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
# OLD: Tagged version
# config['oxkatVersion']="0.3"
# config['oxkatDir']=config['cacheDir']+os.path.sep+"oxkat-%s" % (config['oxkatVersion'])
# config['oxkatURL']="https://github.com/IanHeywood/oxkat/archive/refs/tags/v%s.tar.gz" % (config['oxkatVersion'])
# CURRENT: From Matt's git fork
config['oxkatVersion']="git"
config['oxkatDir']=config['cacheDir']+os.path.sep+"oxkat-%s" % (config['oxkatVersion'])
config['oxkatURL']="https://github.com/mattyowl/oxkat.git"

# Image-processing (source finding scripts from Jonah)
config['catalogScriptsDir']=config['cacheDir']+os.path.sep+"catalog-scripts"

print("Using oxkat version: %s" % (config['oxkatVersion']))
if config['oxkatVersion'] == "git":
    print("Remember to remove %s before running hippoxkatapult if you need to fetch an updated version from the git repository" % (config['oxkatDir']))

# Set-up ----------------------------------------------------------------------------------------------------
dirsToMake=[config['stagingDir'], config['processingDir'], config['productsDir'], config['cacheDir']]
for d in dirsToMake:
    os.makedirs(d, exist_ok = True)

if os.path.exists(config['oxkatDir']) == False:
    topDir=os.getcwd()
    os.chdir(config['cacheDir'])
    if config['oxkatVersion'] != 'git':
        os.system("wget %s" % (config['oxkatURL']))
        os.system("tar -zxvf v%s.tar.gz" % (config['oxkatVersion']))
    else:
        os.system("git clone %s oxkat-%s" % (config['oxkatURL'], config['oxkatVersion']))
    os.chdir(topDir)

if os.path.exists(config['catalogScriptsDir']) == False:
    os.system("git clone https://github.com/mattyowl/Image-processing %s" % (config["catalogScriptsDir"]))

