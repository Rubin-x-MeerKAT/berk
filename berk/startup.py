"""

Startup routines and config for Berk

"""

import os, sys

on_rtd=os.environ.get('READTHEDOCS', None)
if on_rtd is not None:
    os.environ['BERK_ROOT']='.'
    os.environ['BERK_MSCACHE']='MSCache'
    os.environ['BERK_PLATFORM']='chpc'

# Settings are hard-coded for now, but could be put into a YAML config file later ---------------------------
config={}

if "BERK_ROOT" not in os.environ.keys():
    raise Exception("You need to set the BERK_ROOT environment variable - this defines the directory where all work will be done.")
if "BERK_MSCACHE" not in os.environ.keys():
    raise Exception("You need to set the BERK_MSCACHE environment variable - this defines the directory where retrieved measurement sets will be stored.")

os.makedirs(os.environ["BERK_MSCACHE"], exist_ok = True)
os.makedirs(os.environ["BERK_ROOT"], exist_ok = True)

# For using 'collect' to fetch data products from (potentially) multiple locations
if 'BERK_NODES_FILE' in os.environ.keys():
    config['nodesFile']=os.environ['BERK_NODES_FILE']
else:
    config['nodesFile']=None

# We can't just make this about the workload manager as it affects what we feed into oxkat
if 'BERK_PLATFORM' not in os.environ.keys():
    raise Exception("Set BERK_PLATFORM environment variable to either 'hippo' or 'chpc'")
if os.environ['BERK_PLATFORM'] == 'chpc':
    config['workloadManager']='pbs'
elif os.environ['BERK_PLATFORM'] == 'hippo':
    config['workloadManager']='slurm'
else:
    raise Exception("Environment variable BERK_PLATFORM is not set to 'hippo' or 'chpc'")

# Processing and data products will be written in sub-dirs here
config['rootDir']=os.environ["BERK_ROOT"]

# This should eventually be the scratch disk
# config['stagingDir']=config['rootDir']+os.path.sep+"staging"
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
    print("Remember to remove %s before running berk if you need to fetch an updated version from the git repository" % (config['oxkatDir']))

# Set-up ----------------------------------------------------------------------------------------------------
if on_rtd is None:
    dirsToMake=[config['processingDir'], config['productsDir'], config['cacheDir']]
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

