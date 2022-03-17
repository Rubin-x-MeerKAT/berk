"""

This module contains tools for interacting with the MeerKAT archive

"""

import os, sys, glob
from .startup import config
from . import jobs

#------------------------------------------------------------------------------------------------------------
def fetchFromArchive(captureBlockId):
    """Dummy routine (for now). Fetch the requested MeerKAT measurement set, identified using
    captureBlockId (the equivalent of e.g. XMM ObsID), as a tar.gz. Return the path to the .tar.gz.

    """

    pathToTGZ=config['stagingDir']+os.path.sep+"%s_sdp_l0.ms.tar.gz" % (captureBlockId)
    print("archive.fetchFromArchive - dummy retrieve - %s" % (pathToTGZ))
    assert(os.path.exists(pathToTGZ))

    return pathToTGZ

#------------------------------------------------------------------------------------------------------------
def stageMS(captureBlockId):
    """Sets up a measurement set for processing:

        1. Fetch from archive if necessary
        2. Unpack tar.gz if necessary and put measurement set into top-level of staging directory
        3. Return the path to the measurement set itself

    """

    MSPath=config["stagingDir"]+os.path.sep+captureBlockId+"_sdp_l0.ms"
    if os.path.exists(MSPath) == False:
        pathToTGZ=fetchFromArchive(captureBlockId)
        topDir=os.getcwd()
        os.chdir(config['stagingDir'])
        print("staging - find/unpack - have to run as a slurm job and check result")
        import IPython
        IPython.embed()
        sys.exit()
        jobID=jobs.submitJob("tar -zxvf %s" % (pathToTGZ), "unpack")
        os.chdir(topDir)
    assert(os.path.exists(MSPath))




