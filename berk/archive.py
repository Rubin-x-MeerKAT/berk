"""

This module contains tools for interacting with the MeerKAT archive

"""

import os, sys, glob, subprocess
from .startup import config
from . import jobs

#------------------------------------------------------------------------------------------------------------
def stageMS(captureBlockId):
    """Sets up a measurement set for processing:

        1. Fetch from archive if necessary
        2. Unpack tar.gz if necessary and put measurement set into top-level of staging directory
        3. Return the path to the measurement set itself

    """

    # NOTE: MOVE THIS INTO BERK_UNPACK
    MSPath=config["stagingDir"]+os.path.sep+captureBlockId+"_sdp_l0.ms"
    if os.path.exists(MSPath) == False:
        pathToTGZ=fetchFromArchive(captureBlockId)
        topDir=os.getcwd()
        os.chdir(config['stagingDir'])
        for p in os.walk("scratch"):
            if p[1][0] == os.path.split(MSPath)[-1]:
                break
        # NOTE: We'd need a job to do this
        mvCmd="mv %s/%s ." % (p[0], p[1][0])
        os.chdir(topDir)
    assert(os.path.exists(MSPath))




