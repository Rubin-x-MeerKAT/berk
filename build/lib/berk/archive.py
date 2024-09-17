"""

This module contains tools for interacting with the MeerKAT archive

"""

import os, sys, glob, subprocess
from .startup import config
from . import jobs

#------------------------------------------------------------------------------------------------------------
def fetchFromArchive(captureBlockIdLink):
    """Fetch the dataset from the archive and write into the scratch area.

    Args:
        captureBlockIdLink (str): Link copied from https://archive.sarao.ac.za of the form
            https://archive-gw-1.kat.ac.za/captureBlockId/captureBlockId_sdp_l0.full.rdb?token=longTokenString.

    """

    # NOTE: We don't need routine any more
    captureBlockId=captureBlockIdLink.split("https://archive-gw-1.kat.ac.za/")[-1].split("/")[0]
    msPath=os.environ['BERK_MSCACHE']+os.path.sep+"%s_sdp_l0.ms" % (captureBlockId)
    cmd="mvftoms.py %s --flags cam,data_lost,ingest_rfi -o %s" % (captureBlockIdLink, msPath)
    os.system("screen -S fetch-%s -d -m %s" % (cmd))
    print("Fetching %s" % (msPath))

#------------------------------------------------------------------------------------------------------------
def checkFetchComplete(captureBlockId):
    """Check if the measurement set corresponding to the given captureBlockId has been fetched.

    Note:
        This relies on GNU Screen.

    """

    process=subprocess.run(['screen', '-ls'], universal_newlines = True,
                           stdout = subprocess.PIPE, stderr = subprocess.STDOUT)
    if process.stdout.find("fetch-%s" % (captureBlockId)) == -1:
        return False
    else:
        return True

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
        #jobID=jobs.submitJob("tar -zxvf %s" % (pathToTGZ), "unpack")
        for p in os.walk("scratch"):
            if p[1][0] == os.path.split(MSPath)[-1]:
                break
        # NOTE: We'd need a job to do this
        mvCmd="mv %s/%s ." % (p[0], p[1][0])
        os.chdir(topDir)
    assert(os.path.exists(MSPath))




