"""

Tasks that can be performed by Berk

"""

import os
import sys
import subprocess
import glob
import time
import datetime
import astropy.table as atpy
from . import startup, archive, jobs, catalogs, images,  __version__

#------------------------------------------------------------------------------------------------------------
def fetch(captureBlockId):
    """Fetch...

    """

    captureBlockIdLink=captureBlockId
    captureBlockId=captureBlockIdLink.split("https://archive-gw-1.kat.ac.za/")[-1].split("/")[0]
    msPath=os.environ['BERK_MSCACHE']+os.path.sep+"%s_sdp_l0.ms" % (captureBlockId)
    if archive.checkFetchComplete(captureBlockId) == False:
        cmd="mvftoms.py %s --flags cam,data_lost,ingest_rfi -o %s" % (captureBlockIdLink, msPath)
        os.system(cmd)
        # print("Run this command in GNU screen:")
        # print(cmd)
        # print("[yes, this is clunky - but automatically running in screen isn't working at the moment]")
        # os.system("screen -dmS fetch-%s bash -c '%s'" % (captureBlockId, cmd))
        # print("Fetching %s" % (msPath))
    else:
        print("Already fetching %s" % (msPath))
    sys.exit()

#------------------------------------------------------------------------------------------------------------
def store(captureBlockId):
    """Store...

    """

    print("Task 'store' is not implemented yet.")
    sys.exit()

#------------------------------------------------------------------------------------------------------------
def _getBandKey(freqGHz):
    """Return band name ('L' or 'UHF') based on freqGHz.

    """
    if freqGHz > 1.2 and freqGHz < 1.3:
        bandKey='L'
    elif freqGHz > 0.7 and freqGHz < 0.9:
        bandKey='UHF'
    else:
        raise Exception("Not sure what band this is - need to add a band key for it")

    return bandKey

#------------------------------------------------------------------------------------------------------------
def listObservations():
    """List observations available on this machine, and check their processing status with the central list.

    """

    if 'BERK_INFO_FILE' in os.environ.keys():
        tab=atpy.Table().read(os.environ['BERK_INFO_FILE'])
    else:
        print("Set BERK_INFO_FILE environment variable to check processing status of observations against central list.")
        tab=None

    msList=glob.glob(os.environ['BERK_MSCACHE']+os.path.sep+"*_sdp_l0.ms")
    msList.sort()
    print("Downloaded observations available locally [by captureBlockId]:")
    for ms in msList:
        captureBlockId=os.path.split(ms)[-1].split("_")[0]
        status="cached"
        if tab is not None:
            if captureBlockId in tab['captureBlockId']:
                status="processed_and_analysed"
        print("   %s    %s" % (captureBlockId, status))

#------------------------------------------------------------------------------------------------------------
def builddb():
    """Build database...

    """

    # Build global catalog in each band (L-band, UHF)
    globalTabsDict={'L': None, 'UHF': None}
    tabFilesList=glob.glob(startup.config['productsDir']+os.path.sep+'catalogs'+os.path.sep+'*_bdsfcat.fits')
    for t in tabFilesList:
        if t.find("srl_bdsfcat") == -1:
            tab=atpy.Table().read(t)
            freqGHz=tab.meta['FREQ0']/1e9
            bandKey=_getBandKey(freqGHz)
            if globalTabsDict[bandKey] is None:
                globalTabsDict[bandKey]=tab
            else:
                globalTabsDict[bandKey]=atpy.vstack([globalTabsDict[bandKey], tab])

    # To implement: remove duplicates [we may want to make time series?]

    # Sorting, additional meta data, output
    for bandKey in globalTabsDict.keys():
        if globalTabsDict[bandKey] is not None:
            outFileName=startup.config['productsDir']+os.path.sep+"survey_catalog_%s.fits" % (bandKey)
            globalTabsDict[bandKey].sort('DEC')
            globalTabsDict[bandKey].sort('RA')
            globalTabsDict[bandKey].meta['BAND']=bandKey
            globalTabsDict[bandKey].meta['BERKVER']=__version__
            globalTabsDict[bandKey].meta['DATEMADE']=datetime.date.today().isoformat()
            globalTabsDict[bandKey].write(outFileName, overwrite = True)
            catalogs.catalog2DS9(globalTabsDict[bandKey], outFileName.replace(".fits", ".reg"),
                                 idKeyToUse = 'Source_name', RAKeyToUse = 'RA', decKeyToUse = 'DEC')
            print("Wrote %s" % (outFileName))

    # Make image table - centre coords, radius [approx.], RMS, band, image path - UHF and L together.
    # Report command (when we make it) could load and dump some of that info
    outFileName=startup.config['productsDir']+os.path.sep+"images.fits"
    imgFilesList=glob.glob(startup.config['productsDir']+os.path.sep+"images"+os.path.sep+"pbcorr_*.fits")
    statsDictList=[]
    for imgFile in imgFilesList:
        statDict=images.getImagesStats(imgFile)
        captureBlockId=os.path.split(statDict['path'])[-1].split('img_')[-1].split('_sdp')[0]
        statDict['captureBlockId']=captureBlockId
        statDict['path']=statDict['path'].replace(startup.config['productsDir']+os.path.sep, '')
        statDict['band']=_getBandKey(statDict['freqGHz'])
        statsDictList.append(statDict)
    imgTab=atpy.Table()
    for key in statDict.keys():
        arr=[]
        for s in statsDictList:
            arr.append(s[key])
        imgTab[key]=arr
    imgTab.meta['BERKVER']=__version__
    imgTab.meta['DATEMADE']=datetime.date.today().isoformat()
    imgTab.write(outFileName, overwrite = True)
    print("Wrote %s" % (outFileName))

    # Generate survey mask in some format - we'll use that to get total survey area


#------------------------------------------------------------------------------------------------------------
def collect():
    """Collect...

    """

    print("Collecting processed data products...")
    if 'BERK_NODES_FILE' not in os.environ.keys():
        print("You need to set the BERK_NODES_FILE environment variable to use the 'collect' task.")
        sys.exit()
    nodesFilePath=os.environ['BERK_NODES_FILE']
    try:
        stubs=[]
        with open(os.environ['BERK_NODES_FILE'], "r") as inFile:
            for line in inFile.readlines():
                stubs.append(line)
    except:
        import urllib.request  # the lib that handles the url stuff
        stubs=[]
        for line in urllib.request.urlopen(os.environ['BERK_NODES_FILE']):
            l=line.decode('utf-8')
            if l[0] != "#" and len(l) > 3:
                stubs.append(l.strip())

    # Get images
    print("Collecting images...")
    toPath=startup.config['productsDir']+os.path.sep+"images"
    os.makedirs(toPath, exist_ok = True)
    for s in stubs:
        imgPath="processing/*/IMAGES/pbcorr*pcalmask-MFS-image.fits"
        cmd="rsync -avP %s%s %s" % (s, os.path.sep+imgPath, toPath)
        os.system(cmd)

    # Get catalogs
    print("Collecting catalogs...")
    toPath=startup.config['productsDir']+os.path.sep+"catalogs"
    os.makedirs(toPath, exist_ok = True)
    for s in stubs:
        catPath="processing/*/IMAGES/pbcorr_trim_*_pybdsf/*_bdsfcat.fits"
        cmd="rsync -avP %s%s %s" % (s, os.path.sep+catPath, toPath)
        os.system(cmd)

    print("Finished!")
    sys.exit()

#------------------------------------------------------------------------------------------------------------
def process(captureBlockId):
    """Process...

    """

    # Forget staging dir, just do a symbolic link to the MSCache dir
    MSPath=os.environ['BERK_MSCACHE']+os.path.sep+captureBlockId+"_sdp_l0.ms"

    # Setup in processing dir
    MSProcessDir=startup.config['processingDir']+os.path.sep+captureBlockId
    if os.path.exists(MSProcessDir) == True:
        raise Exception("Processing directory %s exists and is not empty - remove it and re-run, if you're sure you don't need its contents." % (MSProcessDir))
    os.makedirs(MSProcessDir)
    os.chdir(MSProcessDir)
    os.system("ln -s %s" % (os.path.abspath(MSPath)))
    oxdirs=['setups', 'tools', 'oxkat', 'data']
    for oxdir in oxdirs:
        os.system("ln -s %s" % (startup.config['oxkatDir']+os.path.sep+oxdir))

    # Generate the oxkat job scripts then spin through + submit them ourselves
    # 1GC
    os.system("python3 setups/1GC.py %s" % os.environ['BERK_PLATFORM'])
    jobCmds=[]
    dependent=[]
    with open("submit_1GC_jobs.sh") as inFile:
        lines=inFile.readlines()
        for line in lines:
            if line.find("sbatch") != -1 and startup.config['workloadManager'] == 'slurm':
                sbatchCmd=line[line.find("sbatch") :].split(" |")[0]
                if sbatchCmd.find("-d afterok:") != -1:
                    sbatchCmd=sbatchCmd.split("}")[-1].strip()
                    dependent.append(True)
                else:
                    sbatchCmd=sbatchCmd.split("sbatch")[-1].strip()
                    dependent.append(False)
                jobCmds.append(sbatchCmd)
            elif line.find("qsub") != -1 and startup.config['workloadManager'] == 'pbs':
                qsubCmd=line[line.find("qsub") :].split(" |")[0]
                if qsubCmd.find("-W depend=afterok") != -1:
                    qsubCmd=qsubCmd.split("}")[-1].strip()
                    dependent.append(True)
                else:
                    qsubCmd=qsubCmd.split("qsub")[-1].strip()
                    dependent.append(False)
                jobCmds.append(qsubCmd)

    jobIDs=[]
    for cmd, dep in zip(jobCmds, dependent):
        if dep == False:
            dependentJobIDs=None
        else:
            dependentJobIDs=jobIDs
        jobName=os.path.split(cmd)[-1]
        jobID=jobs.submitJob(cmd, jobName, dependentJobIDs = dependentJobIDs, workloadManager = startup.config['workloadManager'], cmdIsBatchScript = True)
        jobIDs.append(jobID)

    # Run the FLAG and 2GC setup scripts as a job, then chain them together
    cmd="python3 setups/FLAG.py %s" % (os.environ['BERK_PLATFORM'])
    jobID=jobs.submitJob(cmd, "SETUP_FLAG_JOBS", dependentJobIDs = jobIDs, workloadManager = startup.config['workloadManager'])
    jobIDs.append(jobID)
    cmd="python3 setups/2GC.py %s" % (os.environ['BERK_PLATFORM'])
    jobID=jobs.submitJob(cmd, "SETUP_2GC_JOBS", dependentJobIDs = jobIDs, workloadManager = startup.config['workloadManager'])
    jobIDs.append(jobID)
    cmd="berk_chain %s submit_flag_jobs.sh submit_2GC_jobs.sh" % (startup.config['workloadManager'])
    jobID=jobs.submitJob(cmd, "CHAIN_FLAG+2GC_JOBS", dependentJobIDs = jobIDs, workloadManager = startup.config['workloadManager'])
    print("All jobs submitted")
    sys.exit()

#------------------------------------------------------------------------------------------------------------
def analyse(captureBlockId):
    """Analyse...

    """

    # Setup in processing dir
    MSProcessDir=startup.config['processingDir']+os.path.sep+captureBlockId
    if os.path.exists(MSProcessDir) == False:
        raise Exception("Processing directory %s does not exist - you need to process the data before the analyse task will run." % (MSProcessDir))
    os.chdir(MSProcessDir)
    os.system("ln -s %s" % (startup.config["catalogScriptsDir"]+os.path.sep+"sourcefinding.py"))
    os.system("ln -s %s" % (startup.config["catalogScriptsDir"]+os.path.sep+"catalog_matching.py"))
    os.system("ln -s %s" % (startup.config["catalogScriptsDir"]+os.path.sep+"parsets"))

    # Source finding is fairly lightweight so we put everything in one job script
    # We will have issues with needing to see the internet to fetch cross match catalogs though
    # So we will need to cache NVSS catalogs for a given direction when doing 'stage'
    # OR forget using catalog_matching.py here and just do later with database/catalog scripts
    imgPaths=glob.glob("IMAGES/*pcalmask-MFS-image.fits")
    for i in imgPaths:
        if i.find("pbcorr_trim") != -1:
            continue
        imgPath=os.path.abspath(i)
        cmd="mkat_primary_beam_correct %s -T" % (imgPath)

        imgDir, imgFileName=os.path.split(imgPath)
        pbcorrImgPath=imgDir+os.path.sep+"pbcorr_trim_"+imgFileName
        cmd=cmd+"\npython3 sourcefinding.py c %s -o fits --survey MSS" % (pbcorrImgPath)

        # The bit below isn't going to work on compute nodes - so we may as well move this
        label=imgFileName.split(".ms_")[0].split(".")[0]
        catPath=imgDir+os.path.sep+"pbcorr_trim_"+label+"_pybdsf"+os.path.sep+"pbcorr_trim_"+label+"_bdsfcat.fits"
        cmd=cmd+"\npython3 catalog_matching.py %s NVSS --astro --flux" % (catPath)

        jobID=jobs.submitJob(cmd, 'source-finding-%s' % (imgFileName), dependentJobIDs = None, nodes = 1, tasks = 20, mem = 64000,
                                time = "02:00:00", cmdIsBatchScript = False, workloadManager = startup.config['workloadManager'])
        print("Submitted source finding and analysis job %d" % (jobID))
    sys.exit()



