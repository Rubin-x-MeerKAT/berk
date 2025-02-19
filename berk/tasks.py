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
from . import startup, archive, jobs, catalogs, images, crossmatch,  __version__
import shlex
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord, Longitude

#------------------------------------------------------------------------------------------------------------
def fetch(captureBlockId):
    """Fetch...

    Probably need to run ``berk fetch -o rdb_link`` in GNU screen for this task.

    """

    captureBlockIdLink=captureBlockId
    captureBlockId=captureBlockIdLink.split("https://archive-gw-1.kat.ac.za/")[-1].split("/")[0]
    msPath=os.environ['BERK_MSCACHE']+os.path.sep+"%s_sdp_l0.ms" % (captureBlockId)
    fetchLogPath=os.environ['BERK_MSCACHE']+os.path.sep+"%s_fetch.log" % (captureBlockId)

    if os.path.exists(msPath) is True:
        print("Already fetching %s - if it failed, remove %s and try again" % (captureBlockId, msPath))
        sys.exit()

    cmd="mvftoms.py %s --flags cam,data_lost,ingest_rfi -o %s" % (captureBlockIdLink, msPath)
    cmd_subprocess=shlex.split(cmd)
    with subprocess.Popen(cmd_subprocess, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True) as process:
        with open(fetchLogPath, 'w') as logFile:
            for line in process.stdout:
                logFile.write(line)
                print(line, end='')

            for line in process.stderr:
                logFile.write(line)
                print(line, end='')
        process.wait()

    # Automatic setting of KATSDPTELSTATE_ALLOW_PICKLE=1
    katsString="you can allow the pickles to be loaded by setting KATSDPTELSTATE_ALLOW_PICKLE=1 in the environment"
    with open(fetchLogPath, 'r') as file:
        fetchLogMsg = file.read()
        if katsString in fetchLogMsg:
            if os.path.exists(msPath) == True:
                os.rmdir(msPath)
            os.environ["KATSDPTELSTATE_ALLOW_PICKLE"] = "1"
            with subprocess.Popen(cmd_subprocess, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True) as process:
                with open(fetchLogPath, 'w') as logFile:
                    for line in process.stdout:
                        logFile.write(line)
                        print(line, end='')

                    for line in process.stderr:
                        logFile.write(line)
                        print(line, end='')
                process.wait()

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
        bandKey='Others'
        #raise Exception("%f - Not sure what band this is - need to add a band key for it" %freqGHz)

    return bandKey
    
#------------------------------------------------------------------------------------------------------------
def fixRA(table, racol='RA', wrap_angle=180):
    """Returns table with corrected RA wrap.

    """
    fixTable = table.copy()
    fixTable[racol] = Longitude(table[racol], unit=u.deg, wrap_angle=wrap_angle * u.deg).value
    return fixTable
    
#------------------------------------------------------------------------------------------------------------
def listObservations():
    """List observations available on this machine, and check their processing status with the central list.

    """

    if 'BERK_INFO_FILE' in os.environ.keys():
        tab=atpy.Table().read(os.environ['BERK_INFO_FILE'])
    else:
        print("Set BERK_INFO_FILE environment variable to check processing status of observations against central list.")
        tab=None

    # MJH: Don't hard code in things like this, it is broken for me
    globXmatchTab=None
    # globXmatchLink = 'https://dl.dropbox.com/scl/fi/39gma12peyymhadc2fss2/xmatchCat_DECaLS_r_4p0asec.fits?rlkey=v4bazyqfy294ukt9eqmfaeed7&st=7e7nfcde&dl=0'
    # globXmatchTab=atpy.Table().read(globXmatchLink)

    msList=glob.glob(os.environ['BERK_MSCACHE']+os.path.sep+"*_sdp_l0.ms")
    msList.sort()
    
    print("Downloaded observations available locally [by captureBlockId]:")
    for ms in msList:
        captureBlockId=os.path.split(ms)[-1].split("_")[0]
        status="cached"
        globStr=os.environ['BERK_ROOT']+os.path.sep+'processing'+os.path.sep+captureBlockId+os.path.sep+'IMAGES'+os.path.sep+'img_%s_*_pcalmask-MFS-image.fits' % (captureBlockId)
        if len(glob.glob(globStr)) > 0:
            status="processed"
        if tab is not None:
            if captureBlockId in tab['captureBlockId']:
                status="processed_and_analysed"
        if globXmatchTab is not None:
            if captureBlockId in globXmatchTab['captureBlockId']:
                status="processed_analysed_xmatched-DECaLS"
                
        print("   %s    %s" % (captureBlockId, status))

#------------------------------------------------------------------------------------------------------------
def builddb():
    """Build database...

    """
    
    # Build global catalog in each band (L-band, UHF)
    globalTabsDict={'L': None, 'UHF': None, 'Others': None}
    
    # Fixing RA
    tabFilesList=glob.glob(startup.config['productsDir']+os.path.sep+'catalogs'+os.path.sep+'*_bdsfcat.fits')
    for t in tabFilesList:
        if t.find("srl_bdsfcat") == -1:
            tab=atpy.Table().read(t)
            if any(tab['RA'] < 0.0):
            	#print("\nFixing RA for %s" %t)
            	tab = fixRA(tab, racol='RA', wrap_angle=360)
            freqGHz=tab.meta['FREQ0']/1e9
            bandKey=_getBandKey(freqGHz)
            #tab.meta=None # It'd be good to clear this... but the catalog matching stuff wants many things from here
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
            print("\nWrote %s" % (outFileName))
    
    # Quality flags table
    qualFileName=startup.config['productsDir']+os.path.sep+"qualityFlags.csv"
    if os.path.exists(qualFileName) == True:
        qualTab=atpy.Table().read(qualFileName)
    else:
        qualTab=None
    
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
        
        # plotting images and saving in products/images directory
        imageOutDir = startup.config['productsDir']+os.path.sep+"images"
        images.plotImages(imgFile, imageOutDir, ax_label_deg = True, statsDict=statDict)
        
    imgTab=atpy.Table()
    for key in statDict.keys():
        arr=[]
        for s in statsDictList:
            arr.append(s[key])
        imgTab[key]=arr
    # Update quality flag column - we only use quality column from qualTab
    imgTab['quality']=99
    if qualTab is not None:
        for irow in imgTab:
            mask=irow['path'] == qualTab['path']
            if 'quality' in qualTab.keys():
               #pulling already existing quality and assigning 99 to newly added images 
            	if mask.any():
            	    irow['quality']=qualTab[mask]['quality'][0]
            	else:
            	    irow['quality'] = 99
            	    
    # Output
    imgTab.meta['BERKVER']=__version__
    imgTab.meta['DATEMADE']=datetime.date.today().isoformat()
    imgTab.write(outFileName, overwrite = True)
    print("\nWrote %s" % (outFileName))
    imgTab.write(qualFileName, overwrite = True)

    # Generate survey mask in some format - we'll use that to get total survey area    
    
    # cross-matching # to be worked on
    
    """

    opt_survey = 'DECaLS'
    opt_survey_dr = 'DR10'
    optband_to_match = 'r'
    search_radius_arcsec = 4.0
    
    xmatchDirPath = startup.config['productsDir']+os.path.sep+'xmatches'
    os.makedirs(xmatchDirPath, exist_ok = True)

    globalXmatchTabName=startup.config['productsDir']+os.path.sep+"xmatchCat_%s%s_%s_%sasec.fits" %(opt_survey, opt_survey_dr, optband_to_match, str(search_radius_arcsec).replace('.', 'p'))
    globalXmatchTab=None
    radCatFilesList=glob.glob(startup.config['productsDir']+os.path.sep+'catalogs'+os.path.sep+'*srl_bdsfcat.fits')
    
    for radcat in radCatFilesList:
    	catalogName = radcat.split(os.path.sep)[-1]
    	captureBlockId = catalogName.split('_')[3]
    	targetName = (catalogName.split('_1024ch_')[1]).split('_srl_')[0]
    	
    	# making a subscript for this particular match
    	outSubscript = '%s_%s%s_%sband_%sasec' %(catalogName.replace('.fits',''), opt_survey, opt_survey_dr, optband_to_match, str(search_radius_arcsec).replace(".","p"))
    	
    	radcattab = atpy.Table().read(radcat)
    	freqGHz=radcattab.meta['FREQ0']/1e9
    	radbandname=_getBandKey(freqGHz)
    	
    	xmatchTabName = startup.config['productsDir']+os.path.sep+'xmatches'+os.path.sep+"xmatchTable_%s" %outSubscript+".fits"

    	if os.path.exists(xmatchTabName):
    	    xmatchTab = atpy.Table().read(xmatchTabName)
    	else:
    	    xmatchTab = crossmatch.xmatch_berk(radio_cat=radcat, radio_band=radbandname, xmatchTabOutName=xmatchTabName, opt_survey=opt_survey, opt_survey_dr=opt_survey_dr, opt_mag_col = optband_to_match, search_radius_asec = search_radius_arcsec, makePlots=False, radRACol='RA', radDecCol='DEC', eRadRACol='E_RA', eRadDecCol='E_DEC', outSubscript=outSubscript)
    	    
    	    
    	if(xmatchTab):
    	    if(globalXmatchTab is None):
    	        globalXmatchTab = xmatchTab
    	    else:
    	        
    	        globalXmatchTab = atpy.vstack([globalXmatchTab, xmatchTab])
    	        
    globalXmatchTab.write(globalXmatchTabName, overwrite = True)
    print("\nWrote %s" % (globalXmatchTabName))
    
    """
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
        
    # Get rms
    print("Collecting rms images...")
    toPath=startup.config['productsDir']+os.path.sep+"rms"
    os.makedirs(toPath, exist_ok = True)
    for s in stubs:
        rmsPath="processing/*/IMAGES/pbcorr_trim_*_pybdsf/*_rms.fits"
        cmd="rsync -avP %s%s %s" % (s, os.path.sep+rmsPath, toPath)
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

    ####
    # BEGIN NEW
    # Generate the oxkat job scripts then spin through + submit them ourselves
    # Initial setup
    cmd="python3 setups/0_GET_INFO.py %s" % (os.environ['BERK_PLATFORM'])
    os.system(cmd)
    if os.path.exists("submit_info_job.sh") is True:
        blah
    else:
        raise Exception("Failed to generate submit_info_job.sh")
    jobIDs=[]
    jobID=jobs.submit_job("submit_info_job.sh", "SUBMIT_INFO", dependent_job_ids = None,
                          workload_manager = startup.config['workloadManager'],
                          cmd_is_batch_script = True)
    jobIDs.append(jobID)
    # 1GC, FLAG, 2GC can be chained
    cmd="python3 setups/1GC.py %s" % (os.environ['BERK_PLATFORM'])
    jobID=jobs.submit_job(cmd, "1GC", dependent_job_ids = jobIDs,
                          workload_manager = startup.config['workloadManager'])
    jobIDs.append(jobID)
    cmd="python3 setups/FLAG.py %s" % (os.environ['BERK_PLATFORM'])
    jobID=jobs.submit_job(cmd, "SETUP_FLAG_JOBS", dependent_job_ids = jobIDs,
                          workload_manager = startup.config['workloadManager'])
    jobIDs.append(jobID)
    cmd="python3 setups/2GC.py %s" % (os.environ['BERK_PLATFORM'])
    jobID=jobs.submit_job(cmd, "SETUP_2GC_JOBS", dependent_job_ids = jobIDs,
                          workload_manager = startup.config['workloadManager'])
    jobIDs.append(jobID)
    # Chain
    cmd="berk_chain %s submit_1GC_jobs.sh submit_flag_jobs.sh submit_2GC_jobs.sh"\
        % (startup.config['workloadManager'])
    jobID=jobs.submit_job(cmd, "CHAINED_JOBS", dependent_job_ids = jobIDs,
                          workload_manager = startup.config['workloadManager'])
    print("All jobs submitted")
    sys.exit()
    #### END NEW

    # ####
    # # BEGIN OLD
    # # 1GC
    # os.system("python3 setups/1GC.py %s" % os.environ['BERK_PLATFORM'])
    # jobCmds=[]
    # dependent=[]
    # with open("submit_1GC_jobs.sh") as inFile:
    #     lines=inFile.readlines()
    #     for line in lines:
    #         if line.find("sbatch") != -1 and startup.config['workloadManager'] == 'slurm':
    #             sbatchCmd=line[line.find("sbatch") :].split(" |")[0]
    #             if sbatchCmd.find("-d afterok:") != -1:
    #                 sbatchCmd=sbatchCmd.split("}")[-1].strip()
    #                 dependent.append(True)
    #             else:
    #                 sbatchCmd=sbatchCmd.split("sbatch")[-1].strip()
    #                 dependent.append(False)
    #             jobCmds.append(sbatchCmd)
    #         elif line.find("qsub") != -1 and startup.config['workloadManager'] == 'pbs':
    #             qsubCmd=line[line.find("qsub") :].split(" |")[0]
    #             if qsubCmd.find("-W depend=afterok") != -1:
    #                 qsubCmd=qsubCmd.split("}")[-1].strip()
    #                 dependent.append(True)
    #             else:
    #                 qsubCmd=qsubCmd.split("qsub")[-1].strip()
    #                 dependent.append(False)
    #             jobCmds.append(qsubCmd)
    #
    # jobIDs=[]
    # for cmd, dep in zip(jobCmds, dependent):
    #     if dep == False:
    #         dependentJobIDs=None
    #     else:
    #         dependentJobIDs=jobIDs
    #     jobName=os.path.split(cmd)[-1]
    #     jobID=jobs.submit_job(cmd, jobName, dependent_job_ids = dependentJobIDs,
    #                           workload_manager = startup.config['workloadManager'],
    #                           cmd_is_batch_script = True)
    #     jobIDs.append(jobID)
    #
    # # Run the FLAG and 2GC setup scripts as a job, then chain them together
    # cmd="python3 setups/FLAG.py %s" % (os.environ['BERK_PLATFORM'])
    # jobID=jobs.submit_job(cmd, "SETUP_FLAG_JOBS", dependent_job_ids = jobIDs,
    #                       workload_manager = startup.config['workloadManager'])
    # jobIDs.append(jobID)
    # cmd="python3 setups/2GC.py %s" % (os.environ['BERK_PLATFORM'])
    # jobID=jobs.submit_job(cmd, "SETUP_2GC_JOBS", dependent_job_ids = jobIDs,
    #                       workload_manager = startup.config['workloadManager'])
    # jobIDs.append(jobID)
    # cmd="berk_chain %s submit_flag_jobs.sh submit_2GC_jobs.sh" % (startup.config['workloadManager'])
    # jobID=jobs.submit_job(cmd, "CHAIN_FLAG+2GC_JOBS", dependent_job_ids = jobIDs,
    #                       workload_manager = startup.config['workloadManager'])
    # print("All jobs submitted")
    # sys.exit()
    # # END OLD
    # ####

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

        jobID=jobs.submit_job(cmd, 'source-finding-%s' % (imgFileName), dependent_job_ids = None,
                              nodes = 1, tasks = 20, mem = 64000, time = "02:00:00",
                              cmd_is_batch_script = False,
                              workload_manager = startup.config['workloadManager'])
        print("Submitted source finding and analysis job %d" % (jobID))
    sys.exit()

#------------------------------------------------------------------------------------------------------------
def report():
    """Report...

    """

    # We'll tidy this all up later...

    # Get catalogs if we don't have them already...
    firstPath=startup.config['cacheDir']+os.path.sep+"first_14dec17.fits"
    if os.path.exists(firstPath) == False:
        print("Fetching and caching FIRST catalog [this is only done once]")
        topDir=os.getcwd()
        os.chdir(os.path.abspath(startup.config['cacheDir']))
        os.system("wget http://sundog.stsci.edu/first/catalogs/first_14dec17.fits.gz")
        os.system("gunzip first_14dec17.fits.gz")
        os.chdir(topDir)
    firstTab=atpy.Table().read(firstPath)

    # Compare L-band flux densities with FIRST [which only goes to dec -10 degrees]
    tab=atpy.Table().read(os.path.abspath(startup.config['productsDir'])+os.path.sep+"survey_catalog_L.fits")
    x_tab, x_firstTab, rDeg=catalogs.crossMatch(tab, firstTab, radiusArcmin = 6.0/60)

    # /home/matty/Documents/Rubin_x_MeerKAT/products/survey_catalog_L.fits


    # http://sundog.stsci.edu/first/catalogs/first_14dec17.fits.gz

    import IPython
    IPython.embed()
    sys.exit()

