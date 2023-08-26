"""

This module contains tools for handling catalogs, which are usually :obj:`astropy.table.Table` objects.

"""

from astLib import *
import numpy as np
import operator
import os
import sys
import time
import astropy.table as atpy
from astropy.coordinates import SkyCoord
from astropy.coordinates import match_coordinates_sky
import astropy.io.fits as pyfits
from scipy import ndimage
from . import __version__

# For adding meta data to output
import datetime

#------------------------------------------------------------------------------------------------------------
def catalog2DS9(catalog, outFileName, constraintsList = [], addInfo = [], idKeyToUse = 'name',\
                RAKeyToUse = 'RADeg', decKeyToUse = 'decDeg', color = "cyan", showNames = True,\
                writeBerkInfo = True, coordSys = 'fk5', regionShape = 'point', width = 1):
    """Writes a DS9 region file corresponding to the given catalog. 
    
    Args:
        catalog (:obj:`astropy.table.Table`): An astropy Table where each row represents an object.
        outFileName (:obj:`str`): A file name for the output DS9 region file.
        constraintsList (:obj:`list`, optional): A list of constraints in the same format as used by 
            :func:`selectFromCatalog`.
        addInfo (:obj:`list`, optional): A list of dictionaries with keys named `key` and `fmt` (e.g., 
            ``{'key': "SNR", 'fmt': "%.3f"}``). These will be added to the object label shown in DS9.
        idKeyToUse (:obj:`str`, optional): The name of the key in each object dictionary that defines the 
            object's name. Used to label objects in the DS9 region file.
        RAKeyToUse (:obj:`str`, optional): The name of the key in each object dictionary that contains the 
            RA of the object in decimal degrees.
        decKeyToUse (:obj:`str`, optional): The name of the key in each object dictionary that contains the
            declination of the object in decimal degrees.
        color (:obj:`str`, optional): The color of the plot symbol used by DS9.
        writeBerkInfo (:obj:`bool`, optional): If ``True``, writes a line with the `berk` version and date
            generated at the top of the DS9 .reg file.
        coordSys (:obj:`str`, optional): A string defining the coordinate system used for RA, dec, as 
            understood by DS9.
        
    Returns:
        None
    
    """
    
    cutCatalog=selectFromCatalog(catalog, constraintsList) 
    
    with open(outFileName, "w") as outFile:
        timeStamp=datetime.datetime.today().date().isoformat()
        comment="# DS9 region file"
        if writeBerkInfo == True:
            comment=comment+" generated by Berk (version: %s on %s)\n" % (__version__, timeStamp)
        else:
            comment=comment+"\n"
        outFile.write(comment)
        outFile.write('global dashlist=8 3 width=%d font="helvetica 10 normal" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n' % (width))
        for obj in cutCatalog:
            if len(addInfo) > 0:
                infoString=""
                for d in addInfo:
                    if infoString != "":
                        infoString=infoString+" "
                    if obj[d['key']] != None:
                        infoString=infoString+d['fmt'] % (obj[d['key']])
                    else:
                        infoString=infoString+"%s" % (str(obj[d['key']]))
                infoString=" ["+infoString+"]"
            else:
                infoString=""
            if color == 'key':
                colorString=obj['color']
            else:
                colorString=color
            if showNames == True:
                infoString=str(obj[idKeyToUse])+infoString
            if regionShape == 'point':
                outFile.write("%s;point(%.6f,%.6f) # point=cross color={%s} text={%s}\n" \
                            % (coordSys, obj[RAKeyToUse], obj[decKeyToUse], colorString, infoString))
            elif regionShape == 'circle':
                outFile.write('%s;circle(%.6f,%.6f,360") # color={%s} text={%s}\n' \
                            % (coordSys, obj[RAKeyToUse], obj[decKeyToUse], colorString, infoString))                

#------------------------------------------------------------------------------------------------------------
def makeName(RADeg, decDeg, prefix = 'MKT'):
    """Makes an object name string from the given object coordinates, following the IAU convention.
    
    Args:
        RADeg (:obj:`float`): Right ascension of the object in J2000 decimal degrees.
        decDeg (:obj:`float`): Declination of the object in J2000 decimal degrees.
        prefix (:obj:`str`, optional): Prefix for the object name.
    
    Returns:
        Object name string in the format `prefix JHHMM.m+/-DDMM`.
    
    """
    
    actName=prefix+" J"+_makeRA(RADeg)+_makeDec(decDeg)
    
    return actName

#------------------------------------------------------------------------------------------------------------
def makeLongName(RADeg, decDeg, prefix = "MKT"):
    """Makes a long format object name string from the given object coordinates, following the IAU convention.
    
    Args:
        RADeg (:obj:`float`): Right ascension of the object in J2000 decimal degrees.
        decDeg (:obj:`float`): Declination of the object in J2000 decimal degrees.
        prefix (:obj:`str`, optional): Prefix for the object name.
    
    Returns:
        Object name string in the format `prefix JHHMMSS.s+/-DDMMSS`.
    
    """
    
    actName=prefix+" J"+_makeLongRA(RADeg)+_makeLongDec(decDeg)
    
    return actName
    
#------------------------------------------------------------------------------------------------------------
def _makeRA(myRADeg):
    """Makes RA part of ACT names.
    
    """
    hours=(myRADeg/360)*24
    strHours=("%.10f" % (hours))
    if hours<10:
        sHours="0"+strHours[0]
    else:
        sHours=strHours[:2]
    
    mins=float(strHours[strHours.index("."):])*60
    strMins=("%.10f" % (mins))
    if mins < 10:
        sMins="0"+strMins[:3]
    else:
        sMins=strMins[:4]

    return (sHours+sMins)#[:-2] # Trims off .x as not used in ACT names
        
#------------------------------------------------------------------------------------------------------------
def _makeDec(myDecDeg):
    """Makes dec part of ACT names
    
    """
    
    # Positive
    if myDecDeg>0:
        if myDecDeg<10:
            sDeg="0"+str(myDecDeg)[0]
        else:
            sDeg=str(myDecDeg)[:2]
    
        mins=float(str(myDecDeg)[str(myDecDeg).index("."):])*60
        if mins<10:
            sMins="0"+str(mins)[:1]
        else:
            sMins=str(mins)[:2]
        
        return "+"+sDeg+sMins
    else:
        if myDecDeg>-10:
            sDeg="-0"+str(myDecDeg)[1]
        else:
            sDeg=str(myDecDeg)[:3]
    
        mins=float(str(myDecDeg)[str(myDecDeg).index("."):])*60
        if mins<10:
            sMins="0"+str(mins)[:1]
        else:
            sMins=str(mins)[:2]
        
        return str(sDeg+sMins)

#-------------------------------------------------------------------------------------------------------------
def _makeLongRA(myRADeg):
    """Make a long RA string, i.e. in style of long XCS names
    
    """
    
    hours=(myRADeg/360)*24
    if hours<10:
        sHours="0"+str(hours)[0]
    else:
        sHours=str(hours)[:2]
    
    mins=float(str(hours)[str(hours).index("."):])*60
    if mins<10:
        sMins="0"+str(mins)[0]
    else:
        sMins=str(mins)[:2]
        
    secs=float(str(mins)[str(mins).index("."):])*60
    if secs<10:
        sSecs="0"+str(secs)[:3]
    else:
        sSecs=str(secs)[:4]     
        
    return sHours+sMins+sSecs
        
#-------------------------------------------------------------------------------------------------------------
def _makeLongDec(myDecDeg):
    """Make a long dec sting i.e. in style of long XCS names
    
    """
    # Positive
    if myDecDeg>0:
        if myDecDeg<10:
            sDeg="0"+str(myDecDeg)[0]
        else:
            sDeg=str(myDecDeg)[:2]
    
        mins=float(str(myDecDeg)[str(myDecDeg).index("."):])*60
        if mins<10:
            sMins="0"+str(mins)[:1]
        else:
            sMins=str(mins)[:2]
            
        secs=float(str(mins)[str(mins).index("."):])*60
        if secs<10:
            sSecs="0"+str(secs)[:3]
        else:
            sSecs=str(secs)[:4]         
        
        return "+"+sDeg+sMins+sSecs
    else:
        if myDecDeg>-10:
            sDeg="-0"+str(myDecDeg)[1]
        else:
            sDeg=str(myDecDeg)[:3]
    
        mins=float(str(myDecDeg)[str(myDecDeg).index("."):])*60
        if mins<10:
            sMins="0"+str(mins)[:1]
        else:
            sMins=str(mins)[:2]
            
        secs=float(str(mins)[str(mins).index("."):])*60
        if secs<10:
            sSecs="0"+str(secs)[:3]
        else:
            sSecs=str(secs)[:4]         
        
        return sDeg+sMins+sSecs
        
#-------------------------------------------------------------------------------------------------------------
def selectFromCatalog(catalog, constraintsList):
    """Return a table of objects matching the given constraints from the catalog. 
    
    Args:
        catalog (:obj:`astropy.table.Table`): The catalog from which objects will be selected.
        constraintsList (:obj:`list`): A list of constraints, where each item is a string of the form
            "key < value", "key > value", etc.. Note that the spaces between the key, operator 
            (e.g. '<'), and value are essential.
    
    Returns:
        An astropy Table object.
    
    """
            
    passedConstraint=catalog
    for constraintString in constraintsList:
        key, op, value=constraintString.split()
        passedConstraint=passedConstraint[eval("passedConstraint['%s'] %s %s" % (key, op, value))]

    return passedConstraint

#------------------------------------------------------------------------------------------------------------
def removeDuplicates(tab):
    """Removes duplicate objects from the catalog - keeping the highest SNR detection for each duplicate. 
    This routine is used to clean up the output of MPI runs (where we have overlapping tiles).
    
    Args:
        tab (:obj:`astropy.table.Table`): The object catalog to be checked for duplicates.

    Returns:
        Table with duplicates removed (:obj:`astropy.table.Table`), the number of duplicates found, and a
        list of names for the duplicated objects.
    
    """

    if len(tab) == 1:
        return tab, 1, []
    
    # Find all duplicates
    cat=SkyCoord(ra = tab['RADeg'].data, dec = tab['decDeg'].data, unit = 'deg')
    xIndices, rDeg, sep3d = match_coordinates_sky(cat, cat, nthneighbor = 2)
    mask=np.less(rDeg.value, XMATCH_RADIUS_DEG)
    noDupMask=np.greater_equal(rDeg.value, XMATCH_RADIUS_DEG)
    dupTab=tab[mask]
    noDupTab=tab[noDupMask]
    
    # All duplicates removed?
    if mask.sum() == 0:
        return tab, 0, []
    
    # Much faster
    keepMask=np.zeros(len(dupTab), dtype = bool)
    for i in range(len(dupTab)):
        # NOTE: astCoords does not like atpy.Columns sometimes...
        rDeg=astCoords.calcAngSepDeg(dupTab['RADeg'][i], dupTab['decDeg'][i], dupTab['RADeg'].data, dupTab['decDeg'].data)
        mask=np.less_equal(rDeg, XMATCH_RADIUS_DEG)
        if mask.sum() == 0:	# This ought not to be possible but catch anyway
            bestIndex=i
        else:
            indices=np.where(mask == True)[0]
            bestIndex=indices[np.equal(dupTab['SNR'][mask], dupTab['SNR'][mask].max())][0]
        keepMask[bestIndex]=True
    keepTab=dupTab[keepMask]
    
    keepTab=atpy.vstack([keepTab, noDupTab])
    keepTab.sort('RADeg')
    
    return keepTab, len(dupTab), dupTab['name']

#------------------------------------------------------------------------------------------------------------
def crossMatch(refCatalog, matchCatalog, radiusArcmin = 2.5):
    """Cross matches `matchCatalog` onto `refCatalog` for objects found within some angular radius 
    (specified in arcmin).
    
    Args:
        refCatalog (:obj:`astropy.table.Table`): The reference catalog.
        matchCatalog (:obj:`astropy.table.Table`): The catalog to match onto the reference catalog.
        radiusArcmin (:obj:`float`, optional): Cross-match radius in arcmin.
    
    Returns:
        Cross-matched reference catalog, matchCatalog, and array of angular separation in degrees, for 
        objects in common within the matching radius. The cross matched columns are sorted such that rows in
        each correspond to the matched objects.
    
    """
    
    inTab=refCatalog
    outTab=matchCatalog
    RAKey1, decKey1=getTableRADecKeys(inTab)
    RAKey2, decKey2=getTableRADecKeys(outTab)
    cat1=SkyCoord(ra = inTab[RAKey1].data, dec = inTab[decKey1].data, unit = 'deg')
    xMatchRadiusDeg=radiusArcmin/60.
    cat2=SkyCoord(ra = outTab[RAKey2].data, dec = outTab[decKey2].data, unit = 'deg')
    xIndices, rDeg, sep3d = match_coordinates_sky(cat1, cat2, nthneighbor = 1)
    mask=np.less(rDeg.value, xMatchRadiusDeg)  
    matched_outTab=outTab[xIndices]
    inTab=inTab[mask]
    matched_outTab=matched_outTab[mask]
    rDeg=rDeg.value[mask]
    
    return inTab, matched_outTab, rDeg

#------------------------------------------------------------------------------------------------------------
def removeCrossMatched(refCatalog, matchCatalog, radiusArcmin = 2.5):
    """Cross matches `matchCatalog` onto `refCatalog` for objects found within some angular radius 
    (specified in arcmin), and returns `refCatalog` with the matching entries removed.
    
    Args:
        refCatalog (:obj:`astropy.table.Table`): The reference catalog.
        matchCatalog (:obj:`astropy.table.Table`): The catalog to match onto the reference catalog.
        radiusArcmin (:obj:`float`, optional): Cross-match radius in arcmin.
    
    Returns:
        Cross-matched reference catalog (:obj:`astropy.table.Table`) with matches to `matchCatalog` removed.
        
    """
        
    inTab=refCatalog
    outTab=matchCatalog
    RAKey1, decKey1=getTableRADecKeys(inTab)
    RAKey2, decKey2=getTableRADecKeys(outTab)
    cat1=SkyCoord(ra = inTab[RAKey1].data, dec = inTab[decKey1].data, unit = 'deg')
    xMatchRadiusDeg=radiusArcmin/60.
    cat2=SkyCoord(ra = outTab[RAKey2].data, dec = outTab[decKey2].data, unit = 'deg')
    xIndices, rDeg, sep3d = match_coordinates_sky(cat1, cat2, nthneighbor = 1)
    mask=np.greater(rDeg.value, xMatchRadiusDeg)  
    inTab=inTab[mask]
    
    return inTab
    
#------------------------------------------------------------------------------------------------------------
def getTableRADecKeys(tab):
    """Returns the column names in the table in which RA, dec coords are stored, after trying a few possible
    name variations.
    
    Args:
        tab (:obj:`astropy.table.Table`): The table to search.
        
    Returns:
        Name of the RA column, name of the dec. column
    
    """
    RAKeysToTry=['ra', 'RA', 'RADeg']
    decKeysToTry=['dec', 'DEC', 'decDeg', 'Dec']
    RAKey, decKey=None, None
    for key in RAKeysToTry:
        if key in tab.keys():
            RAKey=key
            break
    for key in decKeysToTry:
        if key in tab.keys():
            decKey=key
            break
    if RAKey is None or decKey is None:
        raise Exception("Couldn't identify RA, dec columns in the supplied table.")
    
    return RAKey, decKey

#------------------------------------------------------------------------------------------------------------
def getCatalogWithinImage(tab, shape, wcs, mask = None):
    """Returns the subset of the catalog with coordinates within the image defined by the given `shape`, 
    `wcs`. Optionally, a `mask` may also be applied.
    
    Args:
        tab (:obj:`astropy.table.Table`): Catalog, as an astropy Table object. Must have columns called 
            'RADeg', 'decDeg' that contain object coordinates in decimal degrees.
        shape (:obj:`list`): Shape of the array corresponding to the image / map.
        wcs (:obj:`astWCS.WCS`): WCS of the image.
        mask (optional, :obj:`np.ndarray`): Mask with same dimensions and WCS as the image. Pixels with
            value = 1 indicate valid area, and pixels with value = 0 are considered to be outside the mask.
            If this is given, the returned catalog will contain only objects in the valid area defined by
            this image mask.
    
    Returns:
        An astropy Table containing the subset of objects within the image.
    
    """
    
    xyCoords=np.array(wcs.wcs2pix(tab['RADeg'].tolist(), tab['decDeg'].tolist()))
    selected=[]
    for i in range(len(tab)):
        x, y=xyCoords[i][0], xyCoords[i][1]
        if np.isnan(x) == True or np.isnan(y) == True:
            selected.append(False)
            continue
        if x >= 0 and x < shape[1]-1 and y >= 0 and y < shape[0]-1:
            if mask is not None:
                if mask[int(round(y)), int(round(x))] == 1:
                    selected.append(True)
                else:
                    selected.append(False)
            else:
                selected.append(True)
        else:
            selected.append(False)
    
    return tab[selected]

#------------------------------------------------------------------------------------------------------------
def addFootprintColumnToCatalog(tab, label, areaMask, wcs):
    """Add `footprint_label` column to the catalog, flagging objects found within the valid area of the given
    mask.

    Args:
        tab (:obj:`astropy.table.Table`): Catalog, as an astropy Table object. Must have columns called
            'RADeg', 'decDeg' that contain object coordinates in decimal degrees.
        label (:obj:`str`): A column named `footprint_label` will be added to the catalog. Objects in the
            catalog that fall within the valid area of the given area mask will have `footprint_label`
            set to True.
        areaMask (:obj:`np.ndarray`): Mask image defining the footprint corresponding to the given WCS.
            Pixels with value = 1 indicate valid area, and pixels with value = 0 are considered to be
            outside the mask.
        wcs (:obj:`astWCS.WCS`): WCS of the area mask that defines the footprint.

    Returns:
        An astropy Table with `footprint_label` column added.

    """

    inMask=np.zeros(len(tab['RADeg'].data), dtype = bool)
    coords=wcs.wcs2pix(tab['RADeg'].data, tab['decDeg'].data)
    coords=np.array(np.round(coords), dtype = int)
    mask1=np.logical_and(coords[:, 0] >= 0, coords[:, 1] >= 0)
    mask2=np.logical_and(coords[:, 0] < areaMask.shape[1], coords[:, 1] < areaMask.shape[0])
    mask=np.logical_and(mask1, mask2)
    inMask[mask]=inMask[mask]+areaMask[coords[:, 1][mask], coords[:, 0][mask]]
    tab['footprint_%s' % (label)]=inMask

    return tab