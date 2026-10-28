"""

This module contains tools for handling catalogs, which are usually :obj:`astropy.table.Table` objects.

"""

from astLib import astWCS, astCoords
import numpy as np
import astropy.table as atpy
from astropy.coordinates import SkyCoord, Longitude
from astropy.coordinates import match_coordinates_sky, search_around_sky
from astropy import units as u
from . import __version__
import os
from astropy.cosmology import FlatLambdaCDM
from scipy.optimize import brentq

# For adding meta data to output
import datetime

#------------------------------------------------------------------------------------------------------------
def listCatalogInFile(catalogName, outFileName):
    """
    Append or create a text file to record the name of a radio catalogue with no optical counterpart.

    Args:
        catalogName (:obj:`str`): Name or identifier of the radio catalogue without an optical match.
        outFileName (:obj:`str`): Path to the output text file where the catalogue name will be recorded.

    Returns:
        :obj:`int`: Returns 0 upon successful write operation.
    """

    with open(outFileName, 'a', encoding='utf8') as outFile:
        outFile.write("%s\n" %(catalogName))
    return 0

#------------------------------------------------------------------------------------------------------------
def IsFieldContinuous(RAValues):
    """
    Check whether a set of right ascension values forms a continuous field
    without wrapping around the 0°/360° boundary.

    Args:
        RAValues (:obj:`array-like`): Array or list of RA values in degrees.

    Returns:
        bool:
            - True if the RA values form a continuous field (no jump > 180°).
            - False if there is a large jump indicating the field wraps around
              the RA=0°/360° boundary.
    """

    diffsRA = np.diff(np.sort(RAValues))
    maxGap = np.max(diffsRA)
    if maxGap > 180.0:
        return False
    else:
        return True

#------------------------------------------------------------------------------------------------------------
def detectWrapAngle(ra_values):
    """
    If the RA range spans the 0/360 boundary, wrap at 180.
    Detection: min < threshold_low AND max > threshold_high.
    """
    ra = np.asarray(ra_values)
    ra = ra % 360
    ra_min, ra_max = ra.min(), ra.max()

    if ra_min < 10 and ra_max > 350:
        return 180.0  # field straddles 0/360 → centre around 0
    else:
        return 360.0  # normal case

#------------------------------------------------------------------------------------------------------------
def fixRA(table, raCol='RA', wrapAngle=None):
    """Returns table with corrected RA wrap.

    Args:
        table (:obj:`~astropy.table.Table`): Input table with RA values.
        raCol (:obj:`str`, optional): Name of the RA column. Default is 'RA'.
        wrapAngle (:obj:`float`, optional): Angle at which to wrap RA in degrees. If None, it is auto-detected from the data.

    Returns:
        :obj:`~astropy.table.Table`: Table with RA values wrapped to [0, 360) range.
    """
    fixTable = table.copy()

    if wrapAngle is None:
        wrapAngle = detectWrapAngle(table[raCol])
        print(f"Auto-detected wrap angle: {wrapAngle} degrees based on RA distribution.")

    fixTable[raCol] = Longitude(table[raCol], unit=u.deg, wrap_angle=wrapAngle * u.deg).value
    if IsFieldContinuous(fixTable[raCol]) is False:
        newWrapAngle = 180.0 if wrapAngle == 360 else 360
        fixTable[raCol] = Longitude(table[raCol], unit=u.deg, wrap_angle=newWrapAngle * u.deg).value
        if IsFieldContinuous(fixTable[raCol]) is False:
            print("\nIssue with RA wrap fixing for table of " %table.meta['INIMAGE'])

    return fixTable

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
        if np.isnan(x) is True or np.isnan(y) is True:
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

#------------------------------------------------------------------------------------------------------------
def calculateRadioLum(fluxJy, redshift, spectralIndex=0.7, cosmology=None):
    """
    Computes the radio luminosity of a source from its observed flux density and redshift,
    applying a K-correction assuming a power-law spectrum.

    Args:
        fluxJy (float or np.ndarray): Observed flux density in Jansky (Jy).
        redshift (float or np.ndarray): Redshift of the source.
        spectralIndex (float): Spectral index alpha, defined such that S_nu ~ nu^{-alpha}.
            Default is 0.7.
        cosmology (astropy.cosmology instance, optional): Cosmology to use for luminosity
            distance calculation. If None, defaults to FlatLambdaCDM with H0=70 km/s/Mpc
            and Om0=0.3.

    Returns:
        float or np.ndarray: Radio luminosity in W/Hz.

    Notes:
        The K-correction applied is (1+z)^{alpha-1}, appropriate for a power-law spectrum
        S_nu ~ nu^{-alpha} observed at a fixed frequency.
    """

    if cosmology is None:
        cosmo = FlatLambdaCDM(H0=70 * u.km / u.s / u.Mpc, Om0=0.3)
    else:
        cosmo = cosmology

    z = redshift
    alpha = spectralIndex

    # Convert flux to W/m^2/Hz
    S_Wm2Hz = (fluxJy * u.Jy).to(u.W / u.m**2 / u.Hz)

    # Luminosity distance
    DL = cosmo.luminosity_distance(z)      # in Mpc
    DL_si = DL.to(u.m)                     # convert to metres

    # Compute luminosity (with K-correction)
    L = 4 * np.pi * DL_si**2 * S_Wm2Hz * (1 + z)**(alpha - 1)

    # Convert to W/Hz
    L_WHz = L.to(u.W / u.Hz).value

    return L_WHz

#------------------------------------------------------------------------------------------------------------
def calculateFluxDensityFromRadLum(luminosityWHz, redshift, spectralIndex=0.7, cosmology=None):
    """
    Computes the observed radio flux density of a source from its radio luminosity
    and redshift, applying the inverse K-correction assuming a power-law spectrum.

    Args:
        luminosityWHz (float or np.ndarray): Radio luminosity in W/Hz.
        redshift (float or np.ndarray): Redshift of the source.
        spectralIndex (float): Spectral index alpha, defined such that S_nu ~ nu^{-alpha}.
            Default is 0.7.
        cosmology (astropy.cosmology instance, optional): Cosmology to use for luminosity
            distance calculation. If None, defaults to FlatLambdaCDM with H0=70 km/s/Mpc
            and Om0=0.3.

    Returns:
        float or np.ndarray: Observed flux density in Jansky (Jy).
    """

    if cosmology is None:
        cosmo = FlatLambdaCDM(H0=70 * u.km / u.s / u.Mpc, Om0=0.3)
    else:
        cosmo = cosmology

    z = redshift
    alpha = spectralIndex

    # Attach units to the input luminosity
    L_si = luminosityWHz * (u.W / u.Hz)

    # Calculate luminosity distance and convert to meters
    DL = cosmo.luminosity_distance(z)
    DL_si = DL.to(u.m)

    # Compute observed flux density in W/m^2/Hz (incorporating inverse K-correction)
    S_Wm2Hz = L_si / (4 * np.pi * DL_si**2 * (1 + z)**(alpha - 1))

    # Convert back to Jansky and extract the numerical value
    fluxJy = S_Wm2Hz.to(u.Jy).value

    return fluxJy

#------------------------------------------------------------------------------------------------------------
def fluxDensityAtRedshift_uJy(z, Lrest_WHz, alpha=0.7, cosmology=None):
    """
    Computes the expected observed flux density of a source at a given redshift,
    assuming a power-law spectrum and applying a K-correction.

    Args:
        z (float or np.ndarray): Redshift at which to evaluate the flux density.
        Lrest_WHz (float or np.ndarray): Rest-frame radio luminosity in W/Hz.
        alpha (float): Spectral index, defined such that S_nu ~ nu^{-alpha}.
            Default is 0.7.
        cosmology (astropy.cosmology instance, optional): Cosmology to use for
            luminosity distance. If None, defaults to FlatLambdaCDM with
            H0=70 km/s/Mpc and Om0=0.3.

    Returns:
        float or np.ndarray: Expected flux density in microjansky (uJy).
    """

    if cosmology is None:
        cosmology = FlatLambdaCDM(H0=70 * u.km / u.s / u.Mpc, Om0=0.3)

    DL = cosmology.luminosity_distance(z)      # in Mpc
    DL_si = DL.to(u.m)                     # convert to metres

    S_Wm2Hz = Lrest_WHz / (4 * np.pi * DL_si**2 * (1 + z)**(alpha - 1))
    S_Jy = S_Wm2Hz * (u.W / u.m**2 / u.Hz).to(u.Jy)
    S_uJy = S_Jy.value*1E6

    return S_uJy

#------------------------------------------------------------------------------------------------------------
def calculateZmaxIterative(galRedshift, galLum_WHz, Slim_uJy, alpha=0.7, zmaxLimit=10.0, step=None, cosmology=None):
    """
    Finds the maximum redshift at which a source of given luminosity would remain
    detectable above a survey flux limit, using a linear forward-stepping search.

    Starting from the source's observed redshift, the flux density is recomputed
    at each step of size `step` and compared to the survey limit. The search stops
    as soon as the flux density drops below Slim_uJy, or once zmaxLimit is reached.

    Args:
        galRedshift (float): Observed redshift of the source. Search begins here
            since zmax >= galRedshift by definition.
        galLum_WHz (float): Rest-frame radio luminosity of the source in W/Hz.
        Slim_uJy (float): Survey flux density limit in microjansky (uJy),
            typically N-sigma * RMS of the image.
        alpha (float): Spectral index, defined such that S_nu ~ nu^{-alpha}.
            Default is 0.7.
        zmax_limit (float): Hard upper bound on redshift to search. If the source
            remains detectable at this redshift, zmax_limit is returned.
            Default is 10.0.
        step (float): Step size in redshift for the iterative search. Default is 0.001.
        cosmology (astropy.cosmology instance, optional): Cosmology to use. If None,
            defaults to FlatLambdaCDM with H0=70 km/s/Mpc and Om0=0.3.

    Returns:
        float: Maximum redshift zmax at which the source flux equals Slim_uJy.
            Returns zmax_limit if the source is detectable across the full search range.
    """

    if cosmology is None:
        cosmology = FlatLambdaCDM(H0=70 * u.km / u.s / u.Mpc, Om0=0.3)

    if step is None:
        step = 0.001

    zIter = galRedshift
    SIter = fluxDensityAtRedshift_uJy(zIter, galLum_WHz, alpha=alpha, cosmology=cosmology)

    while SIter >= Slim_uJy: # and zIter < zmaxLimit:
        zIter += step
        SIter = fluxDensityAtRedshift_uJy(zIter, galLum_WHz, alpha=alpha, cosmology=cosmology)

    return zIter

#------------------------------------------------------------------------------------------------------------
def calculateZmax(galRedshift, galLum_WHz, Slim_uJy, alpha=0.7, zmaxLimit=10.0, cosmology=None):
    """
    Finds the maximum redshift at which a source of given luminosity would remain
    detectable above a survey flux limit, using Brent's root-finding method.

    Args:
        galRedshift (float): Observed redshift of the source. Search begins here
            since zmax >= galRedshift by definition.
        galLum_WHz (float): Rest-frame radio luminosity of the source in W/Hz.
        Slim_uJy (float): Survey flux density limit in microjansky (uJy),
            typically N-sigma * RMS of the image.
        alpha (float): Spectral index, defined such that S_nu ~ nu^{-alpha}.
            Default is 0.7.
        zmax_limit (float): Hard upper bound on redshift to search. If the source
            remains detectable at this redshift, zmax_limit is returned.
            Default is 10.0.
        cosmology (astropy.cosmology instance, optional): Cosmology to use. If None,
            defaults to FlatLambdaCDM with H0=70 km/s/Mpc and Om0=0.3.

    Returns:
        float: Maximum redshift zmax at which the source flux equals Slim_uJy.
            Returns zmax_limit if the source is detectable across the full search range.
    """

    if cosmology is None:
        cosmology = FlatLambdaCDM(H0=70 * u.km / u.s / u.Mpc, Om0=0.3)

    def f(z):
        return fluxDensityAtRedshift_uJy(z, galLum_WHz, alpha=alpha, cosmology=cosmology) - Slim_uJy

    if f(zmaxLimit) >= 0:
        return zmaxLimit

    if f(galRedshift) < 0:
        return galRedshift

    return brentq(f, galRedshift, zmaxLimit)

#------------------------------------------------------------------------------------------------------------
def calculateComovingVolBetweenZ_h3Mpc3(skyArea, zMin, zMax, cosmology=None):
    """
    Computes the comoving volume of a spherical shell between two redshifts,
    scaled to the solid angle subtended by a survey field.

    Args:
        skyArea (float): Sky area of the survey field in square degrees.
        zMin (float): Minimum redshift of the shell.
        zMax (float): Maximum redshift of the shell.
        cosmology (astropy.cosmology instance, optional): Cosmology to use for
            comoving volume calculation. If None, defaults to FlatLambdaCDM
            with H0=70 km/s/Mpc and Om0=0.3.

    Returns:
        float: Comoving volume of the shell subtended by the survey field,
            in units of h^-3 Mpc^3.
    """
    if cosmology is None:
        cosmology = FlatLambdaCDM(H0=70 * u.km / u.s / u.Mpc, Om0=0.3)

    totalSkyArea = (4 * np.pi * u.sr).to(u.deg**2).value # (~ 41252.96 deg^2); full sky = 4pi steradians

    vmax = cosmology.comoving_volume(zMax).value  # Mpc^3
    vmin = cosmology.comoving_volume(zMin).value  # Mpc^3

    volume = (vmax - vmin) * skyArea / totalSkyArea

    volume_h3Mpc3 = volume * (cosmology.h**3)  # Convert to h^-3 Mpc^3

    return volume_h3Mpc3

#------------------------------------------------------------------------------------------------------------
def getUnitlessValues(col):
    """
    Returns unitless numeric values from an astropy Table Column.
    If the column has a unit (not None), returns .value.
    Otherwise, returns the column as-is.
    """
    if hasattr(col, 'unit') and col.unit is not None:
        return col.value
    else:
        return col

#------------------------------------------------------------------------------------------------------------
def removeDuplicateSources(tab, matchRadius_arcsec=6.0, raCol='RA', decCol='DEC', fluxCol='Total_flux', fluxErrCol='E_Total_flux'):
    """
    Removes duplicate sources from a merged catalogue arising from overlapping
    image footprints. For each pair of sources within matchRadius_arcsec,
    retains the detection from the deeper image (lower RMS).

    Args:
        tab (astropy.table.Table): Merged source catalogue with 'RA', 'DEC',
            and 'radCatPath' columns.
        matchRadius_arcsec (float): Matching radius in arcseconds. Should be
            approximately one beam FWHM. Default is 6.0.
        raCol (str): Column name for right ascension in degrees. Default is 'RA_rad'.
        decCol (str): Column name for declination in degrees. Default is 'DEC_rad'.
        fluxCol (str): Column name for total flux density. Default is 'Total_flux'.
        fluxErrCol (str): Column name for flux density uncertainty. Default is 'E_Total_flux'.

    Returns:
        astropy.table.Table: Deduplicated catalogue.
    """
    tab = tab.copy()
    tab['FluxErrRatio'] = tab[fluxCol] / tab[fluxErrCol]

    raArray = getUnitlessValues(tab[raCol])
    decArray = getUnitlessValues(tab[decCol])
    coords = SkyCoord(ra=raArray*u.deg, dec=decArray*u.deg)
    idx1, idx2, sep, _ = search_around_sky(coords, coords, matchRadius_arcsec*u.arcsec)

    remove = np.zeros(len(tab), dtype=bool)
    processed = set()
    for i, j in zip(idx1, idx2):
        if i == j:
            continue
        pair = frozenset([i, j])
        if pair in processed:
            continue
        processed.add(pair)
        if tab['FluxErrRatio'][i] <= tab['FluxErrRatio'][j]:
            remove[i] = True
        else:
            remove[j] = True

    n_dupes = remove.sum()
    print("Removed %d duplicates out of %d within %0.2f arcsec" % (n_dupes, len(tab), matchRadius_arcsec))
    return tab[~remove]

#------------------------------------------------------------------------------------------------------------
def convertFluxFreq(flux1, nu1, nu2, alpha=0.7):
    """
    Converts radio flux densities between two observing frequencies assuming
    a power-law radio spectrum.

    Args:
        flux1 (float or np.ndarray):
            Flux density at the original frequency. Can be a scalar or array.
            Units may be any consistent flux density unit (e.g. Jy, mJy, uJy),
            as the conversion is unit-independent.
        nu1 (float or np.ndarray):
            Original observing frequency in the same units as ``nu2``
            (e.g. Hz, MHz or GHz).
        nu2 (float or np.ndarray):
            Target observing frequency in the same units as ``nu1``.
        alpha (float, optional):
            Radio spectral index, defined such that

                S_nu ∝ nu^{-alpha}

            Default is 0.7.

    Returns:
        float or np.ndarray:
            Flux density at the target frequency, in the same units as
            ``flux1``.
    """

    flux1 = np.asarray(flux1)
    nu1 = np.asarray(nu1)

    flux2 = flux1 * (nu2 / nu1) ** (-alpha)

    return flux2

#------------------------------------------------------------------------------------------------------------
def convertLuminosityFreq(luminosity1, nu1, nu2, alpha=0.7):
    """
    Converts radio luminosities between two frequencies assuming
    a power-law radio spectrum.

    Args:
        luminosity1 (float or np.ndarray):
            Monochromatic radio luminosity at the original frequency
            (typically in W/Hz).
        nu1 (float or np.ndarray):
            Original frequency in the same units as ``nu2``.
        nu2 (float or np.ndarray):
            Target frequency in the same units as ``nu1``.
        alpha (float, optional):
            Radio spectral index, defined such that

                L_nu ∝ nu^{-alpha}

            Default is 0.7.

    Returns:
        float or np.ndarray:
            Monochromatic radio luminosity at the target frequency, in the
            same units as ``luminosity1``.
    """

    luminosity1 = np.asarray(luminosity1)
    nu1 = np.asarray(nu1)

    luminosity2 = luminosity1 * (nu2 / nu1) ** (-alpha)

    return luminosity2