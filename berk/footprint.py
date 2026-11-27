"""

Sky footprint of the MeerKAT images, stored as a MOC (Multi-Order Coverage map; see
https://www.ivoa.net/documents/MOC/), for cutting catalogues to the area covered by the images and for
calculating the total sky area covered.

The footprint is made by rasterising the non-blank (i.e., not NaN) pixels of each image onto HEALPix cells.
A cell is in the footprint if its centre falls on a non-blank pixel of any image. The footprint of several
images is the union of their footprints, so overlapping pointings are only counted once in the area.

The default resolution is HEALPix order 14 (NSIDE = 16384; cells about 12.9 arcsec across), which is much
finer than the size of a pointing, but coarser than the image pixels, so that making the footprint is fast.

"""

import os
import glob
import numpy as np
import astropy.io.fits as pyfits
import astropy.units as u
from astropy.coordinates import Angle
from astropy.wcs import WCS
from astropy_healpix import HEALPix
from . import images

IMAGE_PATTERN = "pbcorr_trim_img_*_pcalmask-MFS-image.fits"
DEFAULT_ORDER = 14
FULL_SKY_SQ_DEG = 4*np.pi*(180/np.pi)**2

#------------------------------------------------------------------------------------------------------------
def _importMOC():
    """Imports mocpy's MOC class (only needed for the footprint).

    """
    try:
        from mocpy import MOC
    except ImportError:
        raise ImportError("The mocpy package is needed for footprints - install it with: pip install mocpy")
    return MOC

#------------------------------------------------------------------------------------------------------------
def _readImage(imgFileName):
    """Returns the 2D image data and celestial WCS of a FITS image (degenerate frequency and Stokes axes are
    removed).

    """

    with pyfits.open(imgFileName, memmap = True) as img:
        header = img[0].header
        data = img[0].data
        while data.ndim > 2:
            data = data[0]
        data = np.array(data)
    wcs = WCS(header).celestial

    return data, wcs

#------------------------------------------------------------------------------------------------------------
def getImageHealpixCells(imgFileName, order = DEFAULT_ORDER):
    """Returns the HEALPix cells (NESTED scheme) whose centres fall on non-blank (not NaN) pixels of the image.

    Args:
        imgFileName (:obj:`str`): Path to a FITS image.
        order (:obj:`int`, optional): HEALPix order (NSIDE = 2**order).

    Returns:
        :obj:`np.ndarray`: HEALPix cell indices (NESTED scheme) at the given order.

    """

    data, wcs = _readImage(imgFileName)
    ny, nx = data.shape

    hp = HEALPix(nside = 2**order, order = 'nested', frame = 'icrs')

    # Candidate cells: those containing a subsample of the non-blank pixels, taken every `step` pixels, where
    # the spacing is at most a third of the cell size (so no cell with its centre on a non-blank pixel is
    # missed, other than perhaps tiny slivers at the edges). Done in blocks of rows to limit memory use.
    # (This is much faster than astropy_healpix's cone search at high resolution.)
    pixScaleDeg = np.sqrt(np.abs(np.linalg.det(wcs.pixel_scale_matrix)))
    step = max(1, int(hp.pixel_resolution.to(u.deg).value/3/pixScaleDeg))
    xs = np.arange(0, nx, step)
    candidates = []
    rowsPerBlock = max(1, 2000000//len(xs))
    ys = np.arange(0, ny, step)
    for i in range(0, len(ys), rowsPerBlock):
        yy, xx = np.meshgrid(ys[i:i+rowsPerBlock], xs, indexing = 'ij')
        ok = np.isfinite(data[yy, xx])
        if ok.sum() == 0:
            continue
        candidates.append(np.unique(hp.skycoord_to_healpix(wcs.pixel_to_world(xx[ok], yy[ok]))))
    if len(candidates) == 0:
        return np.array([], dtype = np.int64)
    cells = np.unique(np.concatenate(candidates))

    # Keep the cells whose centres fall on non-blank pixels
    x, y = wcs.world_to_pixel(hp.healpix_to_skycoord(cells))
    x = np.round(x)
    y = np.round(y)
    inImage = np.isfinite(x) & np.isfinite(y) & (x >= 0) & (x < nx) & (y >= 0) & (y < ny)
    cells = cells[inImage]
    valid = np.isfinite(data[y[inImage].astype(int), x[inImage].astype(int)])

    return cells[valid]

#------------------------------------------------------------------------------------------------------------
def findImages(imageDir, pattern = IMAGE_PATTERN, band = None):
    """Returns the paths of the images in a directory that match the given file name pattern.

    Args:
        imageDir (:obj:`str`): Directory containing the images.
        pattern (:obj:`str`, optional): File name pattern (glob) of the images.
        band (:obj:`str`, optional): If given (e.g., 'L', 'UHF'), only images in this band (found from the
            frequency in the image header) are returned.

    Returns:
        Sorted list of image paths.

    """

    imgPaths = sorted(glob.glob(os.path.join(imageDir, pattern)))
    if band is not None:
        imgPaths = [p for p in imgPaths if getImageBand(p) == band]

    return imgPaths

#------------------------------------------------------------------------------------------------------------
def getImageBand(imgFileName):
    """Returns the band (e.g., 'L', 'UHF', 'S') of an image, found from the frequency in its header.

    """

    with pyfits.open(imgFileName) as img:
        freqGHz = img[0].header['CRVAL3']/1e9

    return images.getBandKey(freqGHz)

#------------------------------------------------------------------------------------------------------------
def makeFootprintMOC(imgPaths, order = DEFAULT_ORDER, verbose = True):
    """Makes the footprint of a set of images, as a MOC (the union of the non-blank areas of all the images).

    Args:
        imgPaths (:obj:`list`): Paths of FITS images.
        order (:obj:`int`, optional): HEALPix order (NSIDE = 2**order) at which the images are rasterised.
        verbose (:obj:`bool`, optional): If True, print progress.

    Returns:
        :obj:`mocpy.MOC`: The footprint.

    """

    MOC = _importMOC()
    if len(imgPaths) == 0:
        raise ValueError("No images given to make the footprint from")

    cellsList = []
    for i, imgPath in enumerate(imgPaths):
        if verbose:
            print("[%d/%d] %s" % (i+1, len(imgPaths), os.path.basename(imgPath)))
        cellsList.append(getImageHealpixCells(imgPath, order = order))
    cells = np.unique(np.concatenate(cellsList)).astype(np.uint64)

    return MOC.from_healpix_cells(cells, np.full(len(cells), order, dtype = np.uint8), max_depth = order)

#------------------------------------------------------------------------------------------------------------
def readFootprintMOC(mocFileName):
    """Reads a footprint MOC from a FITS file.

    """

    MOC = _importMOC()

    return MOC.from_fits(mocFileName)

#------------------------------------------------------------------------------------------------------------
def getFootprintAreaSqDeg(moc):
    """Returns the sky area (square degrees) covered by a footprint MOC.

    """

    return moc.sky_fraction*FULL_SKY_SQ_DEG

#------------------------------------------------------------------------------------------------------------
def getInFootprintMask(tab, moc, RACol = 'RADeg', decCol = 'decDeg'):
    """Returns a boolean array that is True for the objects in a catalogue that lie within the footprint.

    Args:
        tab (:obj:`astropy.table.Table`): Catalogue.
        moc (:obj:`mocpy.MOC`): Footprint.
        RACol, decCol (:obj:`str`, optional): Position columns (decimal degrees).

    Returns:
        :obj:`np.ndarray`: Boolean mask.

    """

    lon = Angle(np.array(tab[RACol], dtype = float), unit = "deg")
    lat = Angle(np.array(tab[decCol], dtype = float), unit = "deg")

    return np.array(moc.contains_lonlat(lon, lat), dtype = bool)

#------------------------------------------------------------------------------------------------------------
def getTabInFootprint(tab, moc, RACol = 'RADeg', decCol = 'decDeg'):
    """Returns the subset of a catalogue that lies within the footprint.

    Args:
        tab (:obj:`astropy.table.Table`): Catalogue.
        moc (:obj:`mocpy.MOC`): Footprint.
        RACol, decCol (:obj:`str`, optional): Position columns (decimal degrees).

    Returns:
        :obj:`astropy.table.Table`: The objects in `tab` within the footprint.

    """

    return tab[getInFootprintMask(tab, moc, RACol = RACol, decCol = decCol)]
