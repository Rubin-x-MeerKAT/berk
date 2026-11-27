.. _Commands:

=============
Berk Commands
=============

Command-line programs in the **Berk** package are described below.


.. _berkCommand:
    
berk
----

.. argparse::
   :filename: ../bin/berk
   :func: makeParser
   :prog: berk
   
   :program:`berk` can perform many tasks, such as fetching data from
   the MeerKAT archive ; calibration and imaging ; cataloging ; and
   distribution of data products.


.. _berkFootprintCommand:

berk_footprint
--------------

.. argparse::
   :filename: ../bin/berk_footprint
   :func: make_parser
   :prog: berk_footprint

   :program:`berk_footprint` makes the sky footprint of a set of **Berk**-processed MeerKAT images
   (by default, files matching ``pbcorr_trim_img_*_pcalmask-MFS-image.fits`` in a directory), stored as a
   `MOC <https://www.ivoa.net/documents/MOC/>`_. The footprint is the union of the non-blank pixels of the
   images, rasterised onto HEALPix cells (order 14, about 12.9 arcsec, by default), so overlapping pointings
   are only counted once in its area. It can then cut any catalogue to the area covered by the images.
