.. _Usage:

=====
Usage
=====

Initial set-up
==============

Environment variables - only needs to be done once.

First set the ``BERK_ROOT`` environment variable to point to where you
want all the processing, and analysis to happen. Set the ``BERK_MSCACHE``
environment variable to point to a directory where retrieved measurement sets
will be stored.

Workflow
========

Using each of the tasks in **Berk**...

The pipeline runs through a single command with the format::

    berk task captureBlockId

Here, ``captureBlockId`` is used to identify the dataset in the MeerKAT archive,
and ``task`` is one of:

``fetch``:
    Fetch data from the archive, and store it where ``berk`` can
    find it. Here, ``captureBlockId`` should actually be a link to a ``.rdb``
    file on the archive, of the form
    https://archive-gw-1.kat.ac.za/captureBlockId/captureBlockId_sdp_l0.full.rdb?token=longTokenString

``process``:
    Calibrate and image the MeerKAT data using ``Oxkat``. At present, this
    runs to the 2GC (self-calibration) stage, producing continuum images.

``analyse``:
    Produce catalogs from the images and run various tests.


