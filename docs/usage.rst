.. _Usage:

=====
Usage
=====


Initial set-up
==============

**Berk** requires the following environment variables to be set:

* ``BERK_ROOT``: Path to the root directory in which **Berk** will work. Inside this directory,
  several other directories will be created (e.g., ``processing``, ``products``).

* ``BERK_MSCACHE``: Path to a directory where downloaded measurement sets will be stored.

* ``BERK_PLATFORM``: Needed for setting the compute cluster environment (workload manager to use etc.).
  Currently, this should be set to either ``'chpc'`` or ``'hippo'``.

* ``BERK_NODES_FILE``: Used only by the ``collect`` task, this should be a path (or URL) to a file
  that specifies ``BERK_ROOT`` directories on remote computers that you wish to ``collect`` from
  (see below).

* ``BERK_INFO_FILE``: Used by the ``list`` task, this should be a URL that points to a file containing
  information on the processed fields (e.g., ``products/images.fits``, as produced by the ``builddb``
  task).

.. _Workflow:

Workflow
========

Everything that the user can ask **Berk** to do is controlled with the :ref:`berkCommand` command.

A typical workflow would consist of the following (each of these tasks must be completed in this order):

``berk fetch -o https://archive-gw-1.kat.ac.za/captureBlockId/captureBlockId_sdp_l0.full.rdb?token=longTokenString``:
    Fetch data from the MeerKAT archive, and store it where **Berk** can
    find it. Here, the ``-o`` argument should be provided with a link to a ``.rdb``
    file on the MeerKAT archive, as indicated above. You will probably need to run this
    task from within, e.g., GNU Screen.

``berk process -o captureBlockId``:
    Submit jobs to calibrate and image the MeerKAT observations (identified by ``captureBlockId``) using ``Oxkat``.
    At present, this runs to the 2GC (self-calibration) stage, producing continuum images.

``berk analyse -o captureBlockId``:
    Produce primary beam corrected images and create catalogs from them using ``PyBDSF``.

``berk collect``:
    Collect data products (primary beam corrected images and catalogs) and copy them into the ``BERK_ROOT/products``
    directory. The ``BERK_NODES`` file is used to specify locations of the ``BERK_ROOT`` directories on remote
    machines (one per line, in the format ``user@hostname:$PATH_TO_BERK_ROOT``), from which data products will be
    collected using ``rsync``.

Other tasks include:

``berk list``:
    List observations (i.e., measurement sets) found in the ``BERK_MSCACHE`` directory, by their ``captureBlockId``.

