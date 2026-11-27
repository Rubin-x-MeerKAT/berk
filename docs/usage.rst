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

Settings for cross-matching
---------------------------

Cross-matching with optical surveys (the ``xmatch`` task, the ``xmatch`` command, or
:mod:`berk.crossmatch`) uses the following settings:

* ``ZCLUSTER_CACHE``: Directory under which optical catalogue queries are cached (in ``zCluster/cache``).
  If not set, ``BERK_ROOT/cache`` is used by the ``xmatch`` task, and ``~/.cache/berk`` otherwise.

* ``DL_USERNAME``, ``DL_PASSWORD``: Data Lab credentials, used to retrieve DECaLS catalogues. These
  are usually not needed: once you have logged in to Data Lab (e.g., with ``datalab login``, or by
  entering your username and password when prompted by **Berk**), the token stored in
  ``~/.datalab`` is re-used for as long as it is valid. If ``DL_USERNAME`` is set, the stored token for
  that user is used. **Berk** will not prompt for a username and password when not running
  interactively (e.g., in a batch job) - so, before submitting an ``xmatch`` job, run ``datalab login``
  (e.g., on a login node).

Each of these can be given as an environment variable or in the **Berk** config file. This is
``~/.config/berk/config.toml`` by default, or can be set with the ``BERK_CONFIG`` environment variable.
Settings go in the ``[crossmatch]`` section, using lower case names, e.g.::

    [crossmatch]
    zcluster_cache = "~/.cache/berk"
    dl_username = "..."

Environment variables take precedence over the config file. As the config file may contain tokens and
passwords, make it readable only by you (``chmod 600 ~/.config/berk/config.toml``). When using
:mod:`berk.crossmatch` from Python, these settings can also be passed as keyword arguments (see
:func:`berk.crossmatch.xmatchRadioOptical`, ``retrieverOptions``).

Reading the config file needs Python 3.11 or later, or the ``tomli`` package.

Rubin catalogues are retrieved from the Rubin Science Platform TAP service, using the access token in
``~/.rsp-tap.token`` (no **Berk** settings are needed). See
https://rsp.lsst.io/guides/auth/creating-user-tokens.html for how to create a token.

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

