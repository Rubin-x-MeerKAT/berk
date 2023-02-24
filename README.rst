**Hippoxkatapult** provides a pipeline for fetching and processing archival
`MeerKAT <https://skaafrica.atlassian.net/wiki/spaces/ESDKB/overview?homepageId=41025669>`_
data on a compute cluster using `Oxkat <https://github.com/IanHeywood/oxkat>`_
for calibration and imaging, and `PyBDSF <https://www.astron.nl/citt/pybdsf/>`_
for generating object catalogs.

* **License:** `GPL v3 <LICENSE>`_
* **Authors:** Matt Hilton
* **Installation:** ``python setup.py install --user``
* **Support:** Please use the `GitHub issues page <https://github.com/mattyowl/hippoxkatapult/issues>`_,
  and/or contact `Matt Hilton <mailto:matt.hilton@mykolab.com>`_.

**Hippoxkatapult** is being developed and deployed on `Hippo <https://astro.ukzn.ac.za/~hippo/>`_,
the High Performance Computing facility at the University of KwaZulu-Natal, with
the aim of using it to perform a serendipitous survey with MeerKAT, and in
turn to generate useful high-level processed MeerKAT data (initially,
continuum images and catalogs) for use by the community.

This project is at an early stage. At the moment, data can be processed in one step
to the 2GC (self-calibration stage), producing continuum images, and catalogs.
Additional steps will be added to perform automated quality checks, organise the
output, and ingest it into a database.


Usage
-----

Set the ``HIPPOXKATAPULT_ROOT`` environment variable to point to where you want
all the storage, processing, and analysis to happen.

The pipeline runs through a single command with the format::

    hippoxkatapult task captureBlockId

Here, ``captureBlockId`` is used to identify the dataset in the MeerKAT archive,
and ``task`` is one of:

``stage``:
    Fetch data from the archive, and unpack it where ``hippoxkatapult`` can find it.
    (actually fetching data is not implemented yet).

``process``:
    Calibrate and image the MeerKAT data using ``Oxkat``. At present, this
    runs to the 2GC (self-calibration) stage, producing continuum images.

``analyse``:
    Produce catalogs from the images and run various tests.


Notes
-----

* You can use ``process`` immediately, if the measurement set is already
  unpacked (or symlinked) in the ``staging`` directory. The measurement
  set needs to have a name in the format ``captureBlockId_*.ms`` for this
  to work.


Requirements
------------

Everything that `Oxkat <https://github.com/IanHeywood/oxkat>`_ needs plus:

* GNU Screen
* Slurm
* katbeam
* python-casacore
* numpyencoder
* Jonah Wagenveld's `Image-processing <https://github.com/JonahDW/Image-processing>`_ scripts, that use `PyBDSF <https://pybdsf.readthedocs.io/en/latest/>`_.

