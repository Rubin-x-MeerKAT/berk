.. image:: https://readthedocs.org/projects/berk/badge/?version=latest

**Berk** is a package for fetching and processing
`MeerKAT <https://skaafrica.atlassian.net/wiki/spaces/ESDKB/overview?homepageId=41025669>`_
data on a compute cluster using `Oxkat <https://github.com/IanHeywood/oxkat>`_
for calibration and imaging, and `PyBDSF <https://www.astron.nl/citt/pybdsf/>`_
for generating object catalogs.

* **Documentation:** https://berk.readthedocs.io
* **License:** `GPL v3 <LICENSE>`_
* **Authors:** Matt Hilton
* **Installation:** ``pip install berk``
* **Support:** Please use the `GitHub issues page <https://github.com/mattyowl/berk/issues>`_,
  and/or contact `Matt Hilton <mailto:matt.hilton@mykolab.com>`_.

**Berk** should run on any compute cluster that runs either Slurm or PBS as a
workload manager. It is is being developed and deployed on `Hippo <https://astro.ukzn.ac.za/~hippo/>`_,
the High Performance Computing facility at the University of KwaZulu-Natal,
and at the `Centre for High Performance Computing <https://www.chpc.ac.za/>`_,
with the aim of using it to perform a serendipitous survey with MeerKAT, and in
turn to generate useful high-level processed MeerKAT data (initially,
continuum images and catalogs) for use by the community.

**Berk** is under active development, and not all documentation may be up to date.
Not all needed features are implemented yet. However, at the moment, data
can be fetched from the archive, processed to images and catalogs, and collected
from multiple compute clusters using only a few shell commands.

**Berk** is named after the faithful servant of `The Thing Upstairs <https://en.wikipedia.org/wiki/The_Trap_Door>`_.
