**Berk** is written in `Python <https://www.python.org/>`_ and is installed in
the standard way. For example, the latest tagged version can be installed using ``pip``:

.. code-block::

   pip install berk

A number of external packages and tools are required - in particular, everything
that Ian Heywood's `Oxkat <https://github.com/IanHeywood/oxkat>`_
and Jonah Wagenveld's `Image-processing <https://github.com/JonahDW/Image-processing>`_ scripts
need (e.g., `PyBDSF <https://pybdsf.readthedocs.io/en/latest/>`_), plus:

* katbeam
* python-casacore
* numpyencoder
* katdal
* rsync

For the moment, not all Python dependencies may be installed automatically by ``pip``.

Note that ``Oxkat`` and the ``Image-processing`` scripts will be installed by
**Berk** when it is first run (specific versions are fetched from GitHub).

If **Berk** has installed correctly, then you should find its command line tools are available, for
example,

.. code-block::

   berk -h

should display a helpful message about the command-line options for the main ``berk`` command.

Note that **Berk** requires some environment variables to be set - see :ref:`Usage`
