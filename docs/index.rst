.. the "raw" directive below is used to hide the title in favor of
   just the logo being visible
.. raw:: html

    <style media="screen" type="text/css">
      h1 {display:none;}
    </style>

.. |br| raw:: html

    <div style="min-height:0.1em;"></div>

***************
sedcreator
***************

.. raw:: html

   <img src="_static/sedcreator_logo.png" onerror="this.src='_static/sedcreator_logo.png'; this.onerror=null;" width="250"/>

**sedcreator** helps you to construct your spectral energy distribution
by providing convinient tools to measure fluxes on any image. It also
provides a set of models to fit your SED with massive star formation models.
It is an open source Python package.


.. Important::
    If you use **sedcreator** for a project that leads to a publication,
    whether directly or as a dependency of another package, please
    include an :doc:`acknowledgment and/or citation <citation>`. Please include also credit
    the python packages that sedcreator depends on.


.. toctree::
   :maxdepth: 1
   :caption: Contents:

Getting Started
===============

.. toctree::
    :maxdepth: 1

    install.rst
    getting_started.rst
    citation.rst
    license.rst

User Documentation
==================

.. toctree::
    :maxdepth: 1

    sedcreator.rst

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`