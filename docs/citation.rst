************
Citation
************

If you use sedcreator for a project that leads to a publication, whether directly or as a dependency of another package, please include the following acknowledgment:

.. code-block:: text

    This research made use of sedcreator (Fedriani et al. 2023, ApJ, 942, 7). ADS url: https://ui.adsabs.harvard.edu/abs/2023ApJ...942....7F/abstract

.. code-block:: text

@ARTICLE{2023ApJ...942....7F,
       author = {{Fedriani}, Rub{\'e}n and {Tan}, Jonathan C. and {Telkamp}, Zoie and {Zhang}, Yichen and {Yang}, Yao-Lun and {Liu}, Mengyao and {De Buizer}, James M. and {Law}, Chi-Yan and {Beltran}, Maria T. and {Rosero}, Viviana and {Tanaka}, Kei E.~I. and {Cosentino}, Giuliana and {Gorai}, Prasanta and {Farias}, Juan and {Staff}, Jan E. and {Whitney}, Barbara},
        title = "{The SOFIA Massive (SOMA) Star Formation Survey. IV. Isolated Protostars}",
      journal = {\apj},
     keywords = {Massive stars, Star formation, Interstellar medium, Jets, Infrared sources, Protostars, Spectral energy distribution, 732, 1569, 847, 870, 793, 1302, 2129, Astrophysics - Astrophysics of Galaxies, Astrophysics - Solar and Stellar Astrophysics},
         year = 2023,
        month = jan,
       volume = {942},
       number = {1},
          eid = {7},
        pages = {7},
          doi = {10.3847/1538-4357/aca4cf},
archivePrefix = {arXiv},
       eprint = {2205.11422},
 primaryClass = {astro-ph.GA},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2023ApJ...942....7F},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}

If you also make use of the SedFitter tool, please also cite the original work for the model grid made by Zhang and Tan (2018).

.. code-block:: text

    @ARTICLE{zhang2018,
        author = {{Zhang}, Yichen and {Tan}, Jonathan C.},
        title = "{Radiation Transfer of Models of Massive Star Formation. IV. The Model Grid and Spectral Energy Distribution Fitting}",
        journal = {\apj},
        keywords = {dust, extinction, ISM: clouds, radiative transfer, stars: formation, stars: massive, Astrophysics - Astrophysics of Galaxies, Astrophysics - Solar and Stellar Astrophysics},
        year = 2018,
        month = jan,
        volume = {853},
        number = {1},
        eid = {18},
        pages = {18},
        doi = {10.3847/1538-4357/aaa24a},
        archivePrefix = {arXiv},
        eprint = {1708.08853},
        primaryClass = {astro-ph.GA},
        adsurl = {https://ui.adsabs.harvard.edu/abs/2018ApJ...853...18Z},
        adsnote = {Provided by the SAO/NASA Astrophysics Data System}
        }
        
sedcreator, in particular SedFluxer, uses multiple resources from Photutils, so please also add acknowledgment their work:

.. code-block:: text

    https://photutils.readthedocs.io/en/stable/citation.html
