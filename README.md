# mosasaurus
Tools for extracting chromatic lightcurves from MultiObject Spectrograph data.

This package is set up so that it should (eventually) be installable via `python setup.py install` and `pip`, but for development purposes I find it most convenient to drop the entire `mosasaurus` directory in my `$PYTHONPATH`, and then I can develop there without having to reinstall the package with every tiny change. The `__init__.py` file in this main directory makes it easier to import in this kludged development mode.

Dependencies include:

+ [zachopy](https://github.com/zkbt/zachopy)
+ [astroquery](https://github.com/astropy/astroquery)
+ Jonathan Irwin's [lfa](https://github.com/mdwarfgeek/lib)
+ (*I think* it no longer depends on pyds9)
+ (Fitting transits will also require [transit](https://github.com/zkbt/transit))

Contributers include:

+ Zach Berta-Thompson
+ Hannah Diamond-Lowe
+ ...
