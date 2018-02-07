# mosasaurus
Tools for extracting chromatic lightcurves from MultiObject Spectrograph data.

This package is set up so that it should (eventually) be installable via `python setup.py install` and `pip`, but for development purposes I find it most convenient to drop the entire `mosasaurus` directory in my `$PYTHONPATH`, and then I can develop there without having to reinstall the package with every tiny change. The `__init__.py` file in this main directory makes it easier to import in this kludged development mode.

*These tools are still under active development. There's no promise they'll do something useful for you if you download them. Since January 2018, the `reorganization` branch is pretty darn far ahead of `master`. It'll be merged in soon.*

Dependencies include:

+ [craftroom](https://github.com/zkbt/craftroom)
+ [astroquery](https://github.com/astropy/astroquery)
+ astropy 1.2.1
+ (Fitting transits will also require [transit](https://github.com/zkbt/transit))

Contributers include:

+ Zach Berta-Thompson
+ Hannah Diamond-Lowe
+ ...
