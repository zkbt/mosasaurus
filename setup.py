#!/usr/bin/env python

# The template for this setup.py came from Tim Morton,
# who I understand took it from Dan F-M. And then Geert
# Barentsen and Christina Hedges helped explain a few
# more neat tips. Thanks all!


import os, sys
from setuptools import setup, find_packages

# Prepare and send a new release to PyPI
if "release" in sys.argv[-1]:
    os.system("python setup.py sdist")
    # uncomment this to test out on test.pypi.com/project/tess-zap
    # os.system("twine upload --repository-url https://test.pypi.org/legacy/ dist/*")
    os.system("twine upload dist/*")
    os.system("rm -rf dist/tesszap*")
    sys.exit()

# a little kludge to be able to get the version number from the __init__.py
if sys.version_info[0] < 3:
    import __builtin__ as builtins
else:
    import builtins
builtins.__MOSASAURUS_SETUP__ = True
import mosasaurus
version = mosasaurus.__version__

# pull the long description from the readmedef readme():
def readme():
    with open('README.md') as f:
        return f.read()

setup(name = "mosasaurus",
    version = version,
    description = "Tools for extracting chromatic light curves from MultiObject Spectra.",
    long_description = readme(),
    author = "Zach Berta-Thompson",
    author_email = "zach.bertathompson@colorado.edu",
    url = "https://github.com/zkbt/mosasaurus",
    packages = find_packages(),
    package_data = {'mosasaurus': [ '../data/LDSS3C/vph-red/*',
                                    '../data/LDSS3C/vph-all/*']},
    include_package_data=True,
    scripts = [],
    classifiers=[
      'Intended Audience :: Science/Research',
      'Programming Language :: Python',
      'Topic :: Scientific/Engineering :: Astronomy'
      ],
    install_requires=['numpy', 'astropy', 'astroquery', 'scipy', 'matplotlib',  'craftroom'], #'emcee', 'corner',
    dependency_links=['git+https://github.com/zkbt/craftroom.git@master#egg=craftroom'],
    zip_safe=False,
    license='MIT',
)
