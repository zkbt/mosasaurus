from setuptools import setup, find_packages
import os,sys

def readme():
    with open('README.md') as f:
        return f.read()

# Hackishly inject a constant into builtins to enable importing of the
# package before the library is built.
import sys
if sys.version_info[0] < 3:
    import __builtin__ as builtins
else:
    import builtins
builtins.__MOSASAURUS_SETUP__ = True
import mosasaurus
version = mosasaurus.__version__

"""
# Publish the library to PyPI.
if "publish" in sys.argv[-1]:
    os.system("python setup.py sdist upload")
    sys.exit()

# Push a new tag to GitHub.
if "tag" in sys.argv:
    os.system("git tag -a {0} -m 'version {0}'".format(version))
    os.system("git push --tags")
    sys.exit()
"""

setup(name = "mosasaurus",
    version = version,
    description = "Tools for extracting chromatic light curves from MultiObject Spectra.",
    long_description = readme(),
    author = "Zachory K. Berta-Thompson",
    author_email = "zach.bertathompson@colorado.edu",
    url = "https://github.com/zkbt/mosasaurus",
    packages = find_packages(),
    package_data = {'':[]},
    scripts = [],
    classifiers=[
      'Intended Audience :: Science/Research',
      'Programming Language :: Python',
      'Topic :: Scientific/Engineering :: Astronomy'
      ],
    install_requires=['astropy>=1.2.1','emcee>=2.0',
                      'numpy>=1.9', 'triangle_plot'],
    zip_safe=False
)
