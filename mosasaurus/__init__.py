__version__ = '0.1'

try:
    __MOSASAURUS_SETUP__
except NameError:
    __MOSASAURUS_SETUP__ = False

if not __MOSASAURUS_SETUP__:
    #__all__ = ['dartmouth','basti','padova',
    #           'Isochrone', 'StarModel', 'BinaryStarModel',
    #           'TripleStarModel']
    from .Reducer import Reducer
