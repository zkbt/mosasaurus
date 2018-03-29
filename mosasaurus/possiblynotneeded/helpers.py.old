from mosasaurus.Reducer import Reducer
from mosasaurus.Cube import Cube


def do(files=None):
    for f in files:
        r = Reducer(f)
        r.reduce()
        c = Cube(f)
        c.populate(remake=False, visualize=False)
        c.movieCube(stride=1)
