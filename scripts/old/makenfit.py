filename = 'wasp94_140801.obs'
import Cube
c = Cube.Cube(filename)
c.makeLCs(binsize=1000)
import transmission
ts = transmission.load(filename, 1000)
ts.createMask()
ts.fitRigid(plot=True)
