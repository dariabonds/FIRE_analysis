import yt
from yt.utilities.math_utils import ortho_find

ds = yt.load('../m12i_res56000_md/snapshot_600.hdf5')
_, c = ds.find_max(('gas', 'density'))
#p = yt.ProjectionPlot(ds, 'x', ('gas', 'density'), center=c, width=(100, 'kpc'))
#p.save()

sp = ds.sphere(c, (25.0, 'kpc'))
L = sp.quantities.angular_momentum_vector()
L, edge1, edge2 = ortho_find(L)
p = yt.OffAxisProjectionPlot(ds, edge2, ('gas', 'density'), center=c, width=(100, 'kpc'), north_vector=L)
p.save()
