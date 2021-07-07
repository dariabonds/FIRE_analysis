import yt
from yt.utilities.math_utils import ortho_find

#load in dataset
ds = yt.load('../m12i_res56000_md/snapshot_600.hdf5')
#find galaxy center
_, c = ds.find_max(('gas', 'density'))
#plot galaxy face on
#p = yt.ProjectionPlot(ds, 'x', ('gas', 'density'), center=c, width=(100, 'kpc'))
#p.save()

#create sphere
sp = ds.sphere(c, (25.0, 'kpc'))
#find max angular momentum
L = sp.quantities.angular_momentum_vector()
#find orthogonal vectors
L, edge1, edge2 = ortho_find(L)
#plot galaxy edge on
#p = yt.OffAxisProjectionPlot(ds, edge2, ('gas', 'density'), center=c, width=(100, 'kpc'), north_vector=L)
#p.save()

#create a phase plot
sp2 = ds.sphere(c, (250.0, 'kpc'))
p = yt.PhasePlot(sp2, ('gas', 'density'), ('gas', 'temperature'), ('gas', 'mass'), weight_field=None)
p.set_unit(('gas', 'mass'), 'Msun')
p.save()
