import yt

ds = yt.load('../m12i_res56000_md/snapshot_600.hdf5')
_, c = ds.find_max(('gas', 'density'))
#p = yt.ProjectionPlot(ds, 'x', ('gas', 'density'), center=c, width=(100, 'kpc'))
#p.save()

sp = ds.sphere(c, (25.0, 'kpc'))
L = sp.quantities.angular_momentum_vector()
p = yt.OffAxisProjectionPlot(ds, L, ('gas', 'density'), center=c, width=(100, 'kpc'))
p.save()
