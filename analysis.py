import yt
from yt.utilities.math_utils import ortho_find

##load in dataset
ds = yt.load('../m12i_res56000_md/snapshot_600.hdf5')
#ds = yt.load('../m12i_res7100_md/snapdir_600/snapshot_600.0.hdf5')

##find galaxy center
_, c = ds.find_max(('gas', 'density'))
#particular coord (x,y,z)
#c = [29338.09863660, 30980.12414340, 32479.90455557]

##plot galaxy face on
#p = yt.ProjectionPlot(ds, 'x', ('gas', 'density'), center=c, width=(100, 'kpc'))
#p = yt.ProjectionPlot(ds, 'x', ('gas', 'density'), center=c, width=(10, 'Mpc'))
#p.save()
#import sys; sys.exit()

##create sphere
#sp = ds.sphere(c, (25.0, 'kpc'))
sp = ds.sphere(c, (10.0, 'kpc'))
##find max angular momentum
L = sp.quantities.angular_momentum_vector()
##find orthogonal vectors
L, edge1, edge2 = ortho_find(L)

##plot galaxy edge on
#p = yt.OffAxisProjectionPlot(ds, edge2, ('gas', 'density'), center=c, width=(100, 'kpc'), north_vector=L)
#p = yt.OffAxisProjectionPlot(ds, edge2, ('gas', 'density'), center=c, width=(50, 'kpc'), north_vector=L)
#p.save()
#import sys; sys.exit()

##create a phase plot
#sp2 = ds.sphere(c, (250.0, 'kpc'))
sp2 = ds.sphere(c, (150.0, 'kpc'))
#p = yt.PhasePlot(sp2, ('gas', 'density'), ('gas', 'temperature'), ('gas', 'mass'), weight_field=None)
#p.set_unit(('gas', 'mass'), 'Msun')
#p.save()

##set sphere center and angular momentum
sp2.set_field_parameter('center', c)
sp2.set_field_parameter('normal', L)

##create particle phase plot
p = yt.PhasePlot(sp2, ('PartType0', 'particle_position_spherical_radius'), ('PartType0', 'particle_position_spherical_phi'), ('gas', 'mass'), weight_field=None)
p.set_unit(('gas', 'mass'), 'Msun')
p.set_unit(('PartType0', 'particle_position_spherical_radius'), 'kpc')
# p.set_unit(('PartType0', 'particle_position_spherical_phi'), 'deg')
# p.set_log(('PartType0', 'particle_position_spherical_phi'), False)

p.save()
