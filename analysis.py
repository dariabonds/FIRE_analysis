import yt
from yt.utilities.math_utils import ortho_find
from math import pi

##set particle filters
def angle_low(pfilter, data):
    #filter = data[(pfilter.filtered_type, 'particle_position_spherical_theta')] > 0 and
    filter = data[(pfilter.filtered_type, 'particle_position_spherical_theta')] < (pi/6)
    return filter

yt.add_particle_filter('angle_low', function=angle_low, filtered_type='PartType0', requires=['particle_position_spherical_theta'])

#def angle_mid(pfilter, data):
#    filter = (data[(pfilter.filtered_type, 'particle_position_spherical_theta')], (>= (pi/6) and <= (pi/4)))
#    return filter

#yt.particle_filter('angle_mid', function=angle_mid, filtered_type='PartType0', requires=['particle_position_spherical_theta'])

#def angle_high(pfilter, data):
#    filter = (data[(pfilter.filtered_type, 'particle_position_spherical_theta')], (>= (pi/4) and <= (pi/2)))
#    return filter

#yt.particle_filter('angle_high', function=angles, filtered_type='PartType0', requires=['particle_position_spherical_theta'])

##define
def _phi(field, data):
    return ((abs(data['PartType0', 'particle_position_spherical_theta']-(pi/2.0)))-(pi/2.0))*(-1) * ds.units.rad

yt.add_field(name=('gas', 'phi'), function=_phi, sampling_type='local', units='rad')

##load in dataset
ds = yt.load('../m12i_res56000_md/snapshot_600.hdf5')
#ds = yt.load('../m12i_res7100_md/snapdir_600/snapshot_600.0.hdf5')

ds.add_particle_filter('angle_low')

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
#p = yt.PhasePlot(sp2, ('PartType0', 'particle_position_spherical_radius'), ('gas', 'phi'), ('gas', 'mass'), weight_field=None)
#p.set_unit(('gas', 'mass'), 'Msun')
#p.set_unit(('PartType0', 'particle_position_spherical_radius'), 'kpc')
#p.set_unit(('gas', 'phi'), 'deg')
#p.set_xlim(1e-2, 1e3)
#p.set_ylim(0,90)
#p.set_log(('gas', 'phi'), False)
#p.save()

##apply angle filters
pl = yt.PhasePlot(sp2, ('angle_low', 'particle_position_spherical_radius'), ('angle_low', 'particle_position_spherical_theta'), ('angle_low', 'mass'), weight_field=None)
pl.set_unit(('angle_low', 'mass'), 'Msun')
pl.set_unit(('angle_low', 'particle_position_spherical_radius'), 'kpc')
#pl.set_log(('angle_low', 'particle_position_spherical_theta'), False)
pl.set_xlim(1e-2, 150)
pl.set_ylim(1e-2, pi)
pl.set_zlim(('angle_low', 'mass'), 5e4, 1e9)
#pl.set_unit(('gas', 'phi'), 'deg')
pl.save()

#pm = yt.PhasePlot(sp, ('angle_mid', 'particle_position_spherical_radius'), ('angle_mid', 'particle_position_spherical_theta'), ('gas', 'mass'), weight_field=None)
#pm.set_unit(('gas', 'mass'), 'Msun')
#pm.set_unit(('angle_mid', 'particle_position_spherical_radius'), 'kpc')
#pm.set_unit(('gas', 'phi'), 'deg')
#pm.save()

#ph = yt.PhasePlot(sp, ('angle_high', 'particle_position_spherical_radius'), ('angle_high', 'particle_position_spherical_theta'), ('gas', 'mass'), weight_field=None)
#ph.set_unit(('gas', 'mass'), 'Msun')
#ph.set_unit(('angle_high', 'particle_position_spherical_radius'), 'kpc')
#ph.set_unit(('gas', 'phi'), 'deg')
#ph.save()
