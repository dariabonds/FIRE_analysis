import yt
import trident
import unyt
import matplotlib.pyplot as plt
from yt.utilities.math_utils import ortho_find
from math import pi

##set particle filters #consider defining deg in function
def angle_low(pfilter, data):
    #filter = data[(pfilter.filtered_type, 'particle_position_spherical_theta')] > 0 and
    filter = data[(pfilter.filtered_type, 'phi')] < 15
    return filter

yt.add_particle_filter('angle_low', function=angle_low, filtered_type='PartType0', requires=['phi'])

def angle_low_mid(pfilter, data):
    filter = (data[(pfilter.filtered_type, 'phi')] > 15) \
            & (data[(pfilter.filtered_type, 'phi')] < 30)
    return filter

yt.add_particle_filter('angle_low_mid', function=angle_low_mid, filtered_type='PartType0', requires=['phi'])

def angle_mid_one(pfilter, data):
    filter = (data[(pfilter.filtered_type, 'phi')] > 30) \
            & (data[(pfilter.filtered_type, 'phi')] < 45)
    return filter

yt.add_particle_filter('angle_mid_one', function=angle_mid_one, filtered_type='PartType0', requires=['phi'])

def angle_mid_two(pfilter, data):
    filter = (data[(pfilter.filtered_type, 'phi')] > 45) \
            & (data[(pfilter.filtered_type, 'phi')] < 60)
    return filter

yt.add_particle_filter('angle_mid_two', function=angle_mid_two, filtered_type='PartType0', requires=['phi'])

def angle_mid_high(pfilter, data):
    filter = (data[(pfilter.filtered_type, 'phi')] > 60) \
            & (data[(pfilter.filtered_type, 'phi')] < 75)
    return filter

yt.add_particle_filter('angle_mid_high', function=angle_mid_high, filtered_type='PartType0', requires=['phi'])

def angle_high(pfilter, data):
    filter = (data[(pfilter.filtered_type, 'phi')] >= 75) \
            & (data[(pfilter.filtered_type, 'phi')] <= 90)
    return filter

yt.add_particle_filter('angle_high', function=angle_high, filtered_type='PartType0', requires=['phi'])

##define angle phi
##phi: radians to degrees
def _phi(field, data):
    return ((abs(data['PartType0', 'particle_position_spherical_theta'] - (pi/2.0))) - (pi/2.0)) \
            * (-1) * (180/pi) * unyt.deg

yt.add_field(name=('PartType0', 'phi'), function=_phi, sampling_type='particle', units='deg')

##load in dataset
ds = yt.load('../m12i_res56000_md/snapshot_600.hdf5')
#ds = yt.load('../m12i_res7100_md/snapdir_600/snapshot_600.0.hdf5')

##add filters to data
ds.add_particle_filter('angle_low')
ds.add_particle_filter('angle_low_mid')
ds.add_particle_filter('angle_mid_one')
ds.add_particle_filter('angle_mid_two')
ds.add_particle_filter('angle_mid_high')
ds.add_particle_filter('angle_high')

##find galaxy center
_, c = ds.find_max(('gas', 'density'))
##specific coord (x,y,z)
#c = [29338.09863660, 30980.12414340, 32479.90455557]

##plot galaxy face on
#p = yt.ProjectionPlot(ds, 'x', ('gas', 'density'), center=c, width=(100, 'kpc'))
#p = yt.ProjectionPlot(ds, 'x', ('gas', 'density'), center=c, width=(10, 'Mpc'))
#p.save()

##create sphere
#sp = ds.sphere(c, (20.0, 'kpc'))
sp1 = ds.sphere(c, (10.0, 'kpc'))
##find max angular momentum
L = sp1.quantities.angular_momentum_vector()
##find orthogonal vectors
L, edge1, edge2 = ortho_find(L)

##plot galaxy edge on
#p = yt.OffAxisProjectionPlot(ds, edge2, ('gas', 'density'), center=c, width=(100, 'kpc'), north_vector=L)
#p = yt.OffAxisProjectionPlot(ds, edge2, ('gas', 'density'), center=c, width=(50, 'kpc'), north_vector=L)
#p.save()

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
#p = yt.PhasePlot(sp2, ('PartType0', 'particle_position_spherical_radius'), ('gas', 'phi'), \
#   ('gas', 'mass'), weight_field=None)
#p.set_unit(('gas', 'mass'), 'Msun')
#p.set_unit(('PartType0', 'particle_position_spherical_radius'), 'kpc')
#p.set_unit(('gas', 'phi'), 'deg')
#p.set_xlim(1e-2, 1e3)
#p.set_ylim(0,90)
#p.set_log(('gas', 'phi'), False)
#p.save()

##apply angle filters
#p1 = yt.PhasePlot(sp2, ('angle_low', 'particle_position_spherical_radius'), \
#    ('angle_low', 'phi'), ('angle_low', 'mass'), weight_field=None)
#p1.set_unit(('angle_low', 'mass'), 'Msun')
#p1.set_unit(('angle_low', 'particle_position_spherical_radius'), 'kpc')
#p1.set_log(('angle_low', 'phi'), False)
#p1.set_xlim(1e-2, 150)
#p1.set_ylim(0, 90)
#p1.set_zlim(('angle_low', 'mass'), 5e4, 1e9)
#p1.save()

#p2 = yt.PhasePlot(sp2, ('angle_low_mid', 'particle_position_spherical_radius'), \
#    ('angle_low_mid', 'phi'), ('angle_low_mid', 'mass'), weight_field=None)
#p2.set_unit(('angle_low_mid', 'mass'), 'Msun')
#p2.set_unit(('angle_low_mid', 'particle_position_spherical_radius'), 'kpc')
#p2.set_log(('angle_low_mid', 'phi'), False)
#p2.set_xlim(1e-2, 150)
#p2.set_ylim(0, 90)
#p2.set_zlim(('angle_low_mid', 'mass'), 5e4, 1e9)
#p2.save()

#p3 = yt.PhasePlot(sp2, ('angle_mid_one', 'particle_position_spherical_radius'), \
#    ('angle_mid_one', 'phi'), ('angle_mid_one', 'mass'), weight_field=None)
#p3.set_unit(('angle_mid_one', 'mass'), 'Msun')
#p3.set_unit(('angle_mid_one', 'particle_position_spherical_radius'), 'kpc')
#p3.set_log(('angle_mid_one', 'phi'), False)
#p3.set_xlim(1e-2, 150)
#p3.set_ylim(0, 90)
#p3.set_zlim(('angle_mid_one', 'mass'), 5e4, 1e9)
#p3.save()

#p4 = yt.PhasePlot(sp2, ('angle_mid_two', 'particle_position_spherical_radius'), \
#    ('angle_mid_two', 'phi'), ('angle_mid_two', 'mass'), weight_field=None)
#p4.set_unit(('angle_mid_two', 'mass'), 'Msun')
#p4.set_unit(('angle_mid_two', 'particle_position_spherical_radius'), 'kpc')
#p4.set_log(('angle_mid_two', 'phi'), False)
#p4.set_xlim(1e-2, 150)
#p4.set_ylim(0, 90)
#p4.set_zlim(('angle_mid_two', 'mass'), 5e4, 1e9)
#p4.save()

#p5 = yt.PhasePlot(sp2, ('angle_mid_high', 'particle_position_spherical_radius'), \
#    ('angle_mid_high', 'phi'), ('angle_mid_high', 'mass'), weight_field=None)
#p5.set_unit(('angle_mid_high', 'mass'), 'Msun')
#p5.set_unit(('angle_mid_high', 'particle_position_spherical_radius'), 'kpc')
#p5.set_log(('angle_mid_high', 'phi'), False)
#p5.set_xlim(1e-2, 150)
#p5.set_ylim(0, 90)
#p5.set_zlim(('angle_mid_high', 'mass'), 5e4, 1e9)
#p5.save()

#p6 = yt.PhasePlot(sp2, ('angle_high', 'particle_position_spherical_radius'), \
#    ('angle_high', 'phi'), ('angle_high', 'mass'), weight_field=None)
#p6.set_unit(('angle_high', 'mass'), 'Msun')
#p6.set_unit(('angle_high', 'particle_position_spherical_radius'), 'kpc')
#p6.set_log(('angle_high', 'phi'), False)
#p6.set_xlim(1e-2, 150)
#p6.set_ylim(0, 90)
#p6.set_zlim(('angle_high', 'mass'), 5e4, 1e9)
#p6.save()
