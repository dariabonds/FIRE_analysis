import yt
import trident
import unyt
import matplotlib.pyplot as plt
from yt.utilities.math_utils import ortho_find
from math import pi

##set particle filters
def angle_I(pfilter, data):
    filter = (data[(pfilter.filtered_type, 'phi')] > 0) \
             & (data[(pfilter.filtered_type, 'phi')] < 15)
    return filter

yt.add_particle_filter('angle_I', function=angle_I, filtered_type='gas', requires=['phi'])

def angle_II(pfilter, data):
    filter = (data[(pfilter.filtered_type, 'phi')] > 15) \
            & (data[(pfilter.filtered_type, 'phi')] < 30)
    return filter

yt.add_particle_filter('angle_II', function=angle_II, filtered_type='gas', requires=['phi'])

def angle_III(pfilter, data):
    filter = (data[(pfilter.filtered_type, 'phi')] > 30) \
            & (data[(pfilter.filtered_type, 'phi')] < 45)
    return filter

yt.add_particle_filter('angle_III', function=angle_III, filtered_type='gas', requires=['phi'])

def angle_IV(pfilter, data):
    filter = (data[(pfilter.filtered_type, 'phi')] > 45) \
            & (data[(pfilter.filtered_type, 'phi')] < 60)
    return filter

yt.add_particle_filter('angle_IV', function=angle_IV, filtered_type='gas', requires=['phi'])

def angle_V(pfilter, data):
    filter = (data[(pfilter.filtered_type, 'phi')] > 60) \
            & (data[(pfilter.filtered_type, 'phi')] < 75)
    return filter

yt.add_particle_filter('angle_V', function=angle_V, filtered_type='gas', requires=['phi'])

def angle_VI(pfilter, data):
    filter = (data[(pfilter.filtered_type, 'phi')] >= 75) \
            & (data[(pfilter.filtered_type, 'phi')] <= 90)
    return filter

yt.add_particle_filter('angle_VI', function=angle_VI, filtered_type='gas', requires=['phi'])

##define angle phi
##phi: radians to degrees
def _phi(field, data):
    return ((abs(data['gas', 'spherical_position_theta'] - (pi/2.0))) - (pi/2.0)) \
            * (-1) * (180/pi) * unyt.deg

yt.add_field(name=('gas', 'phi'), function=_phi, sampling_type='local', units='deg')

##load in dataset
#ds = yt.load('/Users/chummels/scratch/FIRE/m12i_res57000/output/snapshot_570.hdf5')
ds = yt.load('../m12i_res56000_md/snapshot_600.hdf5')
#ds = yt.load('../m12i_res7100_md/snapdir_600/snapshot_600.0.hdf5')

##find galaxy center
_, c = ds.find_max(('gas', 'density'))
##specific coord (x,y,z)
#c = ds.arr([29338.09863660, 30980.12414340, 32479.90455557], "code_length")

##add filters to data
ds.add_particle_filter('angle_I')
ds.add_particle_filter('angle_II')
ds.add_particle_filter('angle_III')
#ds.add_particle_filter('angle_IV')
#ds.add_particle_filter('angle_V')
#ds.add_particle_filter('angle_VI')

##plot galaxy face on
#p = yt.ProjectionPlot(ds, 'x', ('gas', 'density'), center=c, width=(100, 'kpc'))
#p = yt.ProjectionPlot(ds, 'x', ('gas', 'density'), center=c, width=(10, 'kpc'))
#p.save()

##create sphere
sp1 = ds.sphere(c, (10.0, 'kpc'))
##find max angular momentum
L = sp1.quantities.angular_momentum_vector()
##find orthogonal vectors
L, edge1, edge2 = ortho_find(L)

##plot galaxy edge on
#p = yt.OffAxisProjectionPlot(ds, edge2, ('gas', 'density'), center=c, width=(100, 'kpc'), north_vector=L)
#p = yt.OffAxisProjectionPlot(ds, edge2, ('gas', 'density'), center=c, width=(10, 'kpc'), north_vector=L)
#p.save()

##create a phase plot
sp2 = ds.sphere(c, (200.0, 'kpc'))
#p = yt.PhasePlot(sp2, ('gas', 'density'), ('gas', 'temperature'), ('gas', 'mass'), weight_field=None)
#p.set_unit(('gas', 'mass'), 'Msun')
#p.save()

##set sphere center and angular momentum
sp2.set_field_parameter('center', c)
sp2.set_field_parameter('normal', L)

##create particle phase plot
#p = yt.PhasePlot(sp2, ('gas', 'spherical_position_radius'), ('gas', 'phi'), \
#   ('gas', 'mass'), weight_field=None)
#p.set_unit(('gas', 'mass'), 'Msun')
#p.set_unit(('gas', 'spherical_position_radius'), 'kpc')
#p.set_unit(('gas', 'phi'), 'deg')
#p.set_xlim(1e-2, 1e3)
#p.set_ylim(1,90)
#p.set_log(('gas', 'phi'), False)
#p.save()

##apply angle filters
#p1 = yt.PhasePlot(sp2, ('angle_I', 'spherical_position_radius'), \
#    ('angle_I', 'phi'), ('angle_I', 'mass'), weight_field=None)
#p1.set_unit(('angle_I', 'mass'), 'Msun')
#p1.set_unit(('angle_I', 'spherical_position_radius'), 'kpc')
#p1.set_log(('angle_I', 'phi'), False)
#p1.set_xlim(1e-2, 150)
#p1.set_ylim(0, 90)
#p1.set_zlim(('angle_I', 'mass'), 5e4, 1e9)
#p1.save()

#p2 = yt.PhasePlot(sp2, ('angle_II', 'spherical_position_radius'), \
#    ('angle_II', 'phi'), ('angle_II', 'mass'), weight_field=None)
#p2.set_unit(('angle_II', 'mass'), 'Msun')
#p2.set_unit(('angle_II', 'spherical_position_radius'), 'kpc')
#p2.set_log(('angle_II', 'phi'), False)
#p2.set_xlim(1e-2, 150)
#p2.set_ylim(0, 90)
#p2.set_zlim(('angle_II', 'mass'), 5e4, 1e9)
#p2.save()

#p3 = yt.PhasePlot(sp2, ('angle_III', 'spherical_position_radius'), \
#    ('angle_III', 'phi'), ('angle_III', 'mass'), weight_field=None)
#p3.set_unit(('angle_III', 'mass'), 'Msun')
#p3.set_unit(('angle_III', 'spherical_position_radius'), 'kpc')
#p3.set_log(('angle_III', 'phi'), False)
#p3.set_xlim(1e-2, 150)
#p3.set_ylim(0, 90)
#p3.set_zlim(('angle_III', 'mass'), 5e4, 1e9)
#p3.save()

#p4 = yt.PhasePlot(sp2, ('angle_IV', 'particle_position_spherical_radius'), \
#    ('angle_IV', 'phi'), ('angle_IV', 'mass'), weight_field=None)
#p4.set_unit(('angle_IV', 'mass'), 'Msun')
#p4.set_unit(('angle_IV', 'particle_position_spherical_radius'), 'kpc')
#p4.set_log(('angle_IV', 'phi'), False)
#p4.set_xlim(1e-2, 150)
#p4.set_ylim(0, 90)
#p4.set_zlim(('angle_IV', 'mass'), 5e4, 1e9)
#p4.save()

#p5 = yt.PhasePlot(sp2, ('angle_V', 'spherical_position_radius'), \
#    ('angle_V', 'phi'), ('angle_V', 'mass'), weight_field=None)
#p5.set_unit(('angle_V', 'mass'), 'Msun')
#p5.set_unit(('angle_V', 'spherical_position_radius'), 'kpc')
#p5.set_log(('angle_V', 'phi'), False)
#p5.set_xlim(1e-2, 150)
#p5.set_ylim(0, 90)
#p5.set_zlim(('angle_V', 'mass'), 5e4, 1e9)
#p5.save()

#p6 = yt.PhasePlot(sp2, ('angle_VI', 'spherical_position_radius'), \
#    ('angle_VI', 'phi'), ('angle_VI', 'mass'), weight_field=None)
#p6.set_unit(('angle_VI', 'mass'), 'Msun')
#p6.set_unit(('angle_VI', 'spherical_position_radius'), 'kpc')
#p6.set_log(('angle_VI', 'phi'), False)
#p6.set_xlim(1e-2, 150)
#p6.set_ylim(0, 90)
#p6.set_zlim(('angle_VI', 'mass'), 5e4, 1e9)
#p6.save()

##compute bulk velocity
bulk_vel = sp2.quantities.bulk_velocity()
##set new sphere
sp3 = ds.sphere(c, (200.0, 'kpc'))
##set bulk velocity field parameter
sp3.set_field_parameter('bulk_velocity', bulk_vel)

##radial velocitys by angle bins
##(‘angle_’, ‘radial_velocity’)
##(‘angle_’, ‘particle_radial_velocity’)
##(‘angle_’, ‘density’)
##(‘angle_’, ‘Density’)
##similar to p.set_y_units(‘Msun’) ‘g cm**-3’ ‘g/cm**3’

#rp0 = yt.create_profile(sp3, ('gas', 'spherical_position_radius'), ('gas', 'radial_velocity'), \
#    units={('gas', 'spherical_position_radius'): 'kpc'}, logs={('gas', 'spherical_position_radius'): False})

rp1 = yt.create_profile(sp3, ('angle_I', 'spherical_position_radius'), ('angle_I', 'radial_velocity'), weight_field=('angle_I', 'mass'), \
    units={('angle_I', 'spherical_position_radius'): 'kpc'}, logs={('angle_I', 'spherical_position_radius'): False})

rp2 = yt.create_profile(sp3, ('angle_II', 'spherical_position_radius'), ('angle_II', 'radial_velocity'), weight_field=('angle_II', 'mass'), \
     units={('angle_II', 'spherical_position_radius'): 'kpc'}, logs={('angle_II', 'spherical_position_radius'): False})

rp3 = yt.create_profile(sp3, ('angle_III', 'spherical_position_radius'), ('angle_III', 'radial_velocity'), weight_field=('angle_III', 'mass'), \
     units={('angle_III', 'spherical_position_radius'): 'kpc'}, logs={('angle_III', 'spherical_position_radius'): False})

##radial velocity profile
p = plt.figure()
ax = p.add_subplot(111)
#ax.plot(rp1.x.value, rp1[("gas", "radial_velocity")].in_units("km/s").value)
ax.plot(rp1.x.value, rp1[("angle_I", "radial_velocity")].in_units("km/s").value, \
        rp2.x.value, rp2[("angle_II", "radial_velocity")].in_units("km/s").value, \
        rp3.x.value, rp3[("angle_III", "radial_velocity")].in_units("km/s").value)
ax.set_xlabel(r"$\mathrm{r\ (kpc)}$")
ax.set_ylabel(r"$\mathrm{v_r\ (km/s)}$")
#ax.legend(["With Correction"])
ax.legend(["0-15", "15-30", "30-45"])
p.savefig("snapshot_600_radial_velocity_profile_0.png")
#p.savefig("snapshot_600_radial_velocity_profile_I.png")

import sys; sys.exit()

##generate species plots
trident.add_ion_fields(ds, ions=['O VI'], ftype='gas')
trident.add_ion_fields(ds, ions=['Mg II'], ftype='gas')

pO = yt.ProjectionPlot(ds, 'z', 'O_p5_number_density')
pO.save()

pMg = yt.ProjectionPlot(ds, 'z', 'Mg_p1_number_density')
pMg.save()

adO = ds.all_data()
phaseO = yt.PhasePlot(adO, ('gas', 'density'), ('gas', 'temperature'), ('gas', 'O_p5_mass'), weight_field=None, fractional=True)
phaseO.save()

adMg = ds.all_data()
phaseMg = yt.PhasePlot(adMg, ('gas', 'density'), ('gas', 'temperature'), ('gas', 'Mg_p1_mass'), weight_field=None, fractional=True)
phaseMg.save()
