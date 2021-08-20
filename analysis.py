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
#ds = yt.load('../m12i_res56000_md/snapshot_600.hdf5') #ds1
#ds = yt.load('../m12i_res7100_md/snapdir_600/snapshot_600.0.hdf5') #ds2
#ds = yt.load('/mnt/data1/GalaxiesOnFire/metaldiff/m12i_res7100_md/output/snapdir_600/snapshot_600.0.hdf5') #ds2
ds = yt.load('/mnt/data1/GalaxiesOnFire/cr_700/output/snapdir_600/snapshot_600.0.hdf5') #ds3

##find galaxy center
#_, c = ds.find_max(('gas', 'density'))
##specific coord (x,y,z)
#c = ds.arr([29338.09863660, 30980.12414340, 32479.90455557], 'code_length') #ds2
c = ds.arr([29345.27830223, 30997.08859958, 32484.0642261], 'code_length') #ds3

##add filters to data
#ds.add_particle_filter('angle_I')
#ds.add_particle_filter('angle_II')
#ds.add_particle_filter('angle_III')
#ds.add_particle_filter('angle_IV')
#ds.add_particle_filter('angle_V')
#ds.add_particle_filter('angle_VI')

##plot galaxy face on
#p = yt.ProjectionPlot(ds, 'x', ('gas', 'density'), center=c, width=(100, 'kpc'))
#p.save()

##create sphere
sp1 = ds.sphere(c, (10.0, 'kpc'))
##find max angular momentum
L = sp1.quantities.angular_momentum_vector()
##find orthogonal vectors
L, edge1, edge2 = ortho_find(L)

##plot galaxy edge on
#p = yt.OffAxisProjectionPlot(ds, edge2, ('gas', 'density'), center=c, width=(200, 'kpc'), north_vector=L)
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

#p4 = yt.PhasePlot(sp2, ('angle_IV', 'spherical_position_radius'), \
#    ('angle_IV', 'phi'), ('angle_IV', 'mass'), weight_field=None)
#p4.set_unit(('angle_IV', 'mass'), 'Msun')
#p4.set_unit(('angle_IV', 'spherical_position_radius'), 'kpc')
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

##radial velocities by angle bins
##(‘angle_’, ‘radial_velocity’)
#rp1 = yt.create_profile(sp3, ('angle_I', 'spherical_position_radius'), ('angle_I', 'radial_velocity'), weight_field=('angle_I', 'mass'), \
#    units={('angle_I', 'spherical_position_radius'): 'kpc'}, logs={('angle_I', 'spherical_position_radius'): False})

#rp2 = yt.create_profile(sp3, ('angle_II', 'spherical_position_radius'), ('angle_II', 'radial_velocity'), weight_field=('angle_II', 'mass'), \
#    units={('angle_II', 'spherical_position_radius'): 'kpc'}, logs={('angle_II', 'spherical_position_radius'): False})

#rp3 = yt.create_profile(sp3, ('angle_III', 'spherical_position_radius'), ('angle_III', 'radial_velocity'), weight_field=('angle_III', 'mass'), \
#    units={('angle_III', 'spherical_position_radius'): 'kpc'}, logs={('angle_III', 'spherical_position_radius'): False})

#rp4 = yt.create_profile(sp3, ('angle_IV', 'spherical_position_radius'), ('angle_IV', 'radial_velocity'), weight_field=('angle_IV', 'mass'), \
#    units={('angle_IV', 'spherical_position_radius'): 'kpc'}, logs={('angle_IV', 'spherical_position_radius'): False})

#rp5 = yt.create_profile(sp3, ('angle_V', 'spherical_position_radius'), ('angle_V', 'radial_velocity'), weight_field=('angle_V', 'mass'), \
#    units={('angle_V', 'spherical_position_radius'): 'kpc'}, logs={('angle_V', 'spherical_position_radius'): False})

#rp6 = yt.create_profile(sp3, ('angle_VI', 'spherical_position_radius'), ('angle_VI', 'radial_velocity'), weight_field=('angle_VI', 'mass'), \
#    units={('angle_VI', 'spherical_position_radius'): 'kpc'}, logs={('angle_VI', 'spherical_position_radius'): False})

##radial velocity profile
#p = plt.figure()
#ax = p.add_subplot(111)
#ax.plot(rp1.x.value, rp1[("angle_I", "radial_velocity")].in_units("km/s").value, \
#        rp2.x.value, rp2[("angle_II", "radial_velocity")].in_units("km/s").value, \
#        rp3.x.value, rp3[("angle_III", "radial_velocity")].in_units("km/s").value, \
#        rp4.x.value, rp4[("angle_IV", "radial_velocity")].in_units("km/s").value, \
#        rp5.x.value, rp5[("angle_V", "radial_velocity")].in_units("km/s").value, \
#        rp6.x.value, rp6[("angle_VI", "radial_velocity")].in_units("km/s").value)
#ax.set_xlabel(r"$\mathrm{r\ (kpc)}$")
#ax.set_ylabel(r"$\mathrm{v_r\ (km/s)}$")
#ax.legend(["0-15", "15-30", "30-45", "45-60", "60-75", "75-90"])
#ax.set_ylim(-100, 100)
#p.savefig("snapshot_600_radial_velocity_profile_1.1.png")
#p.savefig("snapshot_600_radial_velocity_profile_2.png")
#p.savefig("snapshot_600_radial_velocity_profile_CR_2.png")

##radial densities by angle bins
##(‘angle_’, ‘density’)
#rp1 = yt.create_profile(sp3, ('angle_I', 'spherical_position_radius'), ('angle_I', 'density'), weight_field=('angle_I', 'mass'), \
#    units={('angle_I', 'spherical_position_radius'): 'kpc'}, logs={('angle_I', 'spherical_position_radius'): True})

#rp2 = yt.create_profile(sp3, ('angle_II', 'spherical_position_radius'), ('angle_II', 'density'), weight_field=('angle_II', 'mass'), \
#    units={('angle_II', 'spherical_position_radius'): 'kpc'}, logs={('angle_II', 'spherical_position_radius'): True})

#rp3 = yt.create_profile(sp3, ('angle_III', 'spherical_position_radius'), ('angle_III', 'density'), weight_field=('angle_III', 'mass'), \
#    units={('angle_III', 'spherical_position_radius'): 'kpc'}, logs={('angle_III', 'spherical_position_radius'): True})

#rp4 = yt.create_profile(sp3, ('angle_IV', 'spherical_position_radius'), ('angle_IV', 'density'), weight_field=('angle_IV', 'mass'), \
#    units={('angle_IV', 'spherical_position_radius'): 'kpc'}, logs={('angle_IV', 'spherical_position_radius'): True})

#rp5 = yt.create_profile(sp3, ('angle_V', 'spherical_position_radius'), ('angle_V', 'density'), weight_field=('angle_V', 'mass'), \
#    units={('angle_V', 'spherical_position_radius'): 'kpc'}, logs={('angle_V', 'spherical_position_radius'): True})

#rp6 = yt.create_profile(sp3, ('angle_VI', 'spherical_position_radius'), ('angle_VI', 'density'), weight_field=('angle_VI', 'mass'), \
#    units={('angle_VI', 'spherical_position_radius'): 'kpc'}, logs={('angle_VI', 'spherical_position_radius'): True})

##radial density profile
#p = plt.figure()
#ax = p.add_subplot(111)
#ax.plot(rp1.x.value, rp1[("angle_I", "density")].in_units("g/cm**3").value, \
#        rp2.x.value, rp2[("angle_II", "density")].in_units("g/cm**3").value, \
#        rp3.x.value, rp3[("angle_III", "density")].in_units("g/cm**3").value, \
#        rp4.x.value, rp4[("angle_IV", "density")].in_units("g/cm**3").value, \
#        rp5.x.value, rp5[("angle_V", "density")].in_units("g/cm**3").value, \
#        rp6.x.value, rp6[("angle_VI", "density")].in_units("g/cm**3").value)
#ax.set_yscale('log')
#ax.set_xlabel(r"$\mathrm{r\ (kpc)}$")
#ax.set_ylabel(r"$\mathrm{rho\ (g/cm**3)}$")
#ax.legend(["0-15", "15-30", "30-45", "45-60", "60-75", "75-90"])
#p.savefig("snapshot_600_radial_density_profile_1.png")
#p.savefig("snapshot_600_radial_density_profile_2.png")
#p.savefig("snapshot_600_radial_density_profile_3.png")
#p.savefig("snapshot_600_radial_density_profile_CR.png")

##generate ion species felds
trident.add_ion_fields(ds, ions=['O VI'], ftype='gas')
trident.add_ion_fields(ds, ions=['Mg II'], ftype='gas')
trident.add_ion_fields(ds, ions=['H I'], ftype='gas')

##face on ion field projection plots
#pO = yt.ProjectionPlot(ds, 'z', 'O_p5_number_density', center=c, width=(500, 'kpc'))
#pO.save()

#pMg = yt.ProjectionPlot(ds, 'z', 'Mg_p1_number_density', center=c, width=(500, 'kpc'))
#pMg.save()

#pH = yt.ProjectionPlot(ds, 'z', 'H_p0_number_density', center=c, width=(500, 'kpc'))
#pH.save()

##off axis ion projection plots
pO2 = yt.OffAxisProjectionPlot(ds, edge2, 'O_p5_number_density', center=c, width=(100, 'kpc'), north_vector=L)
pO2.save()

pMg2 = yt.OffAxisProjectionPlot(ds, edge2, 'Mg_p1_number_density', center=c, width=(100, 'kpc'), north_vector=L)
pMg2.save()

pH2 = yt.OffAxisProjectionPlot(ds, edge2, 'H_p0_number_density', center=c, width=(100, 'kpc'), north_vector=L)
pH2.save()

##ion field phase plots
#adO = ds.all_data()
#phaseO = yt.PhasePlot(adO, ('gas', 'density'), ('gas', 'temperature'), ('gas', 'O_p5_mass'), weight_field=None, fractional=False)
#phaseO.save()

#adMg = ds.all_data()
#phaseMg = yt.PhasePlot(adMg, ('gas', 'density'), ('gas', 'temperature'), ('gas', 'Mg_p1_mass'), weight_field=None, fractional=False)
#phaseMg.save()

#adH = ds.all_data()
#phaseH = yt.PhasePlot(adH, ('gas', 'density'), ('gas', 'temperature'), ('gas', 'H_p0_mass'), weight_field=None, fractional=False)
#phaseH.save()

##off axis slice and velocity vectors
#p = yt.SlicePlot(ds, 'z', ('gas', 'density'), center=c, width=(200, 'kpc'), data_source=sp3)
#p.annotate_velocity(plot_args={"headwidth": 5})
#p.save()
