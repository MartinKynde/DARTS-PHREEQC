import matplotlib.pyplot as plt

from model import Model
from darts.engines import redirect_darts_output, value_vector
import numpy as np
from darts.tools.plot_darts import *

redirect_darts_output('log.txt')


m = Model()

# Initialize simulations
m.init()

_, poro, _ = m.evaluate_porosity()
volume = np.array(m.reservoir.mesh.volume, copy=False)
total_pv = np.sum(volume[:m.reservoir.nb]*poro)
av_perm = m.const_perm
print("----------------------------------------------")
print("Total pore volume = %f m3" % total_pv)
print("Injection rate = %f m3/hour" % (m.inj_rate / 24))
print("Average initial permeability = %f mDarcy" % av_perm)
print("----------------------------------------------")
year = 365
NT = 10
runtime = 1000 / NT
m.export_vtk()
for i in range(NT):
    m.run_python(runtime)
    m.export_vtk()

print('\nNegative composition occurrence:', m.physics.acc_flux_etor.counter, '\n')

# Plot porosity
# ----------------------------------------------------------------------------------------------------------------------
poro_init, poro, poro_diff = m.evaluate_porosity()
k = 1

#m.plot_layer_map(poro_init, k, 'Porosity initial')
# m.plot_layer_surface(poro_init, k, 'Porosity initial')
m.plot_layer_map(poro, k, 'Porosity')
# m.plot_layer_surface(poro, k, 'Porosity')

#m.plot_layer_map(poro_diff, k, 'Porosity difference')
# m.plot_layer_surface(poro_diff, k, 'Porosity difference')
# ----------------------------------------------------------------------------------------------------------------------

# Plot other parameters
# ----------------------------------------------------------------------------------------------------------------------
Xm = np.copy(m.physics.engine.X[:m.reservoir.nb*m.physics.n_components])

map_data = Xm[1::m.physics.n_components]
name = 'Solid composition'
m.plot_layer_map(map_data, k, name)

map_data = Xm[0::m.physics.n_components]
name = 'Pressure'
m.plot_layer_map(map_data, k, name)

map_data = Xm[2::m.physics.n_components]
name = 'Ba-composition'
m.plot_layer_map(map_data, k, name)

map_data = Xm[3::m.physics.n_components]
name = 'S-composition'
m.plot_layer_map(map_data, k, name)

map_data = Xm[4::m.physics.n_components]
name = 'Ca-composition'
m.plot_layer_map(map_data, k, name)

map_data = Xm[5::m.physics.n_components]
name = 'O-composition'
m.plot_layer_map(map_data, k, name)

map_data = 1 - (Xm[1::m.physics.n_components] + Xm[2::m.physics.n_components]
                 + Xm[3::m.physics.n_components] + Xm[4::m.physics.n_components] +Xm[5::m.physics.n_components])
name = 'H composition'
m.plot_layer_map(map_data, k, name)
# ----------------------------------------------------------------------------------------------------------------------

#m.plot_geothermal_temp_layer_map(map_data,k,"GeothermalTempLayer")
#m.plot_1d(map_data,k,"GeothermalTempLayer")
time_data = pd.DataFrame.from_dict(m.physics.engine.time_data)



# Print statistics
m.print_timers()
m.print_stat()


