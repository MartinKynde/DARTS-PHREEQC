from conversions import convert_composition, correct_composition, calculate_injection_stream, \
    get_mole_fractions, convert_rate, ml_min2ft3_d
from martin_physics import OwnPhysicsClass
from martin_properties import *
from darts.models.darts_model import DartsModel
from darts.models.reservoirs.struct_reservoir import StructReservoir
from darts.engines import *
from darts.physics import *

from win32com.client import Dispatch
import numpy as np
import pickle
import os
from math import fabs



# Definition of your input parameter data structure, change as you see fit (when you need more constant values, etc.)!!
class MyOwnDataStruct:
    def __init__(self, nc, zmin, temp, stoich_matrix, pressure_init, c_r, kin_fact,  exp_w=1, exp_g=1):
        """
        Data structure class which holds various input parameters for simulation
        :param nc: number of components used in simulation
        :param zmin: actual 0 used for composition (ussualy >0, around some small epsilon)
        :param temp: temperature
        """
        self.num_comp = nc
        self.min_z = zmin
        self.temperature = temp
        self.stoich_matrix = stoich_matrix
        self.exp_w = exp_w
        self.exp_g = exp_g
        self.pressure_init = pressure_init
        self.c_r = c_r
        self.kin_fact = kin_fact

# Actual Model class creation here!
class Model(DartsModel):
    def __init__(self):
        # Call base class constructor
        super().__init__()

        # Measure time spend on reading/initialization
        self.timer.node["initialization"].start()

        # Initialize PHREEQC
        self.phreeqc = Dispatch("IPhreeqcCOM.Object")

        # Static parameters input:
        self.nx = 50
        self.ny = 50
        self.nz = 50
        self.dx = 1  # m
        self.dy = 1  # m
        self.dz = 1  # m
        self.volume = self.dx * self.dy * self.dz
        self.depth = 1250  # m
        self.poro = 0.25  # [-]
        self.temperature = 323.15  # K
        self.pressure_init = 130  # bar

        #poro = np.genfromtxt('poro_LC02.txt')[1:self.nx * self.ny * self.nz + 1]
        #poro = 0.25 + np.random.uniform(-0.01, 0.01, self.nx * self.ny * self.nz)
        poro = 0.25
        self.solid_sat = 0

        trans_mult_exp = 4  # 6
        self.const_perm = 1600
        # self.const_perm = 5e6 * self.poro ** 6

        self.inj_rate = convert_rate(500)   # input: ml/min; output: m3/day

        self.c_r = 4.5e-5

        self.kin_fact = 1

        self.comp_min = 1e-7
        self.obl_min = self.comp_min / 10

        self.reservoir = StructReservoir(self.timer, nx=self.nx, ny=self.ny, nz=self.nz, dx=self.dx, dy=self.dy,
                                         dz=self.dz, permx=self.const_perm, permy=self.const_perm,
                                         permz=self.const_perm/10, poro=self.poro, depth=self.depth)

        # ================================================= Well index =================================================
        # ================================================= Well index =================================================
        d_w = 1.5
        r_w = d_w / 2
        # skin = 0
        # peaceman_rad = 0.28 * np.sqrt(self.dx ** 2 + self.dy ** 2) / (1 ** (1 / 4) + 1 ** (1 / 4))
        # kx = 5e6 * (1 - np.mean(self.solid_sat[0 * self.nx * self.ny::self.ny])) ** 6
        # ky = kx
        # well_index = 2 * np.pi * self.dz * np.sqrt(kx * ky) / (np.log(peaceman_rad / r_w) + skin)
        # well_index = well_index * 0.0085267146719160104986876640419948

        well_index = 5
        # ================================================= Well index =================================================


                # add wells
        self.reservoir.add_well("I1", wellbore_diameter=d_w)
        self.reservoir.add_perforation(well=self.reservoir.wells[-1], i=7, j=25, k=1,
                                               multi_segment=False,
                                               verbose=False, well_radius=r_w, well_index=well_index)

        self.reservoir.add_well("P1", wellbore_diameter=d_w)
        self.reservoir.add_perforation(well=self.reservoir.wells[-1], i=self.nx-7, j=25, k=1,
                                               multi_segment=False,
                                               verbose=False, well_radius=r_w, well_index=well_index)

        # Several parameters here related to components used, OBL limits, and injection composition:

        self.cell_property = ['pressure', 'H2O', 'H+', 'OH-', 'HSO4-', 'H3O+', 'SO4-2', 'BaSO4', 'Ba+2', 'BaOH+', 'Ca+2', 'CaOH+', 'CaSO4',
                              'Solid']
        self.phases = ['liq']
        self.components = ['H2O', 'H+', 'OH-', 'HSO4-', 'H3O+', 'SO4-2', 'BaSO4', 'Ba+2', 'BaOH+', 'Ca+2', 'CaOH+', 'CaSO4', 'Solid']
        self.elements = ['Solid', 'Ba', 'S', 'Ca', 'O', 'H']
        self.num_vars = len(self.elements)
        self.n_points = 101 #201
        self.min_p = 120
        self.max_p = 140
        self.min_z = self.obl_min
        self.max_z = 1 - self.obl_min

        # Rate annihilation matrix
        E = np.array([[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1],
                      [0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0],  # Ba
                      [0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 1, 0],  # S
                      [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0],  # Ca
                      [1, 0, 1, 4, 3, 4, 4, 0, 1, 0, 1, 4, 0],  # O
                      [2, 1, 1, 1, 3, 0, 0, 0, 1, 0, 1, 0, 0]])  # H

        # Several parameters related to kinetic reactions:
        stoich_matrix = np.array([1, -1, -1, 0, -4, 0])

        # Create instance of data-structure for simulation (and chemical) input parameters:
        input_data_struct = MyOwnDataStruct(len(self.elements), self.comp_min, self.temperature, stoich_matrix,
                                            self.pressure_init, self.c_r, self.kin_fact)

        # Create property containers:
        self.property_container = property_container(self.phases, self.components)

        self.property_container.phase_name = ['liq']

        self.property_container.rel_perm_ev = dict([('liq', custom_rel_perm(2))])

        self.property_container.flash_ev = custom_flash(self.temperature, self.phreeqc, self.comp_min)
        self.property_container.init_flash_ev = init_flash(self.temperature, self.phreeqc, self.comp_min)

        self.property_container.kin_rate_ev = custom_kinetic_rate(self.temperature, self.comp_min)

        # Create instance of (own) physics class:
        self.physics = OwnPhysicsClass(self.timer, self.elements, self.n_points, self.min_p, self.max_p,
                                       self.min_z, input_data_struct, self.property_container)

        # Compute injection stream
        ba_inj = 17.55/1000000  # molality from PHREEQC spectation17.55
        so4_inj = 86/1000000
        ca_inj = (6881+1492)/1000000
        water = 1 - ba_inj - so4_inj - ca_inj
        # Define injection stream composition,
        self.inj_stream_components = np.array([water, 0, 0, 0, 0, so4_inj, 0, ba_inj, 0, ca_inj, 0, 0, 0])

        self.inj_stream = convert_composition(self.inj_stream_components, E)
        self.inj_stream = correct_composition(self.inj_stream, self.comp_min)

        # ====================================== Initialize reservoir composition ======================================
        print('\tInitializing compositions...')

        bastart = 18/1000000
        so4start = 69/1000000
        castart = (6909+1564)/1000000
        waterstart = 1 - bastart - so4start - castart
        self.initial_comp_components = np.array([waterstart, 0, 0, 0, 0, so4start, 0, bastart, 0, castart, 0, 0, 0])
        self.solid_frac = []
        self.initial_comp = np.zeros((self.nx * self.ny * self.nz + 4, self.num_vars - 1))

        # Interpolated values of non-solid volume (second value always 0 due to no (5,1) interpolator)
        values = value_vector([0] * 2)

        # Iterate over solid saturation and call interpolator
        composition_full = convert_composition(self.initial_comp_components, E)
        composition = correct_composition(composition_full, self.comp_min)
        init_state = value_vector(np.hstack((self.solid_sat, composition)))

            # Call interpolator
        self.physics.comp_itor.evaluate(init_state, values)

            # Assemble initial composition
        self.solid_frac = values[0]
        initial_comp_with_solid = np.multiply(composition_full, 1 - self.solid_frac)
        initial_comp_with_solid[0] = self.solid_frac
        self.initial_comp[:] = correct_composition(initial_comp_with_solid, self.comp_min)

        proxy_init_data = self.physics.comp_itor.point_data
        with open('proxy_init_data.pkl', 'wb') as fp:
            pickle.dump(proxy_init_data, fp, protocol=2)

        # Define initial composition for wells
        for i in range(self.nx * self.ny * self.nz, self.nx * self.ny * self.nz + 2):
            self.initial_comp[i, :] = np.array(self.inj_stream)

        for i in range(self.nx * self.ny * self.nz + 2, self.nx * self.ny * self.nz + 4):
            self.initial_comp[i, :] = np.array(self.initial_comp[0, :])

        print('\tNegative composition occurrence while initializing:', self.physics.comp_etor.counter, '\n')
        # ====================================== Initialize reservoir composition ======================================

        # Some newton parameters for non-linear solution:
        self.params.first_ts = 1e-8
        self.params.max_ts = 10

        self.params.tolerance_newton = 1e-3
        self.params.tolerance_linear = 1e-4
        self.params.max_i_newton = 20
        self.params.max_i_linear = 50
        self.params.trans_mult_exp = trans_mult_exp

        self.runtime = 1
        self.timer.node["initialization"].stop()

    # Initialize reservoir and set boundary conditions:
    def set_initial_conditions(self):
        uniform_pressure = self.pressure_init
        uniform_composition = np.array(self.initial_comp)
        self.physics.set_uniform_initial_conditions(self.reservoir.mesh, uniform_pressure, uniform_composition)

    def set_boundary_conditions(self):
        # New boundary condition by adding wells:
        for i, w in enumerate(self.reservoir.wells):
            if "I" in w.name:
                w.control = self.physics.new_rate_oil_inj(self.inj_rate, self.inj_stream)
            else:
                w.control = self.physics.new_bhp_prod(130)


    def evaluate_porosity(self):
        # Initial porosity
        poro_init = 0.25 - self.solid_sat

        print('\tEvaluating porosity...')
        if os.path.exists('proxy_results_data.pkl'):
            with open('proxy_results_data.pkl', 'rb') as fp:
                print('\t\tReading previous point data...')
                self.physics.results_itor.point_data = pickle.load(fp)

        poro = np.zeros(self.reservoir.nb)
        values = value_vector([0] * 2)
        Xm = np.copy(self.physics.engine.X[:self.reservoir.nb*self.physics.n_components])

        for i in range(self.reservoir.nb):
            state = value_vector(
                Xm[i * self.physics.n_components:i * self.physics.n_components + self.physics.n_components])
            self.physics.results_itor.evaluate(state, values)
            poro[i] = 0.25 - values[0]

        poro_diff = poro - poro_init

        proxy_results_data = self.physics.results_itor.point_data
        with open('proxy_results_data.pkl', 'wb') as fp:
            pickle.dump(proxy_results_data, fp, protocol=2)

        print('\tNegative composition while evaluating results:', self.physics.results_etor.counter, '\n')
        return poro_init, poro, poro_diff

    def plot_layer_map(self, map_data, k, name):
        """
        Function to plot parameter profile of certain layer.
        :param map_data: data array
        :param k: layer index
        :param name: parameter name
        """
        import plotly
        import plotly.graph_objs as go
        import plotly.graph_objs.layout._colorscale
        import plotly.express as px

        nxny = self.reservoir.nx * self.reservoir.ny
        map_data_2d = map_data[nxny * (k - 1): nxny * k].reshape(self.reservoir.ny, self.reservoir.nx).transpose()
        data = [go.Heatmap(
            z=np.rot90(map_data_2d, 2),colorscale='Turbo')]
        layout = go.Layout(title='%s, layer %d' % (name, k),
                           yaxis=dict(scaleratio=1, scaleanchor='x', title='X, block'),
                           xaxis=dict(title='Y, block'))
        fig = go.Figure(data=data, layout=layout)
        fig.update_layout(paper_bgcolor='rgba(0,0,0,0)',plot_bgcolor='rgba(0,0,0,0)')
        plotly.offline.plot(fig, filename='%s_%d_map.html' % (name, k))

