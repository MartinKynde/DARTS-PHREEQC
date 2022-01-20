from martin_operator_evaluator import my_own_acc_flux_etor, my_own_rate_evaluator, my_own_comp_etor, my_own_results_etor
from martin_properties import *

from darts.engines import *
import numpy as np
from darts.models.physics.physics_base import PhysicsBase

# Define our own operator evaluator class
class OwnPhysicsClass(PhysicsBase):
    def __init__(self, timer, components, n_points, min_p, max_p, min_z, input_data_struct, properties,
                 platform='cpu', itor_type='multilinear', itor_mode='adaptive', itor_precision='d', cache=True):
        # Obtain properties from user input during initialization:
        super().__init__(cache)
        self.timer = timer.node["simulation"]
        self.n_points = n_points
        self.min_p = min_p
        self.max_p = max_p
        self.min_z = min_z
        self.components = components
        self.n_components = len(components)
        self.phases = ['vapor', 'liquid']
        self.n_phases = len(self.phases)
        self.n_vars = self.n_components
        NE = self.n_components
        self.vars = ["p"] + components[:-1]
        self.n_axes_points = index_vector([n_points] * self.n_vars)
        self.n_axes_min = value_vector([min_p] + [min_z] * (self.n_components - 1))
        self.n_axes_max = value_vector([max_p] + [1 - min_z] * (self.n_components - 1))

        # Engine initialization
        # engine_name = eval("engine_nc_kin_dif_cpu%d" % self.nr_components)
        self.engine = eval("engine_super_%s%d_%d" % (platform, self.n_components, self.n_phases))()
        self.n_ops = NE + self.n_phases * NE + self.n_phases + self.n_phases * NE + NE + 3 + 2 * self.n_phases + 1

        # Initialize main evaluator
        self.acc_flux_etor = my_own_acc_flux_etor(input_data_struct, properties)

        # Create actual accumulation and flux interpolator:
        self.acc_flux_itor = self.create_interpolator(self.acc_flux_etor, self.n_vars, self.n_ops, self.n_axes_points,
                                                      self.n_axes_min, self.n_axes_max, platform=platform,
                                                      algorithm=itor_type, mode=itor_mode,
                                                      precision=itor_precision)

        # ==============================================================================================================

        # Create initialization evaluator
        self.comp_etor = my_own_comp_etor(input_data_struct.pressure_init, properties.init_flash_ev)

        self.n_axes_points2 = index_vector([2 * n_points] * self.n_vars)
        self.n_axes_min2 = value_vector([0] + [min_z] * (self.n_components - 1))
        self.n_axes_max2 = value_vector([1] + [1 - min_z] * (self.n_components - 1))

        self.comp_itor = self.create_interpolator(self.comp_etor, self.n_vars, 2, self.n_axes_points2,
                                                  self.n_axes_min2, self.n_axes_max2, platform=platform,
                                                  algorithm=itor_type, mode=itor_mode,
                                                  precision=itor_precision)

        # ==============================================================================================================

        # Create results evaluator
        self.results_etor = my_own_results_etor(input_data_struct, properties)

        # Initialize results interpolator
        self.results_itor = self.create_interpolator(self.results_etor, self.n_vars, self.n_phases, self.n_axes_points2,
                                                     self.n_axes_min, self.n_axes_max, platform=platform,
                                                     algorithm=itor_type, mode=itor_mode,
                                                     precision=itor_precision)

        # ==============================================================================================================

        # Create rate evaluator and interpolator:
        self.rate_etor = my_own_rate_evaluator(properties, input_data_struct.temperature, input_data_struct.c_r)

        self.rate_itor = self.create_interpolator(self.rate_etor, self.n_vars, self.n_phases, self.n_axes_points,
                                                  self.n_axes_min, self.n_axes_max, platform=platform,
                                                  algorithm=itor_type, mode=itor_mode,
                                                  precision=itor_precision)

        self.create_itor_timers(self.acc_flux_itor, 'reservoir interpolation')
        self.create_itor_timers(self.comp_itor, 'comp interpolation')
        self.create_itor_timers(self.results_itor, 'results interpolation')
        self.create_itor_timers(self.rate_itor, 'rate interpolation')

        # define well control factories
        # Injection wells (upwind method requires both bhp and inj_stream for bhp controlled injection wells):
        self.new_bhp_inj = lambda bhp, inj_stream: bhp_inj_well_control(bhp, value_vector(inj_stream))
        self.new_rate_gas_inj = lambda rate, inj_stream: rate_inj_well_control(self.phases, 0, self.n_components,
                                                                               self.n_components, rate,
                                                                               value_vector(inj_stream), self.rate_itor)
        self.new_rate_oil_inj = lambda rate, inj_stream: rate_inj_well_control(self.phases, 1, self.n_components,
                                                                               self.n_components, rate,
                                                                               value_vector(inj_stream), self.rate_itor)
        # Production wells:
        self.new_bhp_prod = lambda bhp: bhp_prod_well_control(bhp)
        self.new_rate_gas_prod = lambda rate: rate_prod_well_control(self.phases, 0, self.n_components,
                                                                     self.n_components,
                                                                     rate, self.rate_itor)
        self.new_rate_oil_prod = lambda rate: rate_prod_well_control(self.phases, 1, self.n_components,
                                                                     self.n_components,
                                                                     rate, self.rate_itor)

    # Define some class methods:
    def init_wells(self, wells):
        for w in wells:
            assert isinstance(w, ms_well)
            w.init_rate_parameters(self.n_components, self.phases, self.rate_itor)

    def set_uniform_initial_conditions(self, mesh, uniform_pressure, uniform_composition):
        assert isinstance(mesh, conn_mesh)

        nb = mesh.n_blocks

        # set inital pressure
        pressure = np.array(mesh.pressure, copy=False)
        pressure.fill(uniform_pressure)

        # set initial composition
        mesh.composition.resize(nb * (self.n_components - 1))
        composition = np.array(mesh.composition, copy=False)
        for c in range(self.n_components - 1):
            composition[c::(self.n_components - 1)] = uniform_composition[:, c]

    def set_boundary_conditions(self, mesh, uniform_pressure, uniform_composition):
        assert isinstance(mesh, conn_mesh)

        # Class methods which can create constant pressure and composition boundary condition:
        pressure = np.array(mesh.pressure, copy=False)
        pressure.fill(uniform_pressure)

        mesh.composition.resize(mesh.n_blocks * (self.n_components - 1))
        composition = np.array(mesh.composition, copy=False)
        for c in range(self.n_components - 1):
            composition[c::(self.n_components - 1)] = uniform_composition[c]
