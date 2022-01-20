import CoolProp.CoolProp as CP


from conversions import bar2pa
from darts.physics import *
from darts.engines import *
import os.path as osp
import numpy as np

physics_name = osp.splitext(osp.basename(__file__))[0]

# Define our own operator evaluator class
class my_own_acc_flux_etor(operator_set_evaluator_iface):
    def __init__(self, input_data, properties, thermal=1):
        super().__init__()  # Initialize base-class
        # Store your input parameters in self here, and initialize other parameters here in self
        self.input_data = input_data
        self.num_comp = input_data.num_comp
        self.min_z = input_data.min_z
        self.temperature = input_data.temperature
        self.exp_w = input_data.exp_w
        self.exp_g = input_data.exp_g
        self.c_r = input_data.c_r
        self.kin_fact = input_data.kin_fact
        self.property = properties
        self.thermal = thermal
        self.counter = 0

    def comp_out_of_bounds(self, vec_composition):
        # Check if composition sum is above 1 or element comp below 0, i.e. if point is unphysical:
        temp_sum = 0
        count_corr = 0
        check_vec = np.zeros((len(vec_composition),))

        for ith_comp in range(len(vec_composition)):
            if vec_composition[ith_comp] < self.min_z:
                vec_composition[ith_comp] = self.min_z
                count_corr += 1
                check_vec[ith_comp] = 1
            elif vec_composition[ith_comp] > 1 - self.min_z:
                vec_composition[ith_comp] = 1 - self.min_z
                temp_sum += vec_composition[ith_comp]
            else:
                temp_sum += vec_composition[ith_comp]

        for ith_comp in range(len(vec_composition)):
            if check_vec[ith_comp] != 1:
                vec_composition[ith_comp] = vec_composition[ith_comp] / temp_sum * (1 - count_corr * self.min_z)
        return vec_composition

    def evaluate(self, state, values):
        """
        Class methods which evaluates the state operators for the element based physics
        :param state: state variables [pres, comp_0, ..., comp_N-1]
        :param values: values of the operators (used for storing the operator values)
        :return: updated value for operators, stored in values
        """

        # Composition vector and pressure from state:
        vec_state_as_np = np.asarray(state)
        pressure = vec_state_as_np[0]

        zc = np.append(vec_state_as_np[1:], 1 - np.sum(vec_state_as_np[1:]))

        # # # ================================ Flash ================================ # # #
        # Check for negative composition occurrence
        if (zc < 0).any():
            self.counter += 1
            zc_copy = np.array(zc)
            zc_copy = self.comp_out_of_bounds(zc_copy)
            flash_state = np.append(pressure, zc_copy)
        else:
            flash_state = np.append(pressure, zc)

        # Perform Flash procedure here:
        x, y, rho_phases, kin_state = self.property.flash_ev.evaluate(flash_state)

        # NOTE: z_CaCO3 = zc[-1] = 1 - V - L (since only CaCO3 appears in solid phase)
        sol = zc[0]

        # Calculate liquid phase fraction:
        liq = 1 - sol

        # Note: officially three phases are present now
        rho_w = rho_phases['aq']
        mu_w = CP.PropsSI('V', 'T', self.temperature, 'P|liquid', bar2pa(pressure), 'Water') * 1000


        # in kmol/m3
        rho_s = (1 + self.c_r * (pressure - 1)) * 4500 / 233.38
        # # # ================================ Flash ================================ # # #

        # Get saturations
        sg = 0
        sw = liq / rho_w / (liq / rho_w + sol / rho_s)
        ss = sol / rho_s / (liq / rho_w + sol / rho_s)

        # Need to normalize to get correct Brook-Corey relative permeability
        sg_norm = sg / (sg + sw)
        sw_norm = sw / (sg + sw)

        kr_w = self.property.rel_perm_ev['liq'].evaluate(sw_norm)

        # all properties are in array, and can be separate
        self.x = np.array([y, x])
        self.rho_m = np.array([rho_w])
        self.kr = np.array([kr_w])
        self.mu = np.array([mu_w])
        self.compr = 1 + self.c_r * (pressure - 1)
        self.sat = np.array([sw_norm])

        # Total density
        rho_t = rho_w * sw + rho_s * ss

        # Kinetic reaction rate
        kin_rate = self.property.kin_rate_ev.evaluate(kin_state, ss, rho_s, self.min_z, self.kin_fact)

        nc = self.input_data.num_comp
        nph = 2
        ne = nc

        #       al + bt        + gm + dlt + chi     + rock_temp por    + gr/cap  + por + temp
        total = ne + ne * nph + nph + ne + ne * nph + 3 + 2 * nph + 1

        for i in range(total):
            values[i] = 0

        vec_state_as_np = np.asarray(state)
        pressure = vec_state_as_np[0]

        nc = self.input_data.num_comp
        nph = 2
        ne = nc

        #       al + bt        + gm + dlt + chi     + rock_temp por    + gr/cap  + por
        total = ne + ne * nph + nph + ne + ne * nph + 3 + 2 * nph + 1

        for i in range(total):
            values[i] = 0


        """ CONSTRUCT OPERATORS HERE """

        """ Alpha operator represents accumulation term: """
        for i in range(self.num_comp):
            values[i] = rho_t * zc[i]

        """ Beta operator represents flux term: """
        for j in range(nph):
            shift = ne + ne * j
            for i in range(nc):
                values[shift + i] = self.x[j][i] * self.rho_m * self.kr / self.mu

        """ Gamma operator for diffusion (same for thermal and isothermal) """
        dif_coef = np.array([0, 1, 1, 1, 1, 1]) * 5.2e-10 * 86400
        
        shift = ne + ne * nph
        for j in range(nph):
            values[shift + j] = self.compr * self.sat 

        """ Chi operator for diffusion """
        shift += nph
        for i in range(nc):
            for j in range(nph):
                values[shift + i * nph + j] = dif_coef[i] * self.rho_m * self.x[j][i]

        """ Delta operator for reaction """
        shift += nph * ne
        for i in range(ne):
            values[shift + i] = self.input_data.stoich_matrix[i] * kin_rate

        """ Gravity and Capillarity operators """
        shift += ne
        # E3-> gravity
        for i in range(nph):
            values[shift + 3 + i] = 0

        # E4-> capillarity
        for i in range(nph):
            values[shift + 3 + nph + i] = 0
        # E5_> porosity
        values[shift + 3 + 2 * nph] = 1-ss



        return 0


class my_own_rate_evaluator(operator_set_evaluator_iface):
    # Simplest class existing to mankind:
    def __init__(self, properties, temperature, c_r):
        # Initialize base-class
        super().__init__()
        self.property = properties
        self.temperature = temperature
        self.c_r = c_r



    def comp_out_of_bounds(self, vec_composition):
        # Check if composition sum is above 1 or element comp below 0, i.e. if point is unphysical:
        temp_sum = 0
        count_corr = 0
        check_vec = np.zeros((len(vec_composition),))

        for ith_comp in range(len(vec_composition)):
            if vec_composition[ith_comp] < self.min_z:
                vec_composition[ith_comp] = self.min_z
                count_corr += 1
                check_vec[ith_comp] = 1
            elif vec_composition[ith_comp] > 1 - self.min_z:
                vec_composition[ith_comp] = 1 - self.min_z
                temp_sum += vec_composition[ith_comp]
            else:
                temp_sum += vec_composition[ith_comp]

        for ith_comp in range(len(vec_composition)):
            if check_vec[ith_comp] != 1:
                vec_composition[ith_comp] = vec_composition[ith_comp] / temp_sum * (1 - count_corr * self.min_z)
        return vec_composition

    def evaluate(self, state, values):
        # Composition vector and pressure from state:


        vec_state_as_np = np.asarray(state)
        pressure = vec_state_as_np[0]

        zc = np.append(vec_state_as_np[1:], 1 - np.sum(vec_state_as_np[1:]))

        # # # ================================ Flash ================================ # # #
        # Check for negative composition occurrence
        if (zc < 0).any():
            self.counter += 1
            zc_copy = np.array(zc)
            zc_copy = self.comp_out_of_bounds(zc_copy)
            flash_state = np.append(pressure, zc_copy)
        else:
            flash_state = np.append(pressure, zc)

        # Perform Flash procedure here:
        x, y, rho_phases, kin_state = self.property.flash_ev.evaluate(flash_state)

        # NOTE: z_CaCO3 = zc[-1] = 1 - V - L (since only CaCO3 appears in solid phase)
        sol = zc[0]

        # Calculate liquid phase fraction:
        liq = 1 - sol

        # Note: officially three phases are present now
        rho_w = rho_phases['aq']
        mu_w = CP.PropsSI('V', 'T', self.temperature, 'P|liquid', bar2pa(pressure), 'Water') * 1000

        # in kmol/m3
        rho_s = (1 + self.c_r * (pressure - 1)) * 4500 / 233.38
        # # # ================================ Flash ================================ # # #

        # Get saturations

        sg = 0
        sw = liq / rho_w / (liq / rho_w + sol / rho_s)
        ss = sol / rho_s / (liq / rho_w + sol / rho_s)

        # Need to normalize to get correct Brook-Corey relative permeability
        sg_norm = sg / (sg + sw)
        sw_norm = sw / (sg + sw)

        kr_w = self.property.rel_perm_ev['liq'].evaluate(sw_norm)

        # all properties are in array, and can be separate
        self.x = np.array([y, x])
        self.rho_m = np.array([rho_w])
        self.kr = np.array([kr_w])
        self.mu = np.array([mu_w])
        self.compr = 1 + self.c_r * (pressure - 1)
        self.sat = np.array([sw_norm])

        # Total density
        rho_t = rho_w * sw + rho_s * ss

        # Easiest example, constant volumetric phase rate:
        values[0] = 0   # vapor phase
        values[1] = (rho_w * kr_w / mu_w) / rho_t   # liquid phase
        # Usually some steps will be executed to estimate unit volumetric flow rate based on current state (when
        # multiplied with a pressure difference one obtains actual volumetric flow rate)
        return 0


# Define our own operator evaluator class for initialization
class my_own_comp_etor(operator_set_evaluator_iface):
    def __init__(self, pressure, flash_ev):
        super().__init__()  # Initialize base-class
        self.pressure = pressure
        self.flash_ev = flash_ev
        self.non_solid_mole = 1000
        self.counter = 0

    def comp_out_of_bounds(self, vec_composition):
        # Check if composition sum is above 1 or element comp below 0, i.e. if point is unphysical:
        temp_sum = 0
        count_corr = 0
        check_vec = np.zeros((len(vec_composition),))

        for ith_comp in range(len(vec_composition)):
            if vec_composition[ith_comp] < self.min_z:
                vec_composition[ith_comp] = self.min_z
                count_corr += 1
                check_vec[ith_comp] = 1
            elif vec_composition[ith_comp] > 1 - self.min_z:
                vec_composition[ith_comp] = 1 - self.min_z
                temp_sum += vec_composition[ith_comp]
            else:
                temp_sum += vec_composition[ith_comp]

        for ith_comp in range(len(vec_composition)):
            if check_vec[ith_comp] != 1:
                vec_composition[ith_comp] = vec_composition[ith_comp] / temp_sum * (1 - count_corr * self.min_z)
        return vec_composition

    def evaluate(self, state, values):
        vec_state_as_np = np.asarray(state)
        zc = np.append(vec_state_as_np[1:], 1 - np.sum(vec_state_as_np[1:]))

        # Check for negative composition occurrence
        if (zc < 0).any():
            self.counter += 1
            zc_copy = np.array(zc)
            zc_copy = self.comp_out_of_bounds(zc_copy)
            flash_state = np.hstack((self.pressure, self.non_solid_mole, zc_copy))
        else:
            flash_state = np.hstack((self.pressure, self.non_solid_mole, zc))

        # Perform Flash procedure here:
        non_sol_volume = self.flash_ev.evaluate(flash_state)    # m3

        # vec_state_as_np[0] is solid saturation (volume fraction)
        total_volume = non_sol_volume / (1 - vec_state_as_np[0])
        sol_volume = total_volume * vec_state_as_np[0]  # m3
        sol_mole = sol_volume * (1 + 4.5e-5 * (self.pressure - 1)) * 4500 / 0.23338
        sol = sol_mole / (sol_mole + self.non_solid_mole)
        values[0] = sol

        return 0


# Define results evaluator
class my_own_results_etor(my_own_acc_flux_etor):
    def __init__(self, input_data, properties):
        super().__init__(input_data, properties)

    def evaluate(self, state, values):
        # Composition vector and pressure from state:
        vec_state_as_np = np.asarray(state)
        pressure = vec_state_as_np[0]
        zc = np.append(vec_state_as_np[1:], 1 - np.sum(vec_state_as_np[1:]))

        # # # ================================ Flash ================================ # # #
        # Check for negative composition occurrence
        if (zc < 0).any():
            self.counter += 1

            zc_copy = np.array(zc)
            for i in range(len(zc_copy)):
                if zc_copy[i] < 0:
                    zc_copy[i] = 0
            zc_copy = np.divide(zc_copy, np.sum(zc_copy))

            flash_state = np.append(pressure, zc_copy)
        else:
            flash_state = np.append(pressure, zc)

        # Perform Flash procedure here:
        _, _, rho_phases, _ = self.property.flash_ev.evaluate(flash_state)

        # NOTE: z_CaCO3 = zc[0] = 1 - V - L (since only CaCO3 appears in solid phase)
        sol = zc[0]

        # Calculate liquid phase fraction:
        liq = 1 - sol

        # Note: officially three phases are present now
        rho_w = rho_phases['aq']
        rho_s = (1 + self.c_r * (pressure - 1)) * 4500 / 233.38
        # # # ================================ Flash ================================ # # #

        # Get saturations

        ss = sol / rho_s / (liq / rho_w + sol / rho_s)

        values[0] = ss
        return 0
