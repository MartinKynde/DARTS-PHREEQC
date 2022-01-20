from conversions import bar2atm
from darts.physics import *
import numpy as np
import warnings


class property_container(property_evaluator_iface):
    def __init__(self, phase_name, component_name):
        super().__init__()
        # This class contains all the property evaluators required for simulation
        self.n_phases = len(phase_name)
        self.nc = len(component_name)
        self.component_name = component_name
        self.phase_name = phase_name
        # Allocate (empty) evaluators
        self.density_ev = []
        self.viscosity_ev = []
        self.rel_perm_ev = []
        self.enthalpy_ev = []
        self.kin_rate_ev = []
        self.flash_ev = 0
        self.state_ev = []
        self.init_flash_ev = []


        self.enthalpy_ev = []
        self.rock_energy_ev = []
        self.enthalpy = np.zeros(self.n_phases)

def comp_out_of_bounds(vec_composition, min_z):
    # Check if composition sum is above 1 or element comp below 0, i.e. if point is unphysical:
    temp_sum = 0
    count_corr = 0
    check_vec = np.zeros((len(vec_composition),))

    for ith_comp in range(len(vec_composition)):
        if vec_composition[ith_comp] < min_z:
            vec_composition[ith_comp] = min_z
            count_corr += 1
            check_vec[ith_comp] = 1
        elif vec_composition[ith_comp] > 1 - min_z:
            vec_composition[ith_comp] = 1 - min_z
            temp_sum += vec_composition[ith_comp]
        else:
            temp_sum += vec_composition[ith_comp]

    for ith_comp in range(len(vec_composition)):
        if check_vec[ith_comp] != 1:
            vec_composition[ith_comp] = vec_composition[ith_comp] / temp_sum * (1 - count_corr * min_z)
    return vec_composition

#  Relative permeability based on Corey model

class custom_rel_perm(property_evaluator_iface):
    def __init__(self, exp, sr=0):
        super().__init__()
        self.exp = exp
        self.sr = sr

    def evaluate(self, sat):
        return (sat - self.sr) ** self.exp


class custom_flash:
    def __init__(self, temperature, phreeqc, comp_min):
        self.temperature = temperature - 273.15
        self.phreeqc = phreeqc
        self.comp_min = comp_min

    @staticmethod
    def get_moles(composition):
        # Assume 1000 mol of solution
        total_mole = 1000
        hydrogen_mole = total_mole * composition[-1]
        oxygen_mole = total_mole * composition[-2]
        carbonate_mole = total_mole * composition[-3]
        sulphate_mole = total_mole * composition[-4]
        barium_mole = total_mole * composition[-5]
        return hydrogen_mole, oxygen_mole, carbonate_mole, sulphate_mole, barium_mole

    def get_composition(self, state):
        comp = state[2:]
        # comp[comp <= self.comp_min / 10] = 0
        comp = np.divide(comp, np.sum(comp))
        return comp

    @staticmethod
    def generate_input(*args):
        hydrogen, oxygen, carbonate, sulphate, barium = args[0], args[1], args[2], args[3], args[4]
        if hydrogen / 2 <= oxygen:
            water_mass = hydrogen / 2 * 0.018016
            hydrogen_mole = 0
            oxygen_mole = oxygen - hydrogen / 2
        else:
            water_mass = oxygen * 0.018016
            hydrogen_mole = hydrogen - 2 * oxygen
            oxygen_mole = 0
        carbonate_mole = carbonate
        sulphate_mole = sulphate
        barium_mole = barium

        return water_mass, hydrogen_mole, oxygen_mole, carbonate_mole, sulphate_mole, barium_mole

    def interpret_results(self):
        header = np.array(self.phreeqc.GetSelectedOutputArray())
        results_array = np.array(self.phreeqc.GetSelectedOutputArray()[2])

        # Interpret aqueous phase
        hydrogen_mole_aq = results_array[5]
        oxygen_mole_aq = results_array[6]
        sulphate_mole_aq = results_array[7]
        barium_mole_aq = results_array[8]
        carbonate_mole_aq = results_array[14]

        
        volume_aq = results_array[9] / 1000  # m3
        total_mole_aq = (hydrogen_mole_aq + oxygen_mole_aq + sulphate_mole_aq + barium_mole_aq + carbonate_mole_aq)  # mol
        rho_aq = total_mole_aq / volume_aq / 1000  # kmol/m3

        x = np.array([0,
                      barium_mole_aq / total_mole_aq,
                      sulphate_mole_aq / total_mole_aq,
                      carbonate_mole_aq / total_mole_aq,
                      oxygen_mole_aq / total_mole_aq,
                      hydrogen_mole_aq / total_mole_aq])

        # Interpret gaseous phase
        y = np.zeros(len(x))


        rho_phases = {'aq':rho_aq}
        # Interpret kinetic parameters
        kin_state = {'SI':results_array[10],
                     'SR':results_array[11],
                     'Act(H+)':results_array[12],
                     'mu':results_array[13]}
        return x, y, rho_phases, kin_state

    def evaluate(self, state):
        pressure_atm = bar2atm(state[0])
        comp = self.get_composition(state)
        hydrogen, oxygen,carbonate, sulphate, barium = self.get_moles(comp)
        water_mass, hydrogen_mole, oxygen_mole, carbonate_mole, sulphate_mole, barium_mole = self.generate_input(hydrogen, oxygen, carbonate,
                                                                                                 sulphate, barium)
        string = ("""
        USER_PUNCH
        -headings	H(mol)		 O(mol)		  S(mol)	   Ba(mol)       Vol_aq   SI            SR            ACT("H+") mu  Ca(mol)
        10 PUNCH	TOTMOLE("H") TOTMOLE("O") TOTMOLE("S(6)") TOTMOLE("Ba") SOLN_VOL SI("Barite") SR("Barite") ACT("H+") mu     TOTMOLE("Ca")

        SELECTED_OUTPUT
            -selected_out    true
            -user_punch      true
            -reset           false
            -high_precision  true
            -gases		     CO2(g) H2O(g)
        SOLUTION 1
            temp      %.2f
            pressure  %.4f
            pH        6.2
            units       mg/L
            Alkalinity 40 as HCO3
            Na      51914
            Cl      98582
            -water    %.10f # kg
        REACTION 1
            H         %.10f
            O		  %.10f
            SO4		  %.10f
            Ba        %.10f
            Ca        %.10f      
            1
        KNOBS
            -convergence_tolerance  1e-12
        END""" % (self.temperature, pressure_atm, water_mass, hydrogen_mole, oxygen_mole, sulphate_mole, barium_mole, carbonate_mole))

        self.phreeqc.LoadDatabase('pitzer.dat')
        self.phreeqc.RunString(string)

        x, y, rho_phases, kin_state = self.interpret_results()
        return x, y, rho_phases, kin_state


#units
#mg / L
#Cl
#98642
#Na
#51882
#Ca
#6909
#Mg
#1564
#Ba
#18
#S(6)
#69
class init_flash(custom_flash):
    def __init__(self, temperature, phreeqc, comp_min):
        super().__init__(temperature, phreeqc, comp_min)

    @staticmethod
    def get_moles(comp_data):
        # Assume 1000 mol of solution
        total_mole = comp_data[-1]
        hydrogen_mole = total_mole * comp_data[-2]
        oxygen_mole = total_mole * comp_data[-3]
        carbonate_mole = total_mole * comp_data[-4]
        sulphate_mole = total_mole * comp_data[-5]
        barium_mole = total_mole * comp_data[-6]
        return hydrogen_mole, oxygen_mole, carbonate_mole, sulphate_mole, barium_mole

    def evaluate(self, state):
        pressure_atm = bar2atm(state[0])
        non_solid_composition = np.divide(state[3:], np.sum(state[3:]))
        comp_data = np.hstack((non_solid_composition, state[1]))
        hydrogen, oxygen, carbonate, sulphate, barium = self.get_moles(comp_data)
        water_mass, hydrogen_mole, oxygen_mole, carbonate_mole, sulphate_mole, barium_mole = self.generate_input(hydrogen, oxygen, carbonate,
                                                                                                 sulphate, barium)
        string = ("""
        USER_PUNCH
        -headings	Vol_aq(L)
        10 PUNCH	SOLN_VOL

        SELECTED_OUTPUT
            -selected_out    true
            -user_punch      true
            -reset           false
            -high_precision  true
        SOLUTION 1
            temp      %.2f
            pressure  %.4f
            pH        6.2
            units       mg/L
            Alkalinity 40 as HCO3
            Na      51914
            Cl      98582
            -water    %.3f # kg
        REACTION 1
            H         %.8f
            O		  %.8f
            SO4		  %.8f
            Ba        %.8f
            Ca        %.10f
            1
        END""" % (self.temperature, pressure_atm, water_mass, hydrogen_mole, oxygen_mole, sulphate_mole, barium_mole, carbonate_mole))

        self.phreeqc.LoadDatabase('pitzer.dat')
        self.phreeqc.RunString(string)


        non_solid_volume = np.array(self.phreeqc.GetSelectedOutputArray()[2]) / 1000  # m3
        return non_solid_volume


import math


class custom_kinetic_rate:
    def __init__(self, temperature, comp_min):
        self.temperature = temperature
        self.comp_min = comp_min

    def evaluate(self, kin_state, solid_saturation, rho_s, min_z, kin_fact):
        # Define constants
        specific_sa = 5e-04  # [m2/mol], default = 0.925

        A = 1.43e-04  # [mol * m-2 * s-1]
        Ea = 25000  # [J * mol-1]
        nH = 0.03
        nI = 0.6
        R = 8.314472  # gas constant [J/mol/Kelvin]
        n = 1

        # Define rate parameters
        sat_index = kin_state['SI']
        sat_ratio = kin_state['SR']
        hydrogen_act = kin_state['Act(H+)']
        ionic_strength = kin_state['mu']


        k = A * math.exp(-Ea / (R * (self.temperature + 273.15)))
        kinetic_rate = k * specific_sa * (solid_saturation * rho_s * 1000) * hydrogen_act ** nH * (1 - sat_ratio) * 10 ** (nI * math.sqrt(ionic_strength))
        # [mol/s]
        # Convert to [kmol/d]
        kinetic_rate = kinetic_rate * 60 * 60 * 24 / 1000
        
        return kinetic_rate


class custom_state(property_evaluator_iface):
    def __init__(self):
        super().__init__()

    def evaluate(self, zc, V):
        solid_present = zc[-1] > 1e-8
        vapor_present = V > 1e-8
        liquid_present = V < 1 - 1e-8 - zc[-1]
        if liquid_present and vapor_present and solid_present:
            # Three phase system:
            state = '111'
        elif liquid_present and vapor_present:
            # Liquid and Vapor:
            state = '110'
        elif liquid_present and solid_present:
            # Liquid and Solid:
            state = '101'
        elif vapor_present and solid_present:
            # Vapor and Solid:
            state = '011'
        elif liquid_present:
            # Single phase liquid:
            state = '100'
        elif vapor_present:
            # Single phase vapor:
            state = '010'
        elif solid_present:
            # Single phase solid:
            state = '001'
        else:
            # Something went wrong:
            state = '000'
        return state
