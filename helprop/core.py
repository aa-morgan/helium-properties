# -*- coding: utf-8 -*-
"""
@author: Alex Morgan
"""

import hsML
import numpy as np

from scipy.constants import h, c, m_e, e, pi, epsilon_0, physical_constants
a_0 = physical_constants['Bohr radius'][0]
E_h = physical_constants['Hartree energy'][0]
au_to_debye = 1/0.393430307

def State(n, L, S, ML):
    J = max(L, S)
    return hsML.State(n, L, S, J, ML)

def energy(state, **kwargs):
    qd = hsML.get_qd(state.n, state.L, state.S, state.J)
    n_eff = state.n - qd
    en = hsML.energy(state.n, n_eff)
    units = select_units(unit_type='energy', units=kwargs.get('units', 'ghz'))
    return en * units

def transition_energy(state_1, state_2, **kwargs):
    en_1 = energy(state_1, units='au')
    en_2 = energy(state_2, units='au')
    transition_en = np.abs(en_1 - en_2)
    units = select_units(unit_type='energy', units=kwargs.get('units', 'ghz'))
    return transition_en * units

def transition_dipole_moment(state_1, state_2, **kwargs):
    """ Calculates the electric transition dipole moment between two states
        Returns the value in Debye
    """
    dM_allow = kwargs.get('dM_allow', [0,+1,-1])
    tdm = np.abs(hsML.stark_int(state_1, state_2, dM_allow=[0,+1,-1]))
    units = select_units(unit_type='electric dipole moment', units=kwargs.get('units', 'debye'))
    return tdm * units

def spontaneous_emission_rate(state_1, state_2, **kwargs):
    if not(allowed_transition(state_1, state_2, **kwargs)):
        return 'Not allowed transition'
    en_1 = energy(state_1, units='j')
    en_2 = energy(state_2, units='j')
    omega_transition = 2*pi*np.abs(en_1 - en_2)/h # E=hf=h(w/2pi), w=2pi*E/h
    
    tdm = transition_dipole_moment(state_1, state_2, units='cm')
    einstein_A = ((2*omega_transition**3)/(3*epsilon_0*h*c**3)) * np.abs(tdm)**2
    units = select_units(unit_type='time', units=kwargs.get('units', 's'))
    return einstein_A / units

def spontaneous_emission_rate_all(state_1, **kwargs):
    basis = hsML.Hamiltonian(n_min=1, n_max=state_1.n, S=state_1.S).basis
    emission_rate_sum = 0.0
    allowed_states = []
    for state_2 in basis:
        if allowed_transition(state_1, state_2, **kwargs):
            allowed_states.append(state_2)
            emission_rate_sum += spontaneous_emission_rate(state_1, state_2, **kwargs)
    return emission_rate_sum, allowed_states

def radiative_lifetime(state_1, state_2, **kwargs):
    if not(allowed_transition(state_1, state_2, **kwargs)):
        return 'Not allowed transition'
    return 1/spontaneous_emission_rate(state_1, state_2, **kwargs)

def radiative_lifetime_all(state, **kwargs):
    emission_rate_sum, allowed_states = spontaneous_emission_rate_all(state, **kwargs)
    return 1/emission_rate_sum, allowed_states

def allowed_transition(state_1, state_2, **kwargs):
    dL = state_2.L - state_1.L
    dS = state_2.S - state_1.S
    dML = state_2.ML - state_1.ML
    dEn = energy(state_2) - energy(state_1)
    if  abs(dL) == 1 and \
        dS == 0 and \
        dML in kwargs.get('dM_allow', [0,+1,-1]):
        return True
    else:
        return False

def select_units(unit_type, units):
    if unit_type == 'energy':
        if str(units).lower() in ['au', 'atomic']:
            return 1.0
        elif str(units).lower() in ['j', 'joules', 'E']:
            return E_h
        elif str(units).lower() in ['cm', 'wavenumber', 'wavenumbers']:
            return (E_h/(h*c))*100
        elif str(units).lower() in ['hz', 'frequency', 'freq', 'f']:
            return (E_h/h)
        elif str(units).lower() in ['ghz']:
            return (E_h/h) * 10**-9
    elif unit_type == 'electric dipole moment':
        if str(units).lower() in ['au', 'atomic']:
            return 1.0
        elif str(units).lower() in ['debye', 'd']:
            return au_to_debye
        elif str(units).lower() in ['coulomb meter', 'cm']:
            return e*a_0
    elif unit_type == 'time':
        if str(units).lower() in ['s']:
            return 1.0
        elif str(units).lower() in ['ms']:
            return 10**3
        elif str(units).lower() in ['us']:
            return 10**6
        elif str(units).lower() in ['ns']:
            return 10**9
    else:
        raise('Units '+str(units)+' not known.')