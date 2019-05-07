# -*- coding: utf-8 -*-
"""
@author: Alex Morgan
"""

import hsz
import numpy as np
from .constants import *

def State(n, L, S, M, basis_type='ml', J=None):
    if basis_type.lower() == 'ml':
        _J = max(L, S)
    elif basis_type.lower() == 'mj':
        _J = J
    else:
        raise Exception("basis_type '{}' not recognised!".format(basis_type))
    return hsz.State(n, L, S, _J, M)

def energy(state, **kwargs):
    qd = hsz.get_qd(state.n, state.L, state.S, state.J)
    n_eff = state.n - qd
    en = hsz.energy(state.n, n_eff)
    return convert_units(value=en, unit_type='energy', units=kwargs.get('units', 'ghz'))

def transition_energy(state_1, state_2, **kwargs):
    en_1 = energy(state_1, units='au')
    en_2 = energy(state_2, units='au')
    transition_en = np.abs(en_1 - en_2)
    return convert_units(value=transition_en, unit_type='energy', units=kwargs.get('units', 'ghz'))

def transition_dipole_moment(state_1, state_2, basis_type='ml', **kwargs):
    """ Calculates the electric transition dipole moment between two states
        Returns the value in Debye
    """
    dM_allow = kwargs.get('dM_allow', [0,+1,-1])
    tdm = np.abs(hsz.stark_interaction(state_1, state_2, basis_type, dM_allow=[0,+1,-1]))
    return convert_units(value=tdm, unit_type='electric dipole moment', units=kwargs.get('units', 'debye'))

def spontaneous_emission_rate(state_1, state_2, **kwargs):
    if not(allowed_transition(state_1, state_2, **kwargs)):
        return 'Not allowed transition'
    en_1 = energy(state_1, units='j')
    en_2 = energy(state_2, units='j')
    omega_transition = 2*pi*np.abs(en_1 - en_2)/h # E=hf=h(w/2pi), w=2pi*E/h
    
    tdm = transition_dipole_moment(state_1, state_2, units='cm')
    einstein_A = ((2*omega_transition**3)/(3*epsilon_0*h*c**3)) * np.abs(tdm)**2
    return convert_units(value=einstein_A, unit_type='time', units=kwargs.get('units', 's'))

def spontaneous_emission_rate_all(state_1, **kwargs):
    basis = hsz.HamiltonianMatrix(n_min=1, n_max=state_1.n, S=state_1.S).basis
    emission_rate_sum = 0.0
    allowed_states = []
    for state_2 in basis.states:
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
    dM = state_2.M - state_1.M
    dEn = energy(state_2) - energy(state_1)
    if  abs(dL) == 1 and \
        dS == 0 and \
        dM in kwargs.get('dM_allow', [0,+1,-1]):
        return True
    else:
        return False
    
def hydrogenic_polarizability(state, **kwargs):
    n = state.n
    k = kwargs.get('k', 0)
    m = state.M
    F_0 = En_h/(e*a_0)
    pol = (1/8)*(1/F_0)**2*n**4*(17*n**2 - 3*k**2 - 9*m**2 + 19)
    return convert_units(value=pol, unit_type='polarizability', units=kwargs.get('units', 'hz/(v/cm)2'))

def relative_hydrogenic_polarizability(state_1, state_2, **kwargs):
    pol_1 = hydrogenic_polarizability(state_1, **kwargs)
    pol_2 = hydrogenic_polarizability(state_2, **kwargs)
    return pol_2 - pol_1

def convert_units(value, unit_type, units):
    if unit_type == 'energy':
        if str(units).lower() in ['au', 'atomic', 'ea']:
            return value * 1.0
        elif str(units).lower() in ['j', 'joules', 'E']:
            return value * En_h
        elif str(units).lower() in ['cm', 'wavenumber', 'wavenumbers']:
            return value * (En_h/(h*c*100))
        elif str(units).lower() in ['hz', 'frequency', 'freq', 'f']:
            return (value * En_h/h)
        elif str(units).lower() in ['khz']:
            return value * (En_h/h) * 10**-3
        elif str(units).lower() in ['mhz']:
            return value * (En_h/h) * 10**-6
        elif str(units).lower() in ['ghz']:
            return value * (En_h/h) * 10**-9
        elif str(units).lower() in ['m']:
            return abs((c/(En_h/h)) / value)
        elif str(units).lower() in ['cm']:
            return abs((c/(En_h/h)) * 10**2 / value)
        elif str(units).lower() in ['mm']:
            return abs((c/(En_h/h)) * 10**3 / value)
        elif str(units).lower() in ['um']:
            return abs((c/(En_h/h)) * 10**6 / value)
        elif str(units).lower() in ['nm']:
            return abs((c/(En_h/h)) * 10**9 / value)
    elif unit_type == 'electric dipole moment':
        if str(units).lower() in ['au', 'atomic']:
            return value * 1.0
        elif str(units).lower() in ['debye', 'd']:
            return value * au_to_debye
        elif str(units).lower() in ['coulomb meter', 'cm']:
            return value * e*a_0
    elif unit_type == 'time':
        if str(units).lower() in ['s', 'sec', 'seconds']:
            return value * 1.0
        elif str(units).lower() in ['ms', 'millisec', 'milliseconds']:
            return value * 10**-3
        elif str(units).lower() in ['us', 'microsec', 'microseconds']:
            return value * 10**-6
        elif str(units).lower() in ['ns', 'nanosec', 'nanoseconds']:
            return value * 10**-9
    elif unit_type == 'polarizability':
        if str(units).lower() in ['hz/(v/m)2']:
            return value * En_h/(2*h)
        elif str(units).lower() in ['khz/(v/m)2']:
            return value * En_h/(2*h) *10**-3
        elif str(units).lower() in ['mhz/(v/m)2']:
            return value * En_h/(2*h) *10**-6
        if str(units).lower() in ['hz/(v/cm)2']:
            return value * En_h/(2*h) * 10**4
        elif str(units).lower() in ['khz/(v/cm)2']:
            return value * En_h/(2*h) *10**-3 * 10**4
        elif str(units).lower() in ['mhz/(v/cm)2']:
            return value * En_h/(2*h) *10**-6 * 10**4
        elif str(units).lower() in ['cm2/v']:
            return value * En_h
    else:
        raise('Units '+str(units)+' not known.')