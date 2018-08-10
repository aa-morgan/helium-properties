`Rydberg Helium Properties`
===============

Calculate properties of Rydberg states in helium using the Numerov method.

M. L. Zimmerman et al., Phys. Rev. A, 20, 2251 (1979)
http://dx.doi.org/10.1103/PhysRevA.20.2251

Install
------- 

Install `helium-properties` using `setuptools`,
```bash
git clone https://github.com/axm108/helium-properties
cd helium-properties
python setup.py install
```

Install [`helium-stark-zeeman`](https://github.com/axm108/helium-stark-zeeman) using `setuptools`,
```bash
git clone https://github.com/axm108/helium-stark-zeeman
cd helium-stark-zeeman
python setup.py install
```

Basic usage
-------
Import libraries,
```python
from helprop import *
import pandas as pd
```
Instantiate `State` objects,
```python
state_1 = State(n=68, L=67, S=1, M=67, basis_type='ML')
state_2 = State(n=69, L=68, S=1, M=68, basis_type='ML')
```
Calculate single-state properties,
```python
energy(state_1, units='Ghz')
radiative_lifetime_all(state_1, units='ms')
```
Calculate two-state properties,
```python
transition_energy(state_1, state_2, units='ghz'),
transition_dipole_moment(state_1, state_2, units='debye'),
radiative_lifetime(state_1, state_2, units='ms')
```

Parameters
-------

#### Class: `State`
| Parameter  | Description | Data type | Required | Default |
| ------------ | ------------ | ------------ | ------------ | ------------ |
| `n_min` | Minimum value of the principle quantum number `$n$`, to allow in the basis. | `Int` | Yes | N/A |
| `n_max` | Maximum value of the principle quantum number `$n$`, to allow in the basis. | `Int` | Yes | N/A |
| `L_max` | Maximum value of the orbital angular momentum quantum number to allow in the basis. `None` means no restriction. | `Int` [or `None`] | No | `None` |
| `S` | Value of the total spin orbital angular momentum. [Singlet: `S=0`, Triplet: `S=1`]. `None` means both.  | `Int` [or `None`] | No | `None` |
| `M` | Single value of the azimuthal quantum number to use in the basis. `None` means no restriction. | `Int` [or `None`] | No | `None` |
| `M_max` | Maximum value of the azimuthal quantum number to allow in the basis. `None` means no restriction. | `Int` [or `None`] | No | `None` |
| `basis_type` | Whether to use the `$|n, \ell, S, M_{\ell} \rangle$`, or `$|n, \ell, S, J, M_J \rangle$` basis. Specify using `'ML'`, or `'MJ'`, respectively. | `String` | No | `'ML'` |

#### Method: `energy`
Energy relative to the ionisation energy.

| Parameter  | Description | Data type | Required | Default |
| ------------ | ------------ | ------------ | ------------ | ------------ |
| `state` | State to calculate the energy of. | `State` | Yes | N/A |
| `units` | Unit of energy ['`atomic'`, `'J'`, `'wavenumber'`, `'Hz'`, `'GHz'`] | `String` | No | `'GHz'` |

#### Method: `radiative_lifetime_all`
Radiative lifetime of a state, considering all possible electric dipole allowed transitions.

| Parameter  | Description | Data type | Required | Default |
| ------------ | ------------ | ------------ | ------------ | ------------ |
| `state` | State to calculate the radiative lifetime of. | `State` | Yes | N/A |
| `units` | Unit of time [`'s'`, `'ms'`, `'us'`, `'ns'`] | `String` | No | `'s'` |

#### Method: `transition_energy`
Transition energy between two states.

| Parameter  | Description | Data type | Required | Default |
| ------------ | ------------ | ------------ | ------------ | ------------ |
| `state_1` | State to calculate the transition energy from. | `State` | Yes | N/A |
| `state_2` | State to calculate the transition energy to. | `State` | Yes | N/A |
| `units` | Unit of energy ['`atomic'`, `'J'`, `'wavenumber'`, `'Hz'`, `'GHz'`]  | `String` | No | `'GHz'` |

#### Method: `transition_dipole_moment`
Transition electric dipole moment between two states.

| Parameter  | Description | Data type | Required | Default |
| ------------ | ------------ | ------------ | ------------ | ------------ |
| `state_1` | State to calculate the transition dipole moment from. | `State` | Yes | N/A |
| `state_2` | State to calculate the transition dipole moment to. | `State` | Yes | N/A |
| `units` | Unit of dipole moment [`'atomic'`, `'debye'`, `'coloumb meter'`] | `String` | No | `'debye'` |

#### Method: `radiative_lifetime`
Radiative lifetime between two electric dipole allowed states.

| Parameter  | Description | Data type | Required | Default |
| ------------ | ------------ | ------------ | ------------ | ------------ |
| `state_1` | State to calculate the radiative lifetime from. | `State` | Yes | N/A |
| `state_2` | State to calculate the radiative lifetime to. | `State` | Yes | N/A |
| `units` | Unit of time [`'s'`, `'ms'`, `'us'`, `'ns'`] | `String` | No | `'s'` |


Version information
-------------------

| Library  | Version |
| ------------ | ------------ |
| `Python`  | 3.6.1 64bit [GCC 4.2.1 Compatible Apple LLVM 6.0 (clang-600.0.57)] |
| `IPython` | 5.3.0 |
| `OS` | Darwin 17.4.0 x86_64 i386 64bit |
| `attr` | 17.4.0 |
| `matplotlib` | 2.0.2 |
| `numba` | 0.35.0 |
| `numpy` | 1.14.3 |
| `scipy` | 1.00.0 |
| `sympy` | 1.0 |
| `tqdm` | 4.15.0 |
| `version_information` | 1.0.3 |