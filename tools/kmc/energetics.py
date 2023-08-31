# module dealing with energetics calculations, mainly hamiltonian initialization

from dataclasses import dataclass
import numba
import numpy as np
import pickle
from miscellaneous import is_close, norm, is_parallel, any_, Deprecated


__available_energetics_types__ = [
    'isotropic_second_nearest',
    'anisotropic_third_nearest',
    'anisotropic_third_nearest_reconstruction',
    'pickled_energetics'
]


@dataclass
class Energetics:

    def get_hamiltonian(self, positions, box_bounds):

        pass
        
    @staticmethod
    def _print_hamiltonian(data, row, col, positions, box_bounds, file_name):
    
        num_sites = len(positions)
        num_non_zero = len(data)
        
        prop = len(data) / len(positions) ** 2
    
        with open(file_name, 'w') as file:
        
            print(f'shape = ({num_sites:.0f}, {num_sites:.0f})\n', file=file)
            print(f'proportion of non-zero elements = {prop * 100}%\n', file=file)
            print('\n#########################################\n', file=file)
        
            print(f'(first index, second index): matrix element, difference vector\n', file=file)
            print('-----------------------------------------\n', file=file)
            for I, J, data in zip(row, col, data):
                
                difference = positions[I] - positions[J]
                
                # account for pbc
                
                difference += -np.rint(difference / box_bounds) * box_bounds
                print(f'({I:.0f}, {J:.0f}): {data}, {difference}\n', file=file)


@numba.njit(cache=True)
def _get_hamiltonian_isotropic_second_nearest(
        first_cutoff,
        first_energy,
        second_cutoff,
        second_energy,
        positions,
        box_bounds
):

    """
    _get_hamiltonian_isotropic_second_nearest returns the hamiltonian matrix in key-value format
    for the isotropic second-nearest energetics style
    values are stored in data
    row indices stored in row
    column indices stored in col
    """

    # define empty lists of indices and matrix elements
    # lists are initialized with integers and floats so numba knows which type to initialize the lists with

    row = [0 for _ in range(0)]
    col = [0 for _ in range(0)]
    data = [0.0 for _ in range(0)]
    lengths = [np.abs(bound) for bound in box_bounds]

    # nested for loop, need to loop through all neighbor pairs

    for i in np.arange(len(positions)):
        first_position = positions[i]
        for j, second_position in enumerate(positions):

            # avoid calculations if positions are the same, and double counting (if i > j)

            if i >= j:
                continue

            # need to account for periodic boundary conditions before distance is calculated

            dr = first_position - second_position
            dr += -np.rint(dr / box_bounds) * box_bounds
            distance = norm(dr)

            # append matrix element based on distance, i < j to avoid double-counting

            if distance < first_cutoff:
                data.append(first_energy)
            elif first_cutoff < distance < second_cutoff:
                data.append(second_energy)
            else:
                continue
                
            row.append(i)
            col.append(j)

    return np.array(data), np.array(row), np.array(col)


@dataclass
class IsotropicSecondNearest(Energetics):

    first_cutoff: float
    first_energy: float
    second_cutoff: float
    second_energy: float
    print_hamiltonian: str = None

    def get_hamiltonian(self, positions, box_bounds):
    
        h_data = _get_hamiltonian_isotropic_second_nearest(
            self.first_cutoff,
            self.first_energy,
            self.second_cutoff,
            self.second_energy,
            positions,
            box_bounds
        )
        
        if self.print_hamiltonian is not None:
        
            self._print_hamiltonian(*h_data, positions, box_bounds, self.print_hamiltonian)

        return self.second_cutoff, h_data
        
        
@numba.njit(cache=True)
def _get_hamiltonian_anisotropic_third_nearest_reconstruction(
        e_1a,
        e_1b,
        e_1c,
        second_nearest,
        third_nearest,
        a,
        b,
        c,
        positions,
        box_bounds
):

    """
    _get_hamiltonian_anisotropic_third_nearest_reconstruction returns the hamiltonian matrix in key-value format
    for the anisotropic third-nearest energetics with reconstruction style
    values are stored in data
    row indices stored in row
    column indices stored in col
    """

    # define empty lists of indices and matrix elements
    # lists are initialized with integers and floats so numba knows which type to initialize the lists with

    row = [0 for _ in range(0)]
    col = [0 for _ in range(0)]
    data = [0.0 for _ in range(0)]
    lengths = [np.abs(bound) for bound in box_bounds]

    distances_second = np.array([norm(a + b), norm(b + c), norm(c + a)])
    distance_third = norm(a + b + c)

    # nested for loop, need to loop through all neighbor pairs

    for i in np.arange(len(positions)):
        first_position = positions[i]
        for j, second_position in enumerate(positions):

            # avoid calculations if positions are the same, and double counting (if i > j)

            if i >= j:
                continue

            # need to account for periodic boundary conditions before distance is calculated
            dr = first_position - second_position
            
            for k, dx in enumerate(dr):
                
                if np.abs(dx) > lengths[k] / 2:
                    dx = lengths[k] - dx
                    dr[k] = dx 

            #dr += -np.rint(dr / box_bounds) * box_bounds
            distance = norm(dr)

            # append matrix element based on distance, i < j to avoid double-counting

            if is_parallel(dr, a) and is_close(distance, norm(a)):
                data.append(e_1a)
            elif is_parallel(dr, b) and is_close(distance, norm(b)):
                data.append(e_1b)
            elif is_parallel(dr, c) and is_close(distance, norm(c)):
                data.append(e_1c)
            elif any_([is_close(distance, x) for x in distances_second]):
                data.append(second_nearest)
            elif is_close(distance, distance_third):
                data.append(third_nearest)
            else:
                continue
                
            row.append(i)
            col.append(j)

    return np.array(data), np.array(row), np.array(col)


@Deprecated(reason='WARNING: __name__ deprecated, brick model no longer supported')
@dataclass
class AnisotropicThirdNearestReconstruction(Energetics):

    e_1a: float
    e_1b: float
    e_1c: float
    second_nearest: float
    third_nearest: float
    a: np.ndarray
    b: np.ndarray
    c: np.ndarray
    print_hamiltonian: str = None

    def get_hamiltonian(self, positions, box_bounds):
    
        relative_tolerance = 0.05

#        max_cutoff = (1.0 + relative_tolerance) * np.sqrt(kwargs['a'] ** 2 + kwargs['b'] ** 2 + kwargs['c'] ** 2)
#        possible_cutoffs = [
#            norm(self.a + self.b),
#            norm(self.b + self.c),
#            norm(self.c + self.a)
#        ]
        
#        possible_cutoffs = [(1.0 + relative_tolerance) * x for x in possible_cutoffs]
#        max_cutoff = max(possible_cutoffs)
        max_cutoff = (1.0 + relative_tolerance) * norm(self.a + self.b + self.c)

        h_data = _get_hamiltonian_anisotropic_third_nearest_reconstruction(
            self.e_1a,
            self.e_1b,
            self.e_1c,
            self.second_nearest,
            self.third_nearest,
            self.a,
            self.b,
            self.c,
            positions,
            box_bounds
        )
        
        if self.print_hamiltonian is not None:
        
            self._print_hamiltonian(*h_data, positions, box_bounds, self.print_hamiltonian)
            
        return max_cutoff, h_data


@numba.njit(cache=True)
def _get_hamiltonian_anisotropic_third_nearest(
        e_1a,
        e_1b,
        e_1c,
        e_2a,
        e_2a_p,
        e_2b,
        e_2b_p,
        e_2c,
        e_2c_p,
        e_31,
        e_32,
        e_33,
        e_34,
        a,
        b,
        c,
        positions,
        box_bounds
):
    """
    _get_hamiltonian_anisotropic_third_nearest returns the hamiltonian matrix in key-value format
    for the anisotropic third-nearest energetics style
    values are stored in data
    row indices stored in row
    column indices stored in col
    """

    # define empty lists of indices and matrix elements
    # lists are initialized with integers and floats so numba knows which type to initialize the lists with

    row = [0 for _ in range(0)]
    col = [0 for _ in range(0)]
    data = [0.0 for _ in range(0)]
    lengths = [np.abs(bound) for bound in box_bounds]

    # nested for loop, need to loop through all neighbor pairs

    for i in np.arange(len(positions)):
        first_position = positions[i]
        for j, second_position in enumerate(positions):

            # avoid calculations if positions are the same, and double counting (if i > j)

            if i >= j:
                continue

            # need to account for periodic boundary conditions before distance is calculated

            dr = first_position - second_position
            
            for k, dx in enumerate(dr):
                
                if np.abs(dx) > lengths[k] / 2:
                    dx = lengths[k] - dx
                    dr[k] = dx 
                    
            distance = norm(dr)

            # append matrix element based on distance, i < j to avoid double-counting

            if is_parallel(dr, a) and is_close(distance, norm(a)):
                data.append(e_1a)
            elif is_parallel(dr, b) and is_close(distance, norm(b)):
                data.append(e_1b)
            elif is_parallel(dr, c) and is_close(distance, norm(c)):
                data.append(e_1c)
            elif is_parallel(dr, b + c) and is_close(distance, norm(b + c)):
                data.append(e_2a)
            elif is_parallel(dr, b - c) and is_close(distance, norm(b - c)):
                data.append(e_2a_p)
            elif is_parallel(dr, c + a) and is_close(distance, norm(c + a)):
                data.append(e_2b)
            elif is_parallel(dr, -c + a) and is_close(distance, norm(-c + a)):
                data.append(e_2b_p)
            elif is_parallel(dr, a + b) and is_close(distance, norm(a + b)):
                data.append(e_2c)
            elif is_parallel(dr, a - b) and is_close(distance, norm(a - b)):
                data.append(e_2c_p)
            elif is_parallel(dr, a + b + c) and is_close(distance, norm(a + b + c)):
                data.append(e_31)
            elif is_parallel(dr, a - b - c) and is_close(distance, norm(a - b - c)):
                data.append(e_32)
            elif is_parallel(dr, a - b + c) and is_close(distance, norm(a - b + c)):
                data.append(e_33)
            elif is_parallel(dr, a + b - c) and is_close(distance, norm(a + b - c)):
                data.append(e_34)
            else:
                continue

            row.append(i)
            col.append(j)

    return np.array(data), np.array(row), np.array(col)


@Deprecated(reason='WARNING: __name__ deprecated, brick model no longer supported')
@dataclass
class AnisotropicThirdNearest(Energetics):

    e_1a: float
    e_1b: float
    e_1c: float
    e_2a: float
    e_2a_p: float
    e_2b: float
    e_2b_p: float
    e_2c: float
    e_2c_p: float
    e_31: float
    e_32: float
    e_33: float
    e_34: float
    a: np.ndarray
    b: np.ndarray
    c: np.ndarray
    print_hamiltonian: str = None

    def get_hamiltonian(self, positions, box_bounds):
    
        relative_tolerance = 0.05

#        max_cutoff = (1.0 + relative_tolerance) * np.sqrt(kwargs['a'] ** 2 + kwargs['b'] ** 2 + kwargs['c'] ** 2)
#        possible_cutoffs = [
#            norm(self.a + self.b),
#            norm(self.b + self.c),
#            norm(self.c + self.a)
#        ]
        
#        possible_cutoffs = [(1.0 + relative_tolerance) * x for x in possible_cutoffs]
#        max_cutoff = max(possible_cutoffs)
        max_cutoff = (1.0 + relative_tolerance) * norm(self.a + self.b + self.c)

        h_data = _get_hamiltonian_anisotropic_third_nearest(
            self.e_1a,
            self.e_1b,
            self.e_1c,
            self.e_2a,
            self.e_2a_p,
            self.e_2b,
            self.e_2b_p,
            self.e_2c,
            self.e_2c_p,
            self.e_31,
            self.e_32,
            self.e_33,
            self.e_34,
            self.a,
            self.b,
            self.c,
            positions,
            box_bounds
        )
        
        if self.print_hamiltonian:
            self._print_hamiltonian(*h_data, positions, box_bounds, self.print_hamiltonian)
            
        return max_cutoff, h_data
        
        
        
@dataclass
class PickledEnergetics(Energetics):

    energetics_loading_file: str
    
    def get_hamiltonian(self, positions, box_bounds):
    
        with open(energetics_loading_file, 'rb') as e_file:
            max_cutoff, h_data = pickle.load(e_file)
            
        if self.print_hamiltonian:
            self._print_hamiltonian(*h_data, positions, box_bounds, self.print_hamiltonian)
            
        return max_cutoff, h_data
        
        
def get_energetics(energetics_type: str, **kwargs) -> Energetics:

    if energetics_type.lower() == 'isotropic_second_nearest':

        energetics_kwargs = {
            'first_cutoff': kwargs['first_cutoff'],
            'first_energy': kwargs['first_energy'],
            'second_cutoff': kwargs['second_cutoff'],
            'second_energy': kwargs['second_energy']
        }
        
        energetics_ = IsotropicSecondNearest(**energetics_kwargs)

    elif energetics_type.lower() == 'anisotropic_third_nearest':

        energetics_kwargs = {
            'e_1a': kwargs['e_1a'],
            'e_1b': kwargs['e_1b'],
            'e_1c': kwargs['e_1c'],
            'e_2a': kwargs['e_2a'],
            'e_2a_p': kwargs['e_2a_p'],
            'e_2b': kwargs['e_2b'],
            'e_2b_p': kwargs['e_2b_p'],
            'e_2c': kwargs['e_2c'],
            'e_2c_p': kwargs['e_2c_p'],
            'e_31': kwargs['e_31'],
            'e_32': kwargs['e_32'],
            'e_33': kwargs['e_33'],
            'e_34': kwargs['e_34'],
            'a': kwargs['a'] * np.array([1, 0, 0]),
            'b': kwargs['b'] * np.array([0, 1, 0]),
            'c': kwargs['c'] * np.array([0, 0, 1])
        }
        
        energetics_ = AnisotropicThirdNearest(**energetics_kwargs)
        
    elif energetics_type.lower() == 'anisotropic_third_nearest_reconstruction':
    
        energetics_kwargs = {
            'e_1a': kwargs['e_1a'],
            'e_1b': kwargs['e_1b'],
            'e_1c': kwargs['e_1c'],
            'second_nearest': kwargs['second_nearest'],
            'third_nearest': kwargs['third_nearest'],
            'a': kwargs['a'] * np.array([1, 0, 0]),
            'b': kwargs['b'] * np.array([0, 1, 0]),
            'c': kwargs['c'] * np.array([0, 0, 1])
        }
        
        energetics_ = AnisotropicThirdNearestReconstruction(**energetics_kwargs)
        
    elif energetics_type.lower() == 'pickled_energetics': 
    
        energetics_kwargs = {
            'energetics_loading_file': kwargs['energetics_loading_file']
        }
        
        energetics_ = PickledEnergetics(energetics_loading_file=kwargs['energetics_loading_file'])

    else:

        raise ValueError(f'valid energetics_type inputs are {__available_energetics_types__}')
        
    try:
        energetics_.print_hamiltonian = kwargs['print_hamiltonian']
    except KeyError:
        pass
        
    return energetics_
    
    
@numba.njit(parallel=True, cache=True)
def calculate_total_energy(data: np.ndarray, row: np.ndarray, col: np.ndarray, types: np.ndarray) -> float:
    """
    calculate_total_energy() calculates the total energy from the non-zero keys of the hamiltonian
    """

    assert data.shape == row.shape == col.shape
    
    return np.sum(data * types[row] * types[col])
