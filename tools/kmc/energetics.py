from dataclasses import dataclass
import numba
import numpy as np
from miscellaneous import is_close, norm, is_parallel, any_


@dataclass
class Energetics:

    def get_hamiltonian(self, positions, box_bounds):

        pass


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

    # nested for loop, need to loop through all neighbor pairs

    for i in np.arange(len(positions)):
        first_position = positions[i]
        for j, second_position in enumerate(positions):

            # avoid calculations if positions are the same

            if i == j:
                continue

            # need to account for periodic boundary conditions before distance is calculated

            dr = first_position - second_position
            dr += -np.rint(dr / box_bounds) * box_bounds
            distance = np.sqrt(np.sum(dr ** 2))

            # append matrix element based on distance, i < j to avoid double-counting

            if distance < first_cutoff and i < j:
                data.append(first_energy)
            elif first_cutoff < distance < second_cutoff and i < j:
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

    def get_hamiltonian(self, positions, box_bounds):

        return _get_hamiltonian_isotropic_second_nearest(
            self.first_cutoff,
            self.first_energy,
            self.second_cutoff,
            self.second_energy,
            positions,
            box_bounds
        )
        
        
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

    distances_second = np.array([norm(a + b), norm(b + c), norm(c + a)])
    distance_third = norm(a + b + c)

    # nested for loop, need to loop through all neighbor pairs

    for i in np.arange(len(positions)):
        first_position = positions[i]
        for j, second_position in enumerate(positions):

            # avoid calculations if positions are the same

            if i == j:
                continue

            # need to account for periodic boundary conditions before distance is calculated

            dr = first_position - second_position
            dr += -np.rint(dr / box_bounds) * box_bounds
            distance = norm(dr)

            # append matrix element based on distance, i < j to avoid double-counting

            if is_parallel(dr, a) and i < j:
                data.append(e_1a)
            elif is_parallel(dr, b) and i < j:
                data.append(e_1b)
            elif is_parallel(dr, c) and i < j:
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

    def get_hamiltonian(self, positions, box_bounds):

        return _get_hamiltonian_anisotropic_third_nearest_reconstruction(
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

    # nested for loop, need to loop through all neighbor pairs

    for i in np.arange(len(positions)):
        first_position = positions[i]
        for j, second_position in enumerate(positions):

            # avoid calculations if positions are the same

            if i == j:
                continue

            # need to account for periodic boundary conditions before distance is calculated

            dr = first_position - second_position
            dr += -np.rint(dr / box_bounds) * box_bounds

            # append matrix element based on distance, i < j to avoid double-counting

            if is_parallel(dr, a) and i < j:
                data.append(e_1a)
            elif is_parallel(dr, b) and i < j:
                data.append(e_1b)
            elif is_parallel(dr, c) and i < j:
                data.append(e_1c)
            elif is_parallel(dr, b + c) and i < j:
                data.append(e_2a)
            elif is_parallel(dr, b - c) and i < j:
                data.append(e_2a_p)
            elif is_parallel(dr, c + a) and i < j:
                data.append(e_2b)
            elif is_parallel(dr, -c + a) and i < j:
                data.append(e_2b_p)
            elif is_parallel(dr, a + b) and i < j:
                data.append(e_2c)
            elif is_parallel(dr, a - b) and i < j:
                data.append(e_2c_p)
            elif is_parallel(dr, a + b + c) and i < j:
                data.append(e_31)
            elif is_parallel(dr, a - b - c) and i < j:
                data.append(e_32)
            elif is_parallel(dr, a - b + c) and i < j:
                data.append(e_33)
            elif is_parallel(dr, a + b - c) and i < j:
                data.append(e_34)
            else:
                continue

            row.append(i)
            col.append(j)

    return np.array(data), np.array(row), np.array(col)


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

    def get_hamiltonian(self, positions, box_bounds):

        return _get_hamiltonian_anisotropic_third_nearest(
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
