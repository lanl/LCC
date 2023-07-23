#!/usr/bin/env python
"""
growth_kmc.py: performs a per-molecule KMC simulation of crystal growth from solution starting with a spherical seed
Usage: python growth_kmc.py input_file.json or ./growth_kmc.py input_file.json
"""

__author__ = 'Jacob Jeffries'
__organization__ = ['Clemson University Department of Materials Science and Engineering', 'Theoretical Division, Los Alamos National Laboratory']
__email__ = ['jwjeffr@clemson.edu', 'jwjeffr@lanl.gov']


import numpy as np
from dataclasses import dataclass
import time
from itertools import product
import logging
import numba
import multiprocessing as mp
from bisect import bisect_right
import sys
import json


# define some parameters
__log_file_name__ = 'kmc.log'
__boltzmann_constant__ = 8.617e-5

# define dummy simulation parameters
# a light simulation is called first, compiling all the expensive function calls
__dummy_simulation_parameters__ = {
    'box_dimensions': [10, 10, 10],
    'num_steps': 10,
    'dump_every': 1,
    'dump_file_name': 'small.dump',
    'initial_radius': 8.0,
    'temperature': 300.0,
    'first_cutoff': 7.0,
    'second_cutoff': 7.5,
    'first_energy': -0.291,
    'second_energy': -0.186
}

try:
    json_file_name = sys.argv[1]
except IndexError:
    raise ValueError('needs an input file - see example_input.json')

with open(sys.argv[1], 'r') as file:
    __real_simulation_parameters__ = json.load(file)


logging.basicConfig(filename=__log_file_name__, filemode='w', format='%(message)s', level=logging.INFO)


@dataclass
class Timer:

    """
    Timer class - can be used as a decorator that times functions
    callback keyword argument dictates where callback statement is recorded
    """

    callback: callable = print

    def __call__(self, function: callable) -> callable:

        def wrapper_func(*args, **kwargs):
            past = time.time()
            result = function(*args, **kwargs)
            present = time.time()
            if self.callback is not None:
                self.callback(f'Function {function.__name__} executed in {present - past} sec')
            return result

        return wrapper_func


@Timer(callback=logging.info)
@numba.njit()
def get_hamiltonian(
        positions: np.ndarray,
        box_bounds: np.ndarray,
        first_cutoff: float,
        second_cutoff: float,
        first_energy: float,
        second_energy: float
) -> tuple:

    """
    get_hamiltonian() returns the hamiltonian matrix in key-value format
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
            distance = np.linalg.norm(dr)

            # append matrix element based on distance, i < j to avoid double-counting

            if distance < first_cutoff and i < j:
                row.append(i)
                col.append(j)
                data.append(first_energy)
            elif first_cutoff < distance < second_cutoff and i < j:
                row.append(i)
                col.append(j)
                data.append(second_energy)
                
    return np.array(data), np.array(row), np.array(col)
    
    
def initialize_simulation(box_dimensions: list) -> tuple:

    """
    initialize_simulation() returns arrays of types, positions, box bounds, and site identifiers
    input is box_dimensions, which is a tuple with (i_max + 1, j_max + 1, k_max + 1)
    """

    # initialize lattice vectors

    a = np.array([9.088, 0.0, 0.0])
    b = np.array([0.0, 9.088, 0.0])
    c = np.array([0.0, 0.0, 6.737])

    # create particle labels
    # convoluted way, but all that matters is that we have 2 * box_dims[0] * box_dims[1] * box_dims[2] unique labels
    # along with corresponding types (initialized to 0) and positions

    particle_labels = product(*(*(np.arange(d) for d in box_dimensions), np.arange(2)))
    types = np.zeros(np.product(box_dimensions) * 2, dtype=np.int_)
    positions = np.zeros((np.product(box_dimensions) * 2, 3))

    # populate positions with i * a + j * b + k * c + l * dr

    offset = np.array([4.54348, 4.54346, 3.36908])
    for particle_index, (i, j, k, l) in enumerate(particle_labels):

        positions[particle_index] = i * a + j * b + c * k + offset * l

    # give each site a unique identifier, starting at 1
    num_sites = len(positions)
    ids = np.arange(1, num_sites + 1)

    # simulation box bounds are extrema of the initialized position array
    box_bounds = np.max(positions, axis=0)
    
    return types, positions, box_bounds, ids
    
    
@numba.njit()
def get_neighbors(position: np.ndarray, positions: np.ndarray, box_bounds: tuple, cutoff: float) -> list:

    """
    get_neighbors() returns a list of identifiers which neighbor a specific position
    """

    # get displacement vectors

    displacement_vectors = [position - x for x in positions]

    ids = []

    for j, disp in enumerate(displacement_vectors):

        # account for periodic boundary conditions
        
        disp += -np.rint(disp / box_bounds) * box_bounds
        distance = np.linalg.norm(disp)
        if distance < cutoff:
            ids.append(j + 1)

    return ids
    
    
@Timer(callback=logging.info)
def get_all_neighbors(positions: np.ndarray, box_bounds: np.ndarray, second_cutoff: float, num_cpus: int = None) -> numba.typed.Dict:

    """
    get_all_neighbors returns a dictionary with information about neighbors
    each key is a site identifier, each value is a list of identifiers that neighbor the key
    """

    # default to using all but one cpu, in case there's a master process or something

    if num_cpus is None:
        num_cpus = mp.cpu_count()

    # get neighbor lists parallelized

    args = [(position, positions, box_bounds, second_cutoff) for position in positions]

    with mp.Pool(num_cpus) as p:
        neighbor_lists = p.starmap(get_neighbors, args)

    # initialize numba dict, needs to be a numba dict to numba can compile it in other functions

    neighbor_ids = numba.typed.Dict.empty(key_type=numba.types.int64, value_type=numba.types.int64[:])

    for index, _ in enumerate(positions):
        ids = neighbor_lists[index]
        neighbor_ids[index + 1] = np.array(ids, dtype=np.int_)

    return neighbor_ids
    
    
@numba.njit(parallel=True)
def calculate_total_energy(data: np.ndarray, row: np.ndarray, col: np.ndarray, types: np.ndarray) -> float:

    """
    calculate_total_energy() calculates the total energy from the non-zero keys of the hamiltonian
    """

    assert data.shape == row.shape == col.shape

    # typical parallelized summation
    
    array_to_sum = np.zeros(data.shape[0])
    
    for i in numba.prange(len(data)):
    
        array_to_sum[i] = data[i] * types[row[i]] * types[col[i]]
        
    return np.sum(array_to_sum)
    
    
@Timer(callback=print)
@numba.njit(parallel=True)
def get_surface_sites(ids: np.ndarray, types: np.ndarray, neighbor_ids: dict) -> tuple:

    """
    get_surface_sites() calculates the identifiers of sites at the interface
    returns both solvent and solid identifiers
    """

    # initialize arrays of identifiers, shrink this down after calculation

    solvent_ids = np.zeros(ids.shape, dtype=np.int_)
    solid_ids = np.zeros(ids.shape, dtype=np.int_)

    # populate arrays with identifiers if two sites that neighbor each other have different types
    
    for i in numba.prange(len(ids)):
    
        id_, type_ = ids[i], types[i]
        neighbor_list = neighbor_ids[id_]
        for neighbor_id in neighbor_list:
            neighbor_index = neighbor_id - 1
            second_type = types[neighbor_index]
            if type_ == second_type:
                continue
            if type_ == 1 and second_type == 0:
                solid_id = id_
                solvent_id = neighbor_id
            elif type_ == 0 and second_type == 1:
                solid_id = neighbor_id
                solvent_id = id_
            else:
                raise ValueError

            solvent_ids[i] = solvent_id
            solid_ids[i] = solid_id

    # turn into sets to get unique IDs
    solvent_ids = list(set(solvent_ids))
    solid_ids = list(set(solid_ids))
    
    # need to exclude the 0

    return np.array(solvent_ids[1:], dtype=np.int_), np.array(solid_ids[1:], dtype=np.int_)
    
    
@numba.njit(parallel=True)
def compute_evap_rates(
        types: np.ndarray,
        solid_sites: np.ndarray,
        initial_energy: float,
        beta: float,
        h_data: tuple
) -> np.ndarray:

    """
    compute_evap_rates() calculates all evaporation rates
    """

    # initialize rates array and some prefactor

    rates = np.zeros(solid_sites.shape[0])
    prefactor = 1e+10

    for i in numba.prange(solid_sites.shape[0]):

        # create an array of new types, swap the site of interest to a zero, calculate change in energy

        site_id = solid_sites[i]
        new_types = types.copy()
        site_index = site_id - 1
        new_types[site_index] = 0
        final_energy = calculate_total_energy(*h_data, new_types)

        barrier = final_energy - initial_energy

        # store arrhenius rate

        rates[i] = prefactor * np.exp(-beta * barrier)

    return rates


@numba.njit(parallel=True)
def compute_adsorp_rates(
        types: np.ndarray,
        solvent_sites: np.ndarray,
        initial_energy: float,
        beta: float
) -> np.ndarray:

    """
    compute_evap_rates() calculates all evaporation rates
    """

    # initialize rates array and some prefactor

    rates = np.zeros(solvent_sites.shape[0])
    prefactor = 1e+10

    for i in numba.prange(solvent_sites.shape[0]):

        # arbitrary barrier, should be solvent barrier diffusion

        barrier = 0.9

        # store arrhenius rate

        rates[i] = prefactor * np.exp(-beta * barrier)

    return rates
    
    
def dump_info(
        step: int,
        t: float,
        types: np.ndarray,
        box_bounds: np.ndarray,
        ids: np.ndarray,
        positions: np.ndarray,
        file
) -> None:

    """
    dump_info() dumps the info at a given step
    only dumps occupied sites
    """

    lines = [
        'ITEM: TIMESTEP',
        f'{step:.0f} {t}',
        'ITEM: NUMBER OF ATOMS',
        f'{np.sum(types == 1):.0f}',
        'ITEM: BOX BOUNDS',
        f'0 {box_bounds[0]}',
        f'0 {box_bounds[1]}',
        f'0 {box_bounds[2]}',
        'ITEM: ATOMS id type x y z'
    ]

    for id_, type_, (x, y, z) in zip(ids, types, positions):
        if type_ == 0:
            continue
        line = f'{id_} {type_} {x} {y} {z}'
        lines.append(line)

    for line in lines:
        print(line, file=file, flush=True)

    logging.info(f'info at step {step:.0f} and time {t:.2E} dumped')
    
    
@Timer(callback=logging.info)
def perform_sim(
        box_dimensions: list,
        dump_every: int,
        initial_radius: float,
        dump_file_name: str,
        temperature: float,
        num_steps: int,
        first_cutoff: float,
        second_cutoff: float,
        first_energy: float,
        second_energy: float
) -> None:

    """
    perform_sim() performs the simulation
    """

    # calculate thermodynamic beta

    beta = 1.0 / (__boltzmann_constant__ * temperature)

    # initialize types, positions, box bounds, identifiers, and the corresponding hamiltonian
    
    types, positions, box_bounds, ids = initialize_simulation(box_dimensions)
    h_data = get_hamiltonian(positions, box_bounds, first_cutoff, second_cutoff, first_energy, second_energy)
    logging.info(f'{len(ids):.0f} sites initialized')

    # start with a spherical crystal at the center

    crystal_origin = np.array(box_bounds) / 2
    for index, position in enumerate(positions):
        if np.linalg.norm(position - crystal_origin) <= initial_radius:
            types[index] = 1
    logging.info(f'{np.sum(types == 1):.0f} sites changed to type 1')

    # calculate neighbor dictionary

    neighbor_ids = get_all_neighbors(positions, box_bounds, second_cutoff)
    logging.info('neighbors calculated')

    # begin simulation

    t = 0

    with open(dump_file_name, 'w') as file:

        for step in range(num_steps):

            # kill the simulation if there's no longer an interface
            # equivalently, if crystal either outgrows the box or completely evaporates
        
            if np.sum(types == 1) == 0 or np.sum(types == 0) == 0:
                logging.info('WARNING: no interface detected. ending simulation')
                break

            # log info

            if step % dump_every == 0:

                dump_info(step, t, types, box_bounds, ids, positions, file)

            # recalculate energy, active sites, and possible rates

            initial_energy = calculate_total_energy(*h_data, types)
            solvent_sites, solid_sites = get_surface_sites(ids, types, neighbor_ids)
            evap_rates = compute_evap_rates(types, solid_sites, initial_energy, beta, h_data)
            adsorp_rates = compute_adsorp_rates(types, solvent_sites, initial_energy, beta)

            # put rates and sites into one array

            rates = np.concatenate((evap_rates, adsorp_rates))
            surface_sites = np.concatenate((solid_sites, solvent_sites))
            if step % dump_every == 0:
                logging.info(f'min rate = {min(rates):.2E}')
                logging.info(f'max rate = {max(rates):.2E}')

            # sort rates in ascending order, sort sites according to how rates were sorted

            sorted_indices = rates.argsort()
            sorted_rates, sorted_sites = rates[sorted_indices], surface_sites[sorted_indices]
            cumulative_function = np.cumsum(sorted_rates)
            total_rate = cumulative_function[-1]

            # pick an event according to KMC algorithm
            # bisect_right performs binary search
            random_draw = np.random.uniform(low=0, high=1)
            event_index = bisect_right(cumulative_function, random_draw * total_rate)

            # update timestep

            second_draw = np.random.uniform(low=0, high=1)
            dt = np.log(1.0 / second_draw) * 1.0 / total_rate
            t += dt

            # swap chosen site

            site_index = sorted_sites[event_index] - 1
            types[site_index] = int(not types[site_index])

        # dump final info

        dump_info(step, t, types, box_bounds, ids, positions, file)
                
    logging.info('run finished')


def main():

    # first perform small simulation

    perform_sim(**__dummy_simulation_parameters__)
    logging.info('dummy simulation completed')

    # perform expensive simulation next

    perform_sim(**__real_simulation_parameters__)
    
    
if __name__ == '__main__':

    main()
