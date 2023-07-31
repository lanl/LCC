#!/usr/bin/env python
"""
growth_kmc.py: performs a per-molecule KMC simulation of crystal growth from solution starting with a spherical seed
Usage: python growth_kmc.py input_file.json or ./growth_kmc.py input_file.json
"""

__author__ = 'Jacob Jeffries'
__organization__ = ['Clemson University Department of Materials Science and Engineering',
                    'Theoretical Division, Los Alamos National Laboratory']
__email__ = ['jwjeffr@clemson.edu', 'jwjeffr@lanl.gov']

import numpy as np
import logging
import numba
import multiprocessing as mp
from bisect import bisect_left
from sys import argv
from json import load
from copy import deepcopy
from functools import partial

import lattices
import energetics
from miscellaneous import Timer, norm


__boltzmann_constant__ = 8.617e-5


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
        distance = norm(disp)
        if distance < cutoff:
            ids.append(j + 1)

    return ids


@Timer(callback=logging.info)
def get_all_neighbors(positions: np.ndarray, box_bounds: np.ndarray, max_cutoff: float,
                      num_cpus: int) -> numba.typed.Dict:
    """
    get_all_neighbors returns a dictionary with information about neighbors
    each key is a site identifier, each value is a list of identifiers that neighbor the key
    """

    # get neighbor lists parallelized
        
    func = partial(get_neighbors, positions=positions, box_bounds=box_bounds, cutoff=max_cutoff)
    
    if num_cpus < 1 or type(num_cpus) is not int:
    
        logging.info('WARNING: number of specified cpus needs to be an integer greater than 1. ending simulation')
    
    try:
    
        with mp.Pool(num_cpus) as p:
            neighbor_lists = p.map(func, positions)
            
    except OSError:
    
        new_num_cores = num_cpus - 1
        
        logging.info(f'WARNING: memory error encountered in get_all_neighbors(). '
                     f'returning number of CPUs in this calculation to {new_num_cores:.0f}')

        neighbors_kwargs = {
            'positions': positions,
            'box_bounds': box_bounds,
            'max_cutoff': max_cutoff,
            'num_cpus': new_num_cores
        }
    
        return get_all_neighbors(**neighbors_kwargs)

    # initialize numba dict, needs to be a numba dict so numba can compile it in other functions

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
def compute_evaporation_rates(
        types: np.ndarray,
        solid_sites: np.ndarray,
        initial_energy: float,
        beta: float,
        h_data: tuple,
        prefactor: float
) -> np.ndarray:
    """
    compute_evaporation_rates() calculates all evaporation rates
    """

    # initialize rates array and some prefactor

    rates = np.zeros(solid_sites.shape[0])

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
def compute_adsorption_rates(
        types: np.ndarray,
        solvent_sites: np.ndarray,
        initial_energy: float,
        beta: float,
        prefactor: float,
        barrier: float
) -> np.ndarray:
    """
    compute_adsorption_rates() calculates all adsorption rates
    """

    # initialize rates array and some prefactor

    rates = np.zeros(solvent_sites.shape[0])

    for i in numba.prange(solvent_sites.shape[0]):

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

    header_lines = [
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
    
    for line in header_lines:
        print(line, file=file, flush=True)

    for id_, type_, (x, y, z) in zip(ids, types, positions):
        if type_ == 0:
            continue
        print(f'{id_} {type_} {x} {y} {z}', file=file, flush=True)

    logging.info(f'info at step {step:.0f} and time {t:.2E} dumped')
    
    
__available_lattice_types__ = ['petn_molecular', 'petn_block']
__available_energetics_types__ = ['isotropic_second_nearest', 'anisotropic_third_nearest']


@Timer(callback=logging.info)
def perform_sim(
        box_dimensions: list,
        dump_every: int,
        initial_radius: float,
        dump_file_name: str,
        temperature: float,
        num_steps: int,
        lattice_type: str,
        energetics_type: str,
        evaporation_prefactor: float,
        adsorption_prefactor: float,
        adsorption_barrier: float,
        num_cpus='all',
        log_file_name='kmc.log',
        **kwargs
) -> None:
    """
    perform_sim() performs the simulation
    """
    
    logging.basicConfig(filename=log_file_name, filemode='w', format='%(message)s', level=logging.INFO)

    # calculate thermodynamic beta

    if num_cpus == 'all':
        num_cpus = mp.cpu_count()
        
    elif type(num_cpus) is int:
        numba.set_num_threads(num_cpus)
        
    else:
        logging.info("num_cpus must be type int or 'all'. ending simulation")
        return

    beta = 1.0 / (__boltzmann_constant__ * temperature)

    # initialize types, positions, box bounds, identifiers, and the corresponding hamiltonian

    if lattice_type.lower() == 'petn_molecular':

        lattice_kwargs = {
            'a': kwargs['a'],
            'b': kwargs['b'],
            'c': kwargs['c'],
            'offset': np.array(kwargs['offset']),
            'dimensions': np.array(box_dimensions)
        }

        lattice = lattices.PETNMolecularLattice(**lattice_kwargs)
        
    elif lattice_type.lower() == 'petn_block':
    
        lattice_kwargs = {
            'a': kwargs['a'],
            'b': kwargs['b'],
            'c': kwargs['c'],
            'dimensions': np.array(box_dimensions)
        }
        
        lattice = lattices.PETNBlockLattice(**lattice_kwargs)

    else:

        raise ValueError(f'valid lattice_type inputs are {__available_lattice_types__}')
        
    types, positions, box_bounds, ids = lattice.initialize_simulation()

    if energetics_type.lower() == 'isotropic_second_nearest':

        energetics_kwargs = {
            'first_cutoff': kwargs['first_cutoff'],
            'first_energy': kwargs['first_energy'],
            'second_cutoff': kwargs['second_cutoff'],
            'second_energy': kwargs['second_energy']
        }

        max_cutoff = kwargs['second_cutoff']
        
        energetics_ = energetics.IsotropicSecondNearest(**energetics_kwargs)
        
    elif energetics_type.lower() == 'anisotropic_third_nearest':
    
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
        
        max_cutoff = np.sqrt(kwargs['a'] ** 2 + kwargs['b'] ** 2 + kwargs['c'] ** 2)
        
        energetics_ = energetics.AnisotropicThirdNearest(**energetics_kwargs)

    else:

        raise ValueError(f'valid energetics_type inputs are {__available_energetics_types__}')

    h_data = energetics_.get_hamiltonian(positions, box_bounds)
    logging.info(f'{len(ids):.0f} sites initialized')

    # start with a spherical crystal at the center

    crystal_origin = np.array(box_bounds) / 2
    for index, position in enumerate(positions):
        if np.linalg.norm(position - crystal_origin) <= initial_radius:
            types[index] = 1
    logging.info(f'{np.sum(types == 1):.0f} sites changed to type 1')

    # calculate neighbor dictionary

    neighbor_ids = get_all_neighbors(positions, box_bounds, max_cutoff, num_cpus=num_cpus)
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
            evaporation_args = (
                types,
                solid_sites,
                initial_energy,
                beta,
                h_data,
                evaporation_prefactor
            )
            evaporation_rates = compute_evaporation_rates(*evaporation_args)
            adsorption_args = (
                types,
                solvent_sites,
                initial_energy,
                beta,
                adsorption_prefactor,
                adsorption_barrier
            )
            adsorption_rates = compute_adsorption_rates(*adsorption_args)

            # put rates and sites into one array

            rates = np.concatenate((evaporation_rates, adsorption_rates))
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
            # np.random.uniform() is [0, 1), but KMC calls for (0, 1]
            random_draw = 1.0 - np.random.uniform(low=0, high=1)
            event_index = bisect_left(cumulative_function, random_draw * total_rate)

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

    # uncomment the line below if you are having memory issues
    # mp.set_start_method('forkserver')

    # parse input file
    
    try:
        json_file_name = argv[1]
    except IndexError:
        raise ValueError('needs an input file - see example_input.json')
    
    with open(json_file_name, 'r') as input_file:
        real_simulation_parameters = load(input_file)
    
    # define dummy simulation parameters
    # same simulation with a smaller geometry and different dump information
    
    dummy_simulation_parameters = deepcopy(real_simulation_parameters)
    dummy_simulation_parameters['box_dimensions'] = [10, 10, 10]
    dummy_simulation_parameters['num_steps'] = 10
    dummy_simulation_parameters['dump_every'] = 1
    dummy_simulation_parameters['dump_file_name'] = 'small.dump'
    dummy_simulation_parameters['initial_radius'] = 8.0
    
    # a light simulation is called first, compiling all the expensive function calls

    perform_sim(**dummy_simulation_parameters)
    logging.info('dummy simulation completed')

    # perform expensive simulation next
    
    perform_sim(**real_simulation_parameters)    
    logging.info('real simulation completed')


if __name__ == '__main__':
    main()
