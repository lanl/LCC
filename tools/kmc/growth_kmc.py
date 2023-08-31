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
import sys
import json
import copy
from functools import partial
import time
import pickle

import lattices
import energetics
import configurations
from miscellaneous import Timer, norm, dump_info


__boltzmann_constant__ = 8.617e-5


@numba.njit(cache=True)
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
    

def _get_all_neighbors(positions: np.ndarray, box_bounds: np.ndarray, max_cutoff: float,
                      num_cpus: int) -> numba.typed.Dict:
    """
    _get_all_neighbors returns a dictionary with information about neighbors
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
        
        logging.info(f'WARNING: OSError encountered in _get_all_neighbors(). likely a memory error. '
                     f'changing number of CPUs in this calculation to {new_num_cores:.0f}')
    
        return _get_all_neighbors(positions, box_bounds, max_cutoff, new_num_cores)

    # initialize numba dict, needs to be a numba dict so numba can compile it in other functions

    neighbor_ids = numba.typed.Dict.empty(key_type=numba.types.int64, value_type=numba.types.int64[:])

    for index, _ in enumerate(positions):
        ids = neighbor_lists[index]
        neighbor_ids[index + 1] = np.array(ids, dtype=np.int_)

    return neighbor_ids


@Timer(callback=logging.info)
def get_all_neighbors(*args, **kwargs) -> numba.typed.Dict:
    """
    get_all_neighbors is a wrapper for _get_all_neighbors
    makes it so neighbor calculation isn't wrapped by Timer() multiple times
    """

    return _get_all_neighbors(*args, **kwargs)


@numba.njit(parallel=True, cache=True)
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
    

@numba.njit(cache=True)
def compute_energy_change(
        flipped_bit: int,
        initial_state: np.ndarray,
        final_state: np.ndarray,
        diag_term: float,
        q_vec: np.ndarray
) -> float:

    self_term = (initial_state[flipped_bit] - final_state[flipped_bit]) * diag_term
    interaction_term = np.dot(final_state[flipped_bit] * final_state - initial_state[flipped_bit] * initial_state, q_vec)

    return self_term + interaction_term
    
    
@numba.njit(parallel=True, cache=True)
def compute_evaporation_rates(
        types: np.ndarray,
        solid_sites: np.ndarray,
        beta: float,
        prefactor: float,
        hamiltonian: np.ndarray
) -> np.ndarray:

    rates = np.zeros(solid_sites.shape[0])
    
    for i in numba.prange(solid_sites.shape[0]):
        
        site_id = solid_sites[i]
        new_types = types.copy()
        site_index = site_id - 1
        new_types[site_index] = 0
        diag_term = hamiltonian[site_index, site_index]
        q_vec = hamiltonian[site_index, :] + hamiltonian[:, site_index]
        barrier = compute_energy_change(site_index, types, new_types, diag_term, q_vec)
        
        rates[i] = prefactor * np.exp(-beta * barrier)
        
    return rates


@numba.njit(parallel=True, cache=True)
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

    # uncomment below if we care about site-dependent adsorption
    
    """

    rates = np.zeros(solvent_sites.shape[0])

    for i in numba.prange(solvent_sites.shape[0]):

        # store arrhenius rate

        rates[i] = prefactor * np.exp(-beta * barrier)

    return rates
    """
    
    return np.ones(solvent_sites.shape[0]) * prefactor * np.exp(-beta * barrier)


@Timer(callback=logging.info)
def perform_sim(
        box_dimensions: list,
        dump_every: int,
        dump_file_name: str,
        temperature: float,
        num_steps: int,
        lattice_type: str,
        energetics_type: str,
        init_configuration_type: str,
        evaporation_prefactor: float,
        adsorption_prefactor: float,
        adsorption_barrier: float,
        num_cpus='all',
        log_file_name='kmc.log',
        calculate_surface_every: int = 1,
        lattice_pickle: str = None,
        energetics_pickle: str = None,
        print_hamiltonian: str = None,
        **kwargs
) -> None:
    """
    perform_sim() performs the simulation
    """
    
    logging.basicConfig(filename=log_file_name, filemode='w', format='%(message)s', level=logging.INFO)

    if num_cpus == 'all':
        num_cpus = mp.cpu_count()
        
    elif type(num_cpus) is int:
        numba.set_num_threads(num_cpus)
        
    else:
        logging.info("num_cpus must be type int or 'all'. ending simulation")
        return

    beta = 1.0 / (__boltzmann_constant__ * temperature)

    # initialize types, positions, box bounds, identifiers, and the corresponding hamiltonian
    
    lattice = lattices.get_lattice(lattice_type=lattice_type.lower(), box_dimensions=box_dimensions, **kwargs)
    lattice_info = lattice.initialize_simulation()
    if lattice_pickle:
        with open(lattice_pickle, 'wb') as l_file:
            pickle.dump(lattice_info, l_file)
    positions, box_bounds, ids = lattice_info
    
    energetics_ = energetics.get_energetics(energetics_type=energetics_type.lower(), **kwargs)
    max_cutoff, h_data = energetics_.get_hamiltonian(positions, box_bounds)
    
    ### get array from h_data
    hamiltonian = np.zeros((len(ids), len(ids)))
    
    for element, i, j in zip(*h_data):
        hamiltonian[i, j] = element
    
    logging.info(f'{len(ids):.0f} sites initialized')
    data, row, col = h_data
    if energetics_pickle:
        with open(energetics_pickle, 'wb') as e_file:
            pickle.dump(h_data, e_file)   
    
    configuration = configurations.get_configuration(init_configuration_type=init_configuration_type.lower(), **kwargs)
    types = configuration.set_types(positions, box_bounds)
    logging.info(f'{np.sum(types == 1):.0f} sites changed to type 1')

    # calculate neighbor dictionary
    
    logging.info(f'neighbor cutoff distance for events = {max_cutoff:.2E}')
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

            initial_energy = energetics.calculate_total_energy(*h_data, types)
            
            if step % calculate_surface_every == 0:
            
                solvent_sites, solid_sites = get_surface_sites(ids, types, neighbor_ids)
            
            past = time.time()
            evaporation_rates = compute_evaporation_rates(
                types,
                solid_sites,
                beta,
                evaporation_prefactor,
                hamiltonian
            )
            adsorption_rates = compute_adsorption_rates(
                types,
                solvent_sites,
                initial_energy,
                beta,
                adsorption_prefactor,
                adsorption_barrier
            )
            rate_calc_time = time.time() - past

            # put rates and sites into one array

            rates = np.concatenate((evaporation_rates, adsorption_rates))
            surface_sites = np.concatenate((solid_sites, solvent_sites))

            # sort rates in ascending order, sort sites according to how rates were sorted

            sorted_indices = rates.argsort()
            sorted_rates, sorted_sites = rates[sorted_indices], surface_sites[sorted_indices]
            cumulative_function = np.cumsum(sorted_rates)
            total_rate = cumulative_function[-1]
            
            if step % dump_every == 0:
            
                min_, q1, median, q3, max_ = np.percentile(sorted_rates, q=[0, 25, 50, 75, 100])
            
                logging.info(f'\tRATE STATISTICS BELOW, N = {len(sorted_rates):.0f}')
                logging.info(f'\t\tmin rate = {min_:.2E}')
                logging.info(f'\t\t25th percentile = {q1:.2E}')
                logging.info(f'\t\t50th percentile = {median:.2E}')
                logging.info(f'\t\t75th percentile = {q3:.2E}')
                logging.info(f'\t\tmax rate = {max_:.2E}')
                logging.info(f'\trate calculation took {rate_calc_time} sec')

            # pick an event according to KMC algorithm
            # bisect_left performs binary search
            # np.random.uniform() is [0, 1), but KMC calls for (0, 1]
            random_draw = 1.0 - np.random.uniform(low=0, high=1)
            event_index = bisect_left(cumulative_function, random_draw * total_rate)

            # update timestep

            second_draw = 1.0 - np.random.uniform(low=0, high=1)
            dt = np.log(1.0 / second_draw) * 1.0 / total_rate
            t += dt

            # swap chosen site

            site_index = sorted_sites[event_index] - 1
            types[site_index] = int(not types[site_index])

        # dump final info

        dump_info(step + 1, t, types, box_bounds, ids, positions, file)

    logging.info('run finished')


def main():

    # uncomment the line below if you are having memory issues
    # mp.set_start_method('forkserver')

    # parse input file
    
    try:
        json_file_name = sys.argv[1]
    except IndexError:
        raise ValueError('needs an input file - see example_molecular.json')
    
    with open(json_file_name, 'r') as input_file:
        real_simulation_parameters = json.load(input_file)
    
    # define dummy simulation parameters
    # same simulation with a smaller geometry and different dump information
    
    dummy_simulation_parameters = copy.deepcopy(real_simulation_parameters)
    dummy_simulation_parameters['box_dimensions'] = [10, 10, 10]
    dummy_simulation_parameters['num_steps'] = 10
    dummy_simulation_parameters['dump_every'] = 1
    dummy_simulation_parameters['dump_file_name'] = f"{real_simulation_parameters['dump_file_name']}.small"
    dummy_simulation_parameters['init_configuration_type'] = 'spherical'
    dummy_simulation_parameters['radius'] = 8.0
    
    # a light simulation is called first, compiling all the expensive function calls
    
    perform_sim(**dummy_simulation_parameters)
    logging.info('dummy simulation completed')

    # perform expensive simulation next
    
    perform_sim(**real_simulation_parameters)
    logging.info('real simulation completed')


if __name__ == '__main__':
    main()
