# module dealing with simulation lattices

from dataclasses import dataclass
import numpy as np
from itertools import product
import pickle


__available_lattice_types__ = [
    'triclinic_four_indices',
    'orthogonal_four_indices',
    'triclinic_three_indices',
    'orthogonal_three_indices',
    'petn_lattice',
    'pickled_lattice'
]


@dataclass
class Lattice:

    def initialize_simulation(self) -> tuple:
    
        """
        initialize_simulation() returns arrays of types, positions, box bounds, and site identifiers
        """

        pass
        
        
@dataclass
class TriclinicFourIndices(Lattice):

    dimensions: np.ndarray
    a: float
    b: float
    c: float
    alpha: float
    beta: float
    gamma: float
    offset: np.ndarray
    
    def initialize_simulation(self):
    
        # create particle labels
        # convoluted way, but all that matters is that we have 2 * box_dims[0] * box_dims[1] * box_dims[2] unique labels
        # along with corresponding types (initialized to 0) and positions

        particle_labels = product(*(*(np.arange(d) for d in self.dimensions), np.arange(2)))
        positions = np.zeros((np.product(self.dimensions) * 2, 3))
        
        # convert angles to radians
        
        alpha, beta, gamma = self.alpha * np.pi / 180.0, self.beta * np.pi / 180.0, self.gamma * np.pi / 180.0
        
        # formulas from http://www.aflowlib.org/prototype-encyclopedia/triclinic_lattice.html

        a = self.a * np.array([1, 0, 0])
        b = self.b * np.array([np.cos(gamma), np.sin(gamma), 0])
        cx = self.c * np.cos(beta)
        cy = self.c * (np.cos(alpha) - np.cos(beta) * np.cos(beta)) / np.sin(gamma)
        cz = np.sqrt(self.c ** 2 - cx ** 2 - cy ** 2)
        c = np.array([cx, cy, cz])
        
        # populate positions with i * a + j * b + k * c + offset * l

        for particle_index, (i, j, k, l) in enumerate(particle_labels):
            positions[particle_index] = i * a + j * b + c * k + self.offset * l

        # give each site a unique identifier, starting at 1
        num_sites = len(positions)
        ids = np.arange(1, num_sites + 1)

        # simulation box bounds are extrema of the initialized position array
        box_bounds = np.max(positions, axis=0)
        
        # need to add some skin for periodic boundary conditions
        
        box_bounds += a + b + c

        return positions, box_bounds, ids
        
        
@dataclass
class OrthogonalFourIndices(Lattice):

    dimensions: np.ndarray
    a: float
    b: float
    c: float
    offset: np.ndarray
    
    def initialize_simulation(self):
    
        # orthogonal lattice is just a special case of the triclinic lattice
    
        triclinic_lattice = TriclinicFourIndices(
            dimensions=self.dimensions,
            a=self.a,
            b=self.b,
            c=self.c,
            alpha=90.0,
            beta=90.0,
            gamma=90.0,
            offset=self.offset
        )
        
        return triclinic_lattice.initialize_simulation()
        
        
@dataclass
class TriclinicThreeIndices(Lattice):

    dimensions: np.ndarray
    a: float
    b: float
    c: float
    alpha: float
    beta: float
    gamma: float
    
    def initialize_simulation(self):
    
        # create particle labels
        # convoluted way, but all that matters is that we have 2 * box_dims[0] * box_dims[1] * box_dims[2] unique labels
        # along with corresponding types (initialized to 0) and positions

        particle_labels = product(*(np.arange(d) for d in self.dimensions))
        positions = np.zeros((np.product(self.dimensions) * 2, 3))
        
        # convert angles to radians
        
        alpha, beta, gamma = self.alpha * np.pi / 180.0, self.beta * np.pi / 180.0, self.gamma * np.pi / 180.0

        # populate positions with i * a + j * b + k * c + l * dr
        
        # http://www.aflowlib.org/prototype-encyclopedia/triclinic_lattice.html

        a = self.a * np.array([1, 0, 0])
        b = self.b * np.array([np.cos(gamma), np.sin(gamma), 0])
        cx = self.c * np.cos(beta)
        cy = self.c * (np.cos(alpha) - np.cos(beta) * np.cos(beta)) / np.sin(gamma)
        cz = np.sqrt(self.c ** 2 - cx ** 2 - cy ** 2)
        c = np.array([cx, cy, cz])

        for particle_index, (i, j, k) in enumerate(particle_labels):
            positions[particle_index] = i * a + j * b + c * k

        # give each site a unique identifier, starting at 1
        num_sites = len(positions)
        ids = np.arange(1, num_sites + 1)

        # simulation box bounds are extrema of the initialized position array
        box_bounds = np.max(positions, axis=0)
        
        # need to add some skin for periodic boundary conditions
        
        box_bounds += a + b + c

        return positions, box_bounds, ids
        
        
@dataclass
class OrthogonalThreeIndices(Lattice):

    dimensions: np.ndarray
    a: float
    b: float
    c: float
    
    def initialize_simulation(self):
        
        # orthogonal lattice is just a special case of the triclinic lattice
    
        triclinic_lattice = TriclinicThreeIndices(
            dimensions=self.dimensions,
            a=self.a,
            b=self.b,
            c=self.c,
            alpha=90.0,
            beta=90.0,
            gamma=90.0
        )
        
        return triclinic_lattice.initialize_simulation()
        
        
@dataclass
class PETNLattice(Lattice):

    dimensions: np.ndarray

    def initialize_simulation(self):
    
        # PETN parameters, just for ease of input file
        
        # one can write lattice_type = "petn_lattice" instead of specifying the
        # params in the input file
    
        petn_lattice = OrthogonalFourIndices(
            dimensions=self.dimensions,
            a=9.088,
            b=9.088,
            c=6.737,
            offset=[4.54348, 4.54346, 3.36908]
        )
        
        return petn_lattice.initialize_simulation()
        
        
@dataclass
class PickledLattice(Lattice):

    lattice_loading_file: str
    
    def initialize_simulation(self):
    
        with open(self.lattice_loading_file, 'rb') as l_file:
            positions, box_bounds, ids = pickle.load(l_file)
            
        return positions, box_bounds, ids
            
        
def get_lattice(lattice_type: str, box_dimensions: list, **kwargs) -> Lattice:

    if lattice_type.lower() == 'triclinic_four_indices':
    
        lattice_kwargs = {
            'a': kwargs['a'],
            'b': kwargs['b'],
            'c': kwargs['c'],
            'alpha': kwargs['alpha'],
            'beta': kwargs['beta'],
            'gamma': kwargs['gamma'],
            'offset': np.array(kwargs['offset']),
            'dimensions': np.array(box_dimensions)
        }
        
        return TriclinicFourIndices(**lattice_kwargs)
        
    elif lattice_type.lower() == 'orthogonal_four_indices':
    
        lattice_kwargs = {
            'a': kwargs['a'],
            'b': kwargs['b'],
            'c': kwargs['c'],
            'offset': np.array(kwargs['offset']),
            'dimensions': np.array(box_dimensions)
        }
        
        return OrthogonalFourIndices(**lattice_kwargs)
        
    elif lattice_type.lower() == 'triclinic_three_indices':
    
        lattice_kwargs = {
            'a': kwargs['a'],
            'b': kwargs['b'],
            'c': kwargs['c'],
            'alpha': kwargs['alpha'],
            'beta': kwargs['beta'],
            'gamma': kwargs['gamma'],
            'dimensions': np.array(box_dimensions)
        }
        
        return TriclinicThreeIndices(**lattice_kwargs)
        
    elif lattice_type.lower() == 'orthogonal_three_indices':
    
        lattice_kwargs = {
            'a': kwargs['a'],
            'b': kwargs['b'],
            'c': kwargs['c'],
            'dimensions': np.array(box_dimensions)
        }
        
        return OrthogonalThreeIndices(**lattice_kwargs)
        
    elif lattice_type.lower() == 'petn_lattice':
    
        return PETNLattice(dimensions=box_dimensions)
        
    elif lattice_type.lower() == 'pickled_lattice':
        
        return PickledLattice(lattice_loading_file=kwargs['lattice_loading_file'])

    else:

        raise ValueError(f'valid lattice_type inputs are {__available_lattice_types__}')
