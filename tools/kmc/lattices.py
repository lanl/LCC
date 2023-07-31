from dataclasses import dataclass
import numpy as np
from itertools import product


@dataclass
class Lattice:

    def initialize_simulation(self):
    
        """
        initialize_simulation() returns arrays of types, positions, box bounds, and site identifiers
        """

        pass


@dataclass
class PETNMolecularLattice(Lattice):

    dimensions: np.ndarray
    a: float = 9.088
    b: float = 9.088
    c: float = 6.737
    offset: np.ndarray = np.array([4.54348, 4.54346, 3.36908])
    a_dir: np.ndarray = np.array([1, 0, 0])
    b_dir: np.ndarray = np.array([0, 1, 0])
    c_dir: np.ndarray = np.array([0, 0, 1])

    def initialize_simulation(self):

        # create particle labels
        # convoluted way, but all that matters is that we have 2 * box_dims[0] * box_dims[1] * box_dims[2] unique labels
        # along with corresponding types (initialized to 0) and positions

        particle_labels = product(*(*(np.arange(d) for d in self.dimensions), np.arange(2)))
        types = np.zeros(np.product(self.dimensions) * 2, dtype=np.int_)
        positions = np.zeros((np.product(self.dimensions) * 2, 3))

        # populate positions with i * a + j * b + k * c + l * dr

        a = self.a * self.a_dir
        b = self.b * self.b_dir
        c = self.c * self.c_dir

        for particle_index, (i, j, k, l) in enumerate(particle_labels):
            positions[particle_index] = i * a + j * b + c * k + self.offset * l

        # give each site a unique identifier, starting at 1
        num_sites = len(positions)
        ids = np.arange(1, num_sites + 1)

        # simulation box bounds are extrema of the initialized position array
        box_bounds = np.max(positions, axis=0)

        return types, positions, box_bounds, ids
        
        
@dataclass
class PETNBlockLattice(Lattice):

    dimensions: np.ndarray
    a: float = 9.088
    b: float = 9.088
    c: float = 6.737
    a_dir: np.ndarray = np.array([1, 0, 0])
    b_dir: np.ndarray = np.array([0, 1, 0])
    c_dir: np.ndarray = np.array([0, 0, 1])
    
    def initialize_simulation(self):
    
        # create particle labels
        # convoluted way, but all that matters is that we have box_dims[0] * box_dims[1] * box_dims[2] unique labels
        # along with corresponding types (initialized to 0) and positions

        particle_labels = product(*(np.arange(d) for d in self.dimensions))
        types = np.zeros(np.product(self.dimensions), dtype=np.int_)
        positions = np.zeros((np.product(self.dimensions), 3))

        # populate positions with i * a + j * b + k * c

        a = self.a * self.a_dir
        b = self.b * self.b_dir
        c = self.c * self.c_dir

        for particle_index, (i, j, k) in enumerate(particle_labels):
            positions[particle_index] = i * a + j * b + c * k

        # give each site a unique identifier, starting at 1
        num_sites = len(positions)
        ids = np.arange(1, num_sites + 1)

        # simulation box bounds are extrema of the initialized position array
        box_bounds = np.max(positions, axis=0)

        return types, positions, box_bounds, ids