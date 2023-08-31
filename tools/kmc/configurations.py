# module dealing with initial configurations

from dataclasses import dataclass
import numpy as np


__available_config_types__ = [
    'spherical',
    'dump_file'
]


@dataclass
class Configuration:

    def set_types(positions, box_bounds):
    
        pass


@dataclass
class Spherical(Configuration):

    radius: float
    
    def set_types(self, positions, box_bounds):
    
        origin = np.array(box_bounds) / 2
        
        types = np.zeros(len(positions))
        
        for index, position in enumerate(positions):
            if np.linalg.norm(position - origin) <= self.radius:
                types[index] = 1
                
        return types
        
        
@dataclass
class DumpFile(Configuration):

    config_loading_file: str
    
    def set_types(self, positions, box_bounds):
    
        types = np.zeros(len(positions))
        
        with open(self.config_loading_file, 'r') as file:
            lines = file.readlines()
            
        starting_index = -1
            
        for index, line in enumerate(lines):
            if 'ITEM: TIMESTEP' in line:
                starting_index = index
                
        block = lines[starting_index:]
        data = np.loadtxt(block, skiprows=9)
        
        positions_in_file = data[:, 2:]
        
        ids = np.array(data[:, 0], dtype=int)
        
        if positions[ids - 1] != positions_in_file:
            raise ValueError('dump file and current lattice are incompatible')
        
        types[ids - 1] = 1
        
        return types
        
        
def get_configuration(init_configuration_type: str, **kwargs) -> Configuration:

    if init_configuration_type.lower() == 'spherical':
    
        configuration = Spherical(radius=kwargs['radius'])
        
    elif init_configuration_type.lower() == 'dump_file':
        
        configuration = DumpFile(config_loading_file=kwargs['config_loading_file'])
        
    else:
    
        raise ValueError('valid initial configuration types are {__available_config_types}')
        
    return configuration
            