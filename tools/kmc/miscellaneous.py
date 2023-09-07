from dataclasses import dataclass
import time
import numba
import numpy as np
import logging


@dataclass
class Deprecated:

    """
    decorator for a class that will log a warning when an instance of the class is created
     
    reason is a string possibly containing __name__, where __name__ will be
    replaced by the class's name
    """

    reason: str
    
    def __call__(self, class_):
    
        @dataclass
        class WrapperClass(class_):
            
            def __new__(cls):
                
                warning_string = self.reason.replace('__name__', class_.__name__)
                logging.info(warning_string)
                return super(WrapperClass, cls).__new__(cls)
                
        return WrapperClass


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
                self.callback(f'Function {function.__name__}() executed in {present - past} sec')
            return result

        return wrapper_func
        
        
@numba.njit(cache=True)
def any_(x):

    for var in x:
        if var:
            return True
            
    return False
    
    
@numba.njit(cache=True)
def all_(x):

    for var in x:
        if not var:
            return False
            
    return True
        
        
@numba.njit(cache=True)      
def is_close(a, b, rtol=1e-05, atol=1e-08):

    return np.abs(a - b) <= (atol + rtol * np.abs(b))
        
        
@numba.njit(cache=True)
def norm(x):

    return np.sqrt(np.sum(x ** 2))
    
    
@numba.njit(cache=True)
def is_parallel(x, y):

    x1, x2, x3 = x
    y1, y2, y3 = y
    
    cross = np.array([x2 * y3 - x3 * y2, x3 * y1 - x1 * y3, x1 * y2 - x2 * y1])

    return is_close(norm(cross), 0.0)
    
    
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