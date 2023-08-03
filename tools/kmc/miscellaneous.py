from dataclasses import dataclass
import time
import numba
import numpy as np


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