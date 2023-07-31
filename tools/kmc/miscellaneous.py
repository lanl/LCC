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
        
        
@numba.njit()      
def is_close(a, b, rtol=1e-05, atol=1e-08, equal_nan=False):

    return np.abs(a - b) <= (atol + rtol * np.abs(b))
        
        
@numba.njit()
def norm(x):

    return np.sqrt(np.sum(x ** 2))
    
    
@numba.njit()
def is_parallel(x, y):

    dot_product = np.sum(x * y)
    
    return not is_close(dot_product, 0.0)