"""Minimal NumPy stand-in for the parts of CuPy that run_NNLS.py uses.

Only for tests: it lets the GPU wrapper code (_GpuBatchedSolver) run on a CPU-only
machine. Arrays are plain NumPy arrays, so semantics match CuPy for this subset.
"""
import types

import numpy as np

__version__ = 'fake'

uint8 = np.uint8
asarray = np.asarray
asfortranarray = np.asfortranarray
linalg = types.SimpleNamespace(norm=np.linalg.norm)


def asnumpy(a):
    return np.asarray(a)


class _Stream:
    def __init__(self, ptds=False):
        self.ptds = ptds

    def use(self):
        return self


class _Device:
    id = 0


class _Runtime:
    TOTAL_MEMORY = 16 * (1 << 30)

    @staticmethod
    def memGetInfo():
        return _Runtime.TOTAL_MEMORY, _Runtime.TOTAL_MEMORY

    @staticmethod
    def getDeviceProperties(device_id):
        return {'name': b'Fake GPU (NumPy)', 'major': 0, 'minor': 0,
                'totalGlobalMem': _Runtime.TOTAL_MEMORY}

    @staticmethod
    def driverGetVersion():
        return 0

    @staticmethod
    def runtimeGetVersion():
        return 0


cuda = types.SimpleNamespace(Stream=_Stream, Device=_Device, runtime=_Runtime,
                             set_allocator=lambda allocator=None: None)
