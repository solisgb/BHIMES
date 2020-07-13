# -*- coding: utf-8 -*-
"""
Created on Sun Jul 12 11:41:37 2020

@author: solis
"""
import numpy as np
import swb

# ======================TEST DATA===========================
ntimestep: int = 1
iamax: float = 1
ia0: float = 0.5
whcmax: float = 5
whc0: float = 1

kdirect: float = 0.
kuz: float = 25
klateral: float = 0.
krunoff: float = 0.

n_array_elements: int = 10
p = np.zeros((n_array_elements,), np.float32)
et = np.zeros((n_array_elements,), np.float32)
rch = np.zeros((n_array_elements,), np.float32)
runoff = np.zeros((n_array_elements,), np.float32)
etr = np.zeros((n_array_elements,), np.float32)

p[1] = 25.
p[2] = 40.
p[3] = 4.

et[:] = 4.

# ======================END TEST DATA========================


def test01():
    storages = np.array([iamax, ia0, whcmax, whc0], np.float32)
    k = np.array([kdirect, kuz, klateral, krunoff], np.float32)

    swb.swb01(ntimestep, storages, k, p, et, rch, runoff, etr)
    print(f'precipitation: {p.sum():0.2f}')
    print(f'et: {et.sum():0.2f}')
    print(f'recharge: {rch.sum():0.2f}')
    print(f'runoff: {runoff.sum():0.2f}')
    print(f'et real: {etr.sum():0.2f}')


if __name__ == "__main__":

    try:
        from datetime import datetime
        from time import time
        import traceback

        now = datetime.now()

        startTime = time()

        test01()

        xtime = time() - startTime
        print(f'El script tardó {xtime:0.1f} s')
    except Exception:
        msg = traceback.format_exc()
        print(f'Exception\n{msg}')
    finally:
        print('\nFin')

