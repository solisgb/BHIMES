# -*- coding: utf-8 -*-
import numpy as np
from numba import jit, float32, int32, void
import sqlite3


def set_test_data():
    """
    sets test data to test the function
    """
    select_r0 = \
    """
    select r0
    from r0
    where lat = ?
    order by "month"
    """

    con = sqlite3.connect('r0.db')
    cur = con.cursor()
    cur.execute(select_r0, (38,))
    r0 = np.array([row[0] for row in cur.fetchall()], np.float32)
    im = np.array([i for i in range(12)] + [i for i in range(12)], np.int32)
    tmin = np.ones((24,), np.float32)
    tmax = np.empty((24,), np.float32)
    tmax[:] = 18.
    tavg = np.empty((24,), np.float32)
    tavg = tmax - tmin
    et = np.empty((24,), np.float32)
    hargreaves_samani(r0, im, tmax, tmin, tavg, et)
    print(f'{et.sum():0.2f}')
    np.savetxt('et_test.csv', et, delimiter=',')


@jit(void(float32[:], int32[:], float32[:], float32[:], float32[:],
          float32[:]), nopython=True)
def hargreaves_samani(r0, im, tmax, tmin, tavg, et):
    """
    calculo d ela etp por el metodo de hargreaves-samani
    r0: radiación extraterrestre mm
    im: mes de la observación -1
    tmax, tmin, tavg: tamperatura máxima, míima y media en grados C
    etp: etp calculad en mm
    """
    for i in range(et.size):
        et[i] = 0.0023 * (tavg[i] + 17.78) + r0[im[i]] \
        * (tmax[i] - tmin[i])**0.5


if __name__ == "__main__":

    try:
        from datetime import datetime
        from time import time
        import traceback

        now = datetime.now()

        startTime = time()

        set_test_data()

        xtime = time() - startTime
        print(f'El script tardó {xtime:0.1f} s')

    except Exception:
        msg = traceback.format_exc()
        print(f'Exception\n{msg}')
    finally:
        print('\nFin')
