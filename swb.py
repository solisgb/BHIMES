# -*- coding: utf-8 -*-
"""
Created on Sun Jul 12 12:37:35 2020

@author: solis
"""
from numba import jit, void, float32, int32
import numpy as np

#    @jit(void(int32, float32[:], float32[:], float32[:],
#              float32[:], float32[:], float32[:], float32[:]),
#        nopython=True)
def swb01(ntimestep, storages, k, p, et, rch, runoff, etr):
    """
    soil water balance in a temporal data serie
    parameters:
        ntimestep: number of iterations to work out water balance in
            the soil
        storages: array of water storages values -look at initialization-
        k: array of k values -look at initialization-
        p: precipitacion mm
        et: potential evapotranspiration mm
        rch: recharge mm -output-
        runoff: runoff mm -output-
        etr: real evapotranspiration mm -output-

    iamax, whcmax: max values of ia and whc (data)
    ia1, whc1: initial values of ia and whc, then change in function
    """
    iamax = storages[0]
    ia1 = storages[1]
    whcmax = storages[2]
    whc1 = storages[3]
    for i in range(p.size):
        ia1, p1 = _storage_input(iamax, ia1, p[i])  # ia
        ia1, etr[i] = _storage_output(iamax, ia1, et[i])
        et1 = et[i] - etr[i]

        nts = _nstep_set(ntimestep, p1, whcmax, whc1)
        if nts > 1:
            kdirect = k[0] / nts
            kuz = k[1] / nts
            klateral = k[2] / nts
            krunoff = k[3] / nts
            p1 = p1 / nts
            et1 = et1 / nts
        else:
            kdirect = k[0]
            kuz = k[1]
            klateral = k[2]
            krunoff = k[3]

        for i_time_step in range(nts):
            pi = p1
            eti = et1

            if kdirect > 0. and pi > 0.:  # direct recharge
                rch[i] += min(kdirect, pi)
                pi -= rch[i]

            whc1, pi = _storage_input(whcmax, whc1, pi)

            # if soil is full of water and there is precipitation excess
            if pi > 0. and abs(whc1-whcmax) < 0.0001:
                r1 = min(kuz, pi)
                pi -= r1

                if krunoff > 0. and pi > 0.:
                    r2 = min(krunoff, pi)
                    pi -= r2
                else:
                    r2 = 0.

                r3 = r1 + r2

                if klateral > 0 and r3 > 0.:
                    r4 = min(klateral, r3)
                    r3 -= r4

                rch[i] += r3

                if pi > 0.:
                    runoff[i] += pi
                    pi = 0.  # pedagogical assignment

            whc1, eti = _storage_output(whcmax, whc1, eti)
            x = et1 - eti
            etr[i] += x


def _coef_initial_wstorages(self):
    """
    sets coefs. in function of water initial condition
        then
    """
    if self.initial_condition == self._initial_conditions[0]:
        ia0, whc0 = 0.1, 0.1
    elif self.initial_condition == self._initial_conditions[1]:
        ia0, whc0 = 0.5, 0.5
    else:
        ia0, whc0 = 0.9, 0.9
    return ia0, whc0


@jit(nopython=True)
def _storage_input(smax, scurrent, wi):
    """
    water balance in a storage
    parameters:
        smax: max water storage
        scurrent: current water storage (at the beginning)
        wi: water input
    output:
        we: water excess
        sfinal: storage at the end of the water balance
    """
    sdry = smax - scurrent
    if wi >= sdry:
        we = wi - sdry
        sfinal = smax
    else:
        we = 0.
        sfinal = scurrent + wi
    return sfinal, we


@jit(nopython=True)
def _storage_output(whcmax, whc0, pwr):
    """
    soil water release -et-
    parameters:
        whcmax: max water holding content
        whc0: initial water holding content
        pwr: potential water release -max water release-
    output:
        shc1: final water holding content
        wr: water release
    """
    xwr = pwr * whc0 / whcmax
    if xwr >= whc0:
        wr = whc0
        whc1 = 0.
    else:
        wr = xwr
        whc1 = whc0 - xwr
    return whc1, wr


@jit(int32(int32, float32, float32, float32), nopython=True)
def _nstep_set(nstep0, p, whcmax, whc0):
    """
    sets nstep according to p, et, whc
    """
    if nstep0 == 1:
        return 1

    minp = 1.
    if p < minp:
        return 1

    sdry = whcmax - whc0
    if p < sdry:
        return 1
    else:
        for i in range(nstep0, 1, -1):
            x = sdry / i
            if x >= minp:
                return i
    return i


class Parameter__sensivity():


    def __init__(self, kuzmax, kuzmin, nkuz,
                 whcmax, whcmin, nwhc):
        self.kuzmax = kuzmax
        self.kuzmin = kuzmin
        self.nkuz = nkuz
        self.whcmax = whcmax
        self.whcmin = whcmin
        self.nwhc = nwhc


    def swb01_parameter_sensivity(self,ntimestep, storages, k, p, et,
                                  rch, runoff, etr):
        from inspect import function
        import sqlite3

        fname = function.func_name

        drop_table = f'drop table if exists {fname}'

        create_table = \
        f"""
        create table if not exists {fname} (
            fid integer,
            ntimestep integer,
            iamax real,
            ia0 real,
            whcmax real,
            whc0 real,
            kdirect real,
            kuz real,
            klateral real,
            krunoff real,
            p real,
            np real,
            et real,
            rch real,
            runoff real,
            etr real,
            primary key(fid)
        )
        """

        insert = f"""
        insert into {fname}
        (fid, ntimestep, iamax, ia0, whcmax, whc0, kdirect, kuz, klateral,
         krunoff, p, np, et, rch, runoff, etr)
        values (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """

        psum = p.sum()
        np = np.count_nonzero(p > 0.)
        etavg = np.average(et)

        con = sqlite3.connect(function.func_name + '.db')
        cur = con.cursor()
        cur.execute(drop_table)
        cur.execute(create_table)

        x1 = (self.kuzmax - self.kuzmin) / self.nkuz
        for i in range(self.nkuz):
            kuzmax = self.kuzmin + (x1 * i)
            xk = np.copy(k)
            xk[1] = kuzmax
            x2 = (self.whcmax - self.whcmin) / self.nwhc
            for j in range(self.nwhc):
                whc = self.whcmin + (j * x2)
                xstorages = np.copy(storages)
                xstorages[2] = whc
                swb01(ntimestep, xstorages, xk, p, et, rch, runoff, etr)

                cur.execute(insert, (ntimestep,
                                     xstorages[0], xstorages[1],
                                     xstorages[2], xstorages[3],
                                     xk[0], xk[1], xk[2], xk[3],
                                     p.sum(), ))