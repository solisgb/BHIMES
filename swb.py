# -*- coding: utf-8 -*-
"""
Created on Sun Jul 12 12:37:35 2020

@author: solis
"""
from os.path import join
import sqlite3

from numba import jit, void, float32, int32
import numpy as np


def swb24(whcmax, whcr, whc0, kuz, kdirect, exp, p, nh, et, wd, runoff, etr):
    """
    water soil balance of daily data
    the water balance is calculated in an intradaily basis controlled by
        the constant ndaily_div; for instance, if ndaily_div is 24,
        the balance is calculated each hour. Be carefull, if you modify
        ndaily_div you must modify the nh values accordingly
    args
    in
        whxmax: soil max water holding content -whc- mm
        whcr: residual whc mm
        whc0: initial whc mm
        kuz: soil satured permeability mm/d
        exp: empirically deduced exponent
        p: precipitation -water input- mm/d
        nh: number of daily divisions to calculate the balance
        et: evapotranspiration mm/d
    out
        wd: water drained
        runoff: runoff
        etr: real et
    returns
        balance error in day i hour j
        i: element i (day) where a balance error is raised
        j: element j (hour) where a balance error is raised
        if no error balance is raised it returns 0., -1, -1
    """
    no_error = 0.
    ndaily_div = 24
    no_ij = -1
    if whcmax < 0.1:
        wd.fill(0.)
        runoff[:] = p[:]
        etr.fill(0.)
        return no_error, no_ij, no_ij

    tiny = 0.0001
    kuzh = kuz / ndaily_div
    p_nd = np.empty((ndaily_div), np.float32)
    wd_nd = np.empty((ndaily_div), np.float32)
    runoff_nd = np.empty((ndaily_div), np.float32)
    etr_nd = np.empty((ndaily_div), np.float32)
    whce = whcmax - whcr
    for i in range(p.size):
        p_nd.fill(0.)
        if p[i] > 0.:
            p1h = p[i] / nh[i]
            p_nd[0:nh[i]] = p1h
        et1h = et[i] / ndaily_div
        for j in range(p_nd.size):
            whc1 = whc0 + p_nd[j]
            whc2 = min(whcmax, whc1)
            runoff_nd[j] = whc1 - whc2
            x1 = whc2 - whcr
            wd1 = kuzh * (x1 / whce)**exp
            if p_nd[j] > 0:
                etr1 = 0.
            else:
                etr1 = min(x1, et1h * x1 / whce)
            out1 = wd1 + etr1
            out2 = min(x1, out1)
            if out1 > out2:
                wd1 = wd1 * wd1 / out1
                etr1 = etr1 * etr1 / out1
            wd_nd[j] = wd1
            etr_nd[j] = etr1
            whc3 = whc2 - wd1 - etr1
            balan = p_nd[j] - wd_nd[j] - runoff_nd[j] - etr_nd[j] + whc0 - whc3
            if abs(balan) > tiny:
                return balan, i, j
            whc0 = whc3

        wd[i] = wd_nd.sum()
        runoff[i] = runoff_nd.sum()
        etr[i] = etr_nd.sum()

    return no_error, no_ij, no_ij


def swb24001(whcmax, whcr, whc0, kuz, kdirect, exp, p, nh, et, wd, runoff,
             etr):
    """
    deprecated
    hourly water soil balance of daily data
    args
    in
        whxmax: soil max water holding content -whc- mm
        whcr: residual whc mm
        whc0: initial whc mm
        kuz: soil satured permeability mm/d
        exp: empirically deduced exponent
        p: precipitation -water input- mm/d
        nh: number of daily divisions with rain -for instance,
            if ndaily_div == 24 we are considering hours-
        et: evapotranspiration mm/d
    out
        wd: water drained
        runoff: runoff
        etr: real et
    returns
        balance error in day i hour j
        i: element i (day) where a balance error is raised
        j: element j (hour) where a balance error is raised
        if no error balance is raised it returns 0., -1, -1
    """
    no_error = 0.
    ndaily_div = 24
    no_ij = -1
    if whcmax < 0.1:
        wd.fill(0.)
        runoff[:] = p[:]
        etr.fill(0.)
        return no_error, no_ij, no_ij

    tiny = 0.0001
    kuzh = kuz / ndaily_div
    ph = np.empty((ndaily_div), np.float32)
    wdh = np.empty((ndaily_div), np.float32)
    runoffh = np.empty((ndaily_div), np.float32)
    etrh = np.empty((ndaily_div), np.float32)
    whce = whcmax - whcr
    for i in range(p.size):
        ph.fill(0.)
        if p[i] > 0.:
            p1h = p[i] / nh[i]
            ph[0:nh[i]] = p1h
        et1h = et[i] / ndaily_div
        for j in range(ph.size):
            whc1 = whc0 + ph[j]
            whc2 = min(whcmax, whc1)
            runoffh[j] = whc1 - whc2
            x1 = whc2 - whcr
            x2 = kuzh * (x1 / whce)**exp
            wdh[j] = min(x1, x2)
            whc3 = whc2 - wdh[j]
            if ph[j] > 0:
                etrh[j] = 0.
            else:
                x1 = whc3 - whcr
                etrh[j] = min(x1, et1h * x1 / whce)
            whc3 -= etrh[j]
            balan = ph[j] - wdh[j] - runoffh[j] - etrh[j] + whc0 - whc3
            if abs(balan) > tiny:
                return balan, i, j
            whc0 = whc3

        wd[i] = wdh.sum()
        runoff[i] = runoffh.sum()
        etr[i] = etrh.sum()

    return no_error, no_ij, no_ij


@jit(nopython=True)
def swb01(ntimestep, storages, k, p, et, rch, runoff, etr):
    """
    A new soil water balance function in a temporal data serie
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
    returns
        2 values, depending on water balance in each i loop
        integer: if a balance problem occurs it returns the i value in the loop
            in witch water balance fails; else returns -1 at the end of the
            function
        float: water balance
        the function can be runned with jit activated; at the moment numba
            doesn't suport to raise excepcions, so I need to return the
            pair integer, float to signal an error in the balance

    ia: initial abstraction
    whc: water holding content
    iamax, whcmax: max values of ia and whc (data)
    ia1, whc1: initial values of ia and whc, then change in function
    """
    no_error = 0.
    no_i = -1
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
            whc_initial = whc1
            pi = p1
            eti = et1

            if kdirect > 0. and pi > 0.:  # direct recharge
                rd = min(kdirect, pi)
                pi -= rd
            else:
                rd = 0.

            whc1, pi = _storage_input(whcmax, whc1, pi)

            # if soil is full of water and there is precipitation excess
            if abs(whc1-whcmax) < 0.0001:
                ruz = min(kuz, pi)
                pi -= ruz
                iruz = 1  # flag
            else:
                ruz = 0.
                iruz = 0

            if krunoff > 0. and pi > 0.:
                rr = min(krunoff, pi)
                pi -= rr
            else:
                rr = 0.

            rnf = pi  # used in the water balance
            runoff[i] += pi
            pi = 0.  # pedagogical assignment

            whc1, etr_soil = _storage_output(whcmax, whc1, eti)
            etr[i] += etr_soil

            if iruz == 0:
                whc1, rdr = _ugw_drainage(whcmax, whc1, kuz)
            else:
                rdr = 0.

            rt2aq = rd + ruz + rr + rdr

            if klateral > 0.:
                rl = min(klateral, rt2aq)
                rt2aq -= rl
            else:
                rl = 0.

            rch[i] += rt2aq

            balan = p1 - rt2aq - rl - rnf - etr_soil + whc_initial - whc1
            if abs(balan) > 0.001:
                return balan, i

    return no_error, no_i


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
        whc1: final water holding content
        wr: water release
    """
    tiny = 0.00001
    if whcmax < tiny:
        return 0., 0.
    xwr = pwr * whc0 / whcmax
    if xwr >= whc0:
        wr = whc0
        whc1 = 0.
    else:
        wr = xwr
        whc1 = whc0 - xwr
    return whc1, wr


@jit(nopython=True)
def _ugw_drainage(whcmax, whc0, kuz, exp: float = 12.):
    """
    under ground water drainage
    parameters:
        whcmax: max water holding content
        whc0: initial water holding content
        kuz: permeability unsatured zone
    output:
        whc0: final water holding content
        wd: water drained
    """
    tiny = 0.00001
    if whcmax < tiny:
        return 0., 0.
    n = 25
    whc = np.zeros((n), np.float32)
    whc[0] = whc0
    wd = np.zeros((n), np.float32)
    for i in range(1, n):
        k = kuz * (whc[i-1] / whcmax)**exp
        wd[i] = min(whc[i-1], k)
        whc[i] = whc[i-1] - wd[i]
    return whc[-1], wd.sum()


@jit(nopython=True)
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


@jit(void(float32[:], int32[:], float32[:], float32[:], float32[:],
          float32[:]), nopython=True)
def hargreaves_samani(r0, im, tmax, tmin, tavg, et):
    """
    et by hargreaves-samani
    param
    r0: extraterrestrial radiation mm
    im: for a vector of daily observation dates, im has the month of each
        date -1; ie, for the date 1970-08-22 -> 8 - 1 = 7
    tmax, tmin, tavg: máx, mín, average temperature ºC
    etp (output): et mm
    """
    for i in range(et.size):
        et[i] = 0.0023 * (tavg[i] + 17.78) + r0[im[i]] \
        * (tmax[i] - tmin[i])**0.5


class Parameter__sensivity():
    """
    runs soil water balance function in a given range of parameter values
    the results are stored in an csv or sqlite file with the function name;
    parameters implemented: whc and kuz
    """


    def __init__(self, dir_out: str, param: dict, description: str='test'):
        """
        dir_out: directory with data
        param: param to test
            for each element: key parameter name, values a sequence of 3
            element(min value, max value, num of tries between min and max)
        """
        self.dir_out = dir_out
        self.param = param
        self.description = description


    def pmin(self, par_name: str):
        return min(self.param[par_name][0:2])


    def pmax(self, par_name: str):
        return max(self.param[par_name][0:2])

    def n(self, par_name: str):
        return self.param[par_name][2]


    def delta_get(self, par_name: str):
        return (self.pmax(par_name) - self.pmin(par_name)) / self.n(par_name)


    def swb01_parameter_sensivity(self,ntimestep, storages, k, p, et,
                                  rch, runoff, etr, output_type='csv'):
        """
        calls to swb01 function using a range of predefined values of
            whc and kuz arameters
        to modify the parameters or add more, it's need necessary
            to modify the function
        """

        if output_type not in ('csv', 'sqlite'):
            raise ValueError(f'{output_type} is not a valid output_type')

        dbname = join(self.dir_out, 'swb01_parameter_sensivity'+'.db')
        tname = 'swb01'

        drop_table = f'drop table if exists {tname}'

        create_table = \
        f"""
        create table if not exists {tname} (
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
            ndata integer,
            psum real,
            npgt0 real,
            etsum real,
            rchsum real,
            runoffsum real,
            etrsum real,
            primary key(fid)
        )
        """

        insert = f"""
        insert into {tname}
        (fid, ntimestep, iamax, ia0, whcmax, whc0, kdirect, kuz, klateral,
         krunoff, ndata, psum, npgt0, etsum, rchsum, runoffsum, etrsum)
        values (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """

        ndata = p.size
        psum = p.sum()
        npgt0 = np.count_nonzero(p > 0.)
        etsum = np.sum(et)

        n = 0
        x1 = self.delta_get('kuz')
        x2 = self.delta_get('whc')
        toxy = []

        if output_type == 'csv':
            fo = open(join(self.dir_out, 'swb01_sensivity_analisys.csv'), 'w')
            fo.write('fid,ntimestep,iamax,ia0,whcmax,whc0,kdirect,kuz,'
                     'klateral,krunoff,ndata,psum,npgt0,etsum,rchsum,'
                     'runoffsum,etrsum\n')
        else:
            con = sqlite3.connect(dbname)
            cur = con.cursor()
            cur.execute(drop_table)
            cur.execute(create_table)

        for i in range(self.n('kuz')):
            xk = np.copy(k)
            xk[1] = self.pmin('kuz') + (x1 * i)
            for j in range(self.n('whc')):
                n += 1
                print(f'{n:n}')
                xstorages = np.copy(storages)
                xstorages[2] = self.pmin('whc') + (j * x2)
                rch[:] = 0.
                runoff[:] = 0.
                etr[:] = 0.
                ier, xer = swb01(ntimestep, xstorages, xk, p, et, rch, runoff,
                                 etr)
                if ier >= 0:
                    raise ValueError(f'Balance error in iter {ier:n}',
                                     f'Balance error {xer:0.2f}'
                                     f'{i:n} call loop')

                if output_type != 'csv':
                    cur.execute(insert, (n, ntimestep,
                                     float(xstorages[0]), float(xstorages[1]),
                                     float(xstorages[2]), float(xstorages[3]),
                                     float(xk[0]), float(xk[1]), float(xk[2]),
                                     float(xk[3]),
                                     ndata, float(psum), int(npgt0),
                                     float(etsum),
                                     float(rch.sum()), float(runoff.sum()),
                                     float(etr.sum())))
                else:
                    fo.write(f'{n:n},{ntimestep:n},{xstorages[0]:0.2f},'
                             f'{xstorages[1]:0.2f},{xstorages[2]:0.2f},'
                             f'{xstorages[3]:0.2f},{xk[0]:0.2f},{xk[1]:0.2f},'
                             f'{xk[2]:0.2f},{xk[3]:0.2f},'
                             f'{ndata:n},{psum:0.2f},{npgt0:n},{etsum:0.2f},'
                             f'{rch.sum():0.2f},{runoff.sum():0.2f},'
                             f'{etr.sum():0.2f}\n')

                toxy.append([xstorages[2], xk[1], rch.sum(), runoff.sum(),
                             etr.sum()])
                np.savetxt(join(self.dir_out, f'{n:n}_{tname}.csv'),
                           np.transpose([p, et, rch, runoff, etr]),
                           delimiter=',', fmt='%0.1f',
                           header='p,et,recarge,runoff,etr')

        if output_type != 'csv':
            con.commit()
            con.close()
        else:
            fo.close()

        toxy = np.array(toxy, np.float32)
        for i, item in enumerate(('recharge', 'runoff', 'etr')):
            _contour(f'{self.description} {item}', toxy[:, 0], toxy[:, 1],
                     toxy[:, i+2], 'whc', 'kuz',
                     join(self.dir_out, f'swb01_contour_{item}'))


def _contour(title, x, y, z, xlabel, ylabel, dst, scale: float=1.0):
    """
    3D representarion of the results
    """
    import matplotlib.pyplot as plt

    fig, ax1 = plt.subplots()

    ax1.tricontour(x*scale, y*scale, z*scale, levels=14, linewidths=0.5,
                   colors='k')
    cntr1 = ax1.tricontourf(x, y, z, levels=14, cmap="RdBu_r")

    fig.colorbar(cntr1, ax=ax1)
    ax1.plot(x, y, 'ko', ms=3)
    ax1.set_title(title)
    ax1.set_xlabel(xlabel)
    ax1.set_ylabel(ylabel)

    plt.subplots_adjust(hspace=0.5)
    plt.show()
    fig.savefig(dst)
    plt.close('all')


@jit(nopython=True)
def nhours_generator01(p, nh_gt1):
    """
    A generator of number of hours of rain
    """
    nh = np.zeros((p.size), np.int32)
    j = -1
    for i in range(p.size):
        if p[i] > 0:
            if p[i] <= 1.:
                nh[i] = 1
            else:
                j += 1
                nh[i] = nh_gt1[j]
    return nh


def swb01_001(ntimestep, storages, k, p, et, rch, runoff, etr):
    """
    Deprecated
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
    returns
        2 values, depending on water balance in each i loop
        integer: if a balance problem occurs it returns the i value in the loop
            in witch water balance fails; else returns -1 at the end of the
            function
        float: water balance
        the function can be runned with jit activated; at the moment numba
            doesn't suport to raise excepcions, so I need to return the
            pair integer, float to signal an error in the balance

    ia: initial abstraction
    whc: water holding content
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
            whc_initial = whc1
            pi = p1
            eti = et1

            if kdirect > 0. and pi > 0.:  # direct recharge
                rd = min(kdirect, pi)
                rch[i] += rd
                pi -= rd
            else:
                rd = 0.

            whc1, pi = _storage_input(whcmax, whc1, pi)

            # if soil is full of water and there is precipitation excess
            if pi > 0. and abs(whc1-whcmax) < 0.0001:
                ruz = min(kuz, pi)
                pi -= ruz

                if krunoff > 0. and pi > 0.:
                    rr = min(krunoff, pi)
                    pi -= rr
                else:
                    rr = 0.

                rt = ruz + rr

                if klateral > 0 and rt > 0.:
                    rl = min(klateral, rt)
                    rt -= rl
                else:
                    rl = 0.

                rch[i] += rt

                if pi > 0.:
                    rnf = pi
                    runoff[i] += pi
                    pi = 0.  # pedagogical assignment
                else:
                    rnf = 0.
            else:
                rt = rl = rnf = 0.

            whc1, etr_soil = _storage_output(whcmax, whc1, eti)
            etr[i] += etr_soil
            balan = p1 - rd - rt - rl - rnf - etr_soil + whc_initial - whc1
            if abs(balan) > 0.001:
                return i, balan
    return -1, balan