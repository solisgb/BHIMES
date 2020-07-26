# -*- coding: utf-8 -*-
"""
Created on Sun Jul 12 11:41:37 2020

@author: solis

use simple data sets to call some functions swb module for testing purposes
set your own data or keep the proposed in the test in __main__
changes in testing functions wold have been made after testing and
    therefore some of the functions written in this module can not
    work without modifications
"""
import numpy as np
import swb


def test01(ntimestep, storages, k, p, et, rch, runoff, etr):
    """
    single execution of swb01
    """

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

    storages = np.array([iamax, ia0, whcmax, whc0], np.float32)
    k = np.array([kdirect, kuz, klateral, krunoff], np.float32)

    swb.swb01(ntimestep, storages, k, p, et, rch, runoff, etr)

    print(f'precipitation: {p.sum():0.2f}')
    print(f'et: {et.sum():0.2f}')
    print(f'recharge: {rch.sum():0.2f}')
    print(f'runoff: {runoff.sum():0.2f}')
    print(f'et real: {etr.sum():0.2f}')


def test02():
    """
    sensivity analisys of kuz and shc parameters in func swb01
        this function has its own parameters, set then before call the func
    """
    dir_out = r'H:\tmp'
    kuzmin = 5
    kuzmax = 50
    nkuz = 10
    whcmin = 5
    whcmax = 50
    nwhc = 10

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

    storages = np.array([iamax, ia0, whcmax, whc0], np.float32)
    k = np.array([kdirect, kuz, klateral, krunoff], np.float32)

    ps = swb.Parameter__sensivity(dir_out,
                                  {'kuz': (kuzmin, kuzmax, nkuz),
                                   'whc': (whcmin, whcmax, nwhc)})

    ps.swb01_parameter_sensivity(ntimestep, storages, k, p, et,
                                 rch, runoff, etr)


def test03():
    def ugw_drainage(whcmax, whcr, whc0, kuz, exp, winput, et):
        """
        args
        whxmax: max water holding content mm
        whcr: residual whc
        whc0: initial whc mm
        kuz: satured permeability mm/h
        exp: empirically deduced exponent
        winput: water input mm/h
        et: evapotranspiration mm/h
        output
        whc3: whc at the end
        wd: water drained
        runoff: runoff
        etr: real et
        """
        tiny = 0.00001
        if whcmax < tiny:
            return 0., 0., winput, 0.
        whc1 = whc0 + winput
        whc2 = min(whcmax, whc1)
        runoff = whc1 - whc2
        wd = kuz * ((whc2 - whcr) / (whcmax - whcr))**exp
        wd = min(whc2, wd)
        whc3 = whc2 - wd
        if winput > 0:
            etr = 0.
        else:
            etr = min(whc3, et * whc3 / whcmax)
        whc3 -= etr
        balan = winput - wd - runoff - etr + whc0 - whc3
        if balan > tiny:
            raise ValueError(f'error de balance {balan}:0f')
        return whc3, wd, runoff, etr

    whcmax = 25
    whcr = 1
    whc0 = whcmax * 0.5
    kuz0 = 10
    exp = 0.5
    n = 24
    et = np.empty((n), np.float32)
    et.fill(0.5)
    p = np.zeros((n), np.float32)
    p[4] = 15.

    kuz = kuz0/n
    whc = np.zeros((n+1), np.float32)
    whc[0] = whc0
    wr = np.zeros((n), np.float32)
    runoff = np.zeros((n), np.float32)
    etr = np.zeros((n), np.float32)
    t = np.array([i for i in range(n)], np.float32)
    for i in range(n):
        whc[i+1], wr[i], runoff[i], etr[i] = \
        ugw_drainage(whcmax, whcr, whc[i], kuz, exp, p[i], et[i])
    title = f'whcmax:{whcmax:0.0f}, whc0{whc0:0.0f}, kuz{kuz0:0.0f}, exp:{exp:0.1f}'
    whc_var = whc[1:] - whc[0:-1]

    xyu = [[t, p, 'p'], [t, et,'et'], [t, whc_var,'whc_var'] ]
    xyl = [[t, wr, 'wr'], [t, runoff,'runoff'], [t, etr,'etr']]

    xy2(xyu, xyl, title)
    print(f'p: {p.sum():0.1f}, et: {et.sum():0.1f}')
    print(f'wr: {wr.sum():0.1f}, runoff: {runoff.sum():0.1f}, ',
          f'etr: {etr.sum():0.1f}, ',
          f'whc final: {whc[-1]:0.1f}')


def xy(x, y1, leg1, y2, leg2, y3, leg3, title, xlabel, ylabel):
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots()
    ax.plot(x, y1, label=leg1)
    ax.plot(x, y2, label=leg2)
    ax.plot(x, y3, label=leg3)
    ax.legend()
    ax.set(xlabel=xlabel, ylabel=ylabel, title=title)
    ax.grid()
    plt.show()


def xy2(xy1, xy2, title):
    import matplotlib.pyplot as plt
    fig, _ = plt.subplots()

    plt.suptitle(title)
    plt.subplots_adjust(hspace=0.2, bottom=0.16, top=0.87)

    ax1 = plt.subplot2grid((2, 1), (0, 0))
    ax2 = plt.subplot2grid((2, 1), (1, 0), sharex=ax1) # , sharex=ax1
    ax1.set_ylabel('mm')
    ax2.set_ylabel('mm')

    # subplot superior
    for serie in xy1:
        ax1.plot(serie[0], serie[1], label=serie[2])
        ax1.legend()
    ax1.grid()

    # subplot inferior
    for serie in xy2:
        ax2.plot(serie[0], serie[1], label=serie[2])
        ax2.legend()
    ax2.grid()

    fig.savefig('balance.png')
    plt.show()


if __name__ == "__main__":

    try:
        from datetime import datetime
        from time import time
        import traceback

        # call to functions (change them as you needs)

        now = datetime.now()

        startTime = time()

        test03()

        xtime = time() - startTime
        print(f'El script tardó {xtime:0.1f} s')
    except Exception:
        msg = traceback.format_exc()
        print(f'Exception\n{msg}')
    finally:
        print('\nFin')

