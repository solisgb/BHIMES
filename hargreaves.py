# -*- coding: utf-8 -*-
from numba import jit


@jit(nopython=True)
def hargreaves_samani_01(r0, im, tmax, tmin, tavg, etp):
    """
    calculo d ela etp por el metodo de hargreaves-samani
    r0: radiación extraterrestre mm
    im: mes de la observación -1
    tmax, tmin, tavg: tamperatura máxima, míima y media en grados C
    etp: etp calculad en mm
    """
    for i in range(len(im)):
        etp[i] = 0.0023 * (tavg[i] + 17.78) + r0[im[i]] \
            * (tmax[i] - tmin[i])**0.5
