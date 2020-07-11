# -*- coding: utf-8 -*-
"""
Created on Sun Jun 21 12:50:31 2020

@author: solis

to execute, complete parameters in bhimes.xml and set project and proc values
project: projet in bhimes.xml -several are possible-
proc: calculation procedure; implemented ('basic', 'swb01')
annual_graphs: if True save png files
"""
import littleLogging as logging

project: str = 'BDA202006'
proc: str = 'basic'
annual_graphs: bool = False

if __name__ == "__main__":

    try:
        from datetime import datetime
        from time import time
        import traceback
        from bhimes import BHIMES

        now = datetime.now()

        startTime = time()

        b = BHIMES(project, proc)
        b.aquifer_upsert_from_file()  # controlled in xml
        b.outcrop_upsert_from_file()
        b.met_upsert_from_file01()
        b.swb01()

        if annual_graphs:
            b.save_annual_graphs()  # xy recharge, runoff & ret
            b.save_annual_data_graphs()  # xy p, tmax, tmin, tavg
            b.save_annual_eth_graphs()  # xy pet (hargreaves)

        xtime = time() - startTime
        print(f'El script tardó {xtime:0.1f} s')

    except ValueError:
        msg = traceback.format_exc()
        logging.append(f'ValueError exception\n{msg}')
    except ImportError:
        msg = traceback.format_exc()
        logging.append(f'ImportError exception\n{msg}')
    except Exception:
        msg = traceback.format_exc()
        logging.append(f'Exception\n{msg}')
    finally:
        logging.dump()
        print('\nFin')
