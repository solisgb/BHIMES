# -*- coding: utf-8 -*-
"""
Created on Sun Jun 21 12:50:31 2020

@author: solis

Before running, set carefully parameter values in bhimes.xml and set project,
     et_proc & annual_graphs in this module

Parameter in this module
procedure: funtion name where soil water balance is made: swb01 or swb24
project: project in bhimes.xml
et_proc: et calculation procedure; implemented ('basic', 'hargreaves')
annual_graphs: if True save png files
"""
import littleLogging as logging

procedure: str = 'swb24'
project: str = 'DHS_QCC'
et_proc: str = 'hargreaves'
annual_graphs: bool = False

if __name__ == "__main__":

    try:
        from datetime import datetime
        from time import time
        import traceback
        from bhimes import BHIMES

        now = datetime.now()

        startTime = time()

        b = BHIMES(project, et_proc)

        # you can comment these lines after the data has been uploaded
        # or you can set upsert atributes to 0 in bhimes.xml
        b.aquifer_upsert_from_file()
        b.outcrop_upsert_from_file()
        b.met_upsert_from_file01()

        # run the swb
        b.swb(procedure)

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
