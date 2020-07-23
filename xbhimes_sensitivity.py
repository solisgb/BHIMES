# -*- coding: utf-8 -*-
"""
Created on Sun Jun 21 12:50:31 2020

@author: solis

before running, set carefully parameter values in bhimes.xml and set project &
     et_proc in this module

to run this option you need the sensivity element in bhimes.xml
    to be defined

Parameter in this module
project: project in bhimes.xml
et_proc: et calculation procedure; implemented ('basic', 'hargreaves')
"""
import littleLogging as logging

project: str = 'DHS_QCC'
et_proc: str = 'hargreaves'
output: str = r'H:\off\phdhs2127\recarga\cc_q\output_serie_completa'


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
        # or you can set upsert atribute to 0 in bhimes.xml
        b.aquifer_upsert_from_file()
        b.outcrop_upsert_from_file()
        b.met_upsert_from_file01()

        # sensitivity
        b.swb01_sensitivity(output)

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
