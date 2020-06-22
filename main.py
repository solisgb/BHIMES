# -*- coding: utf-8 -*-
"""
Created on Sun Jun 21 12:50:31 2020

@author: solis
"""
import littleLogging as logging

db = r'H:\off\balan\ch_202006.db'
org1 = r'H:\off\balan\acuiferos.txt'
org2 = r'H:\off\balan\afloramientos.txt'
org3 = r'H:\off\balan\met_centroids.txt'

if __name__ == "__main__":

    try:
        from datetime import datetime
        from time import time
        import traceback
        from bhimes import BHIMES

        now = datetime.now()

        startTime = time()

        b = BHIMES()
        b.create_db(db)
        b.aquifer_upsert_from_file(db, org1, separator='\t')
        b.outcrop_upsert_from_file(db, org2, separator='\t')
        b.met_upsert_from_file01(db, org3, separator='\t')

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
