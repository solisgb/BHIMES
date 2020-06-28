# -*- coding: utf-8 -*-
"""
Created on Sun Jun 21 12:51:28 2020

@author: solis

aquifer data
    fid integer
    xc real, x polygon centroid
    yc real, y polygon centroid
    y4326 real, latitude epgs 4326
    area real, polygon area m2
    name text(100), aquifer name

outcrop data and parameters
    fid integer,
    aquifer integer, aquifer identifier
    lito text(250), lithology
    period text(100), geologic period, age
    perme text(50), permeability mm/d
    area real, polygon area m2
    ia real, initial abstraccion capacity (mainly vegetation) mm
    whc real, soil water holding capacity mm
    kdirect real, direct permeability mm/d
    kuz real, unsatured zone pemeability mm/d
    klateral real, lateral permeability from whc mm/d
    krunoff real, superficial runoff infiltration mm/d

meteorologycal data related with aquifer -in this version-
    fid integer
    date text(10)
    p real, precipitation in date mm
    tmin real, minimum temperature Celsius deg
    tmax real, minimum temperature Celsius deg
    tavg real, average temperature Celsius deg
"""
import sqlite3
from numba import jit
import littleLogging as logging


class BHIMES():
    """
    Balance HIdroMEteorológico en el Suelo
    Soil water balance
    """
    _initial_conditions = ('dry', 'normal', 'wet')  # don't change the order

    def __init__(self, project: str, model: str = 'basic',
                 xml_org: str='bhimes.xml'):
        """
        dst, project database
        tables_2drop, drop table in the list -confirnation required-
            if you drop a table you delete its data without backup them
        """
        self.model = model
        self._read_params(xml_org, project)
        self._create_db()


    def _read_params(self, xml_org: str, project: str):
        import xml.etree.ElementTree as ET

        MAX_TIME_STEP = 24
        tree = ET.parse(xml_org)
        root = tree.getroot()
        prj = None
        for element in root.findall('project'):
            if element.get('name') == project:
                prj = element
                break
        if not prj:
            raise ValueError(f'No se encuentra el project {project}')

        self.description = prj.find('description').text
        self.db = prj.find('db').text
        self.file_aquifers = prj.find('file_aquifers')
        self.file_outcrops = prj.find('file_outcrops')
        self.file_met = prj.find('file_met')
        self.select_aquifers = prj.find('select_aquifers')
        self.select_outcrops = prj.find('select_outcrops')
        self.select_met = prj.find('select_met')
        self.initial_condition = prj.find('initial_condition').text.lower()
        if self.initial_condition not in self._initial_conditions:
            logging.append(f'initial condition changed from ' +\
                           '{self.initial_condition} to normal')
            self.initial_condition = 'normal'
        self.time_step = int(prj.find('time_step').text)
        if self.time_step < 1:
            self.time_step = 1
        elif self.time_step > MAX_TIME_STEP:
            self.time_step = MAX_TIME_STEP
        self.et_avg = \
        [float(item) for item in prj.find('et_avg').text.split(',')]


    def _create_db(self):
        """
        create a sqlite db
        """

        stm1 = \
        """
        create table if not exists aquifer(
            fid integer,
            xc real,
            yc real,
            y4326 real,
            area real,
            name text(100),
            primary key (fid)
        )
        """

        stm2 = \
        """
        create table if not exists outcrop(
            fid integer primary key,
            aquifer integer,
            lito text(250),
            period text(100),
            perme text(50),
            area real,
            ia real,
            whc real,
            kdirect real,
            kuz real,
            klateral real,
            krunoff real,
            foreign key (fid) references aquifer(fid)
        )
        """

        stm3 = \
        """
        create index if not exists outcrop_aquifer
            on outcrop(aquifer)
        """

        stm4 = \
        """
        create table if not exists met(
            fid integer,
            date text(10),
            p real,
            tmin real,
            tmax real,
            tavg real,
            primary key (fid, date)
            foreign key (fid) references aquifer(fid)
        )
        """
        con = sqlite3.connect(self.db)
        con.execute("PRAGMA foreign_keys = 1")
        cur = con.cursor()
        for stm in (stm1, stm2, stm3, stm4):
            try:
                cur.execute(stm)
            except:
                con.close()
                raise ValueError(f'Error al ejecutar {stm}')
        cur.execute('attach database ? as allen', ('r0.db',))
        con.commit()
        con.close()


    def aquifer_upsert_from_file(self):
        """
        inserts or update aquifer data from text file -utf8-
        """

        select1 = \
        """
        select fid from aquifer where fid=?
        """
        update1 = \
        """
        update aquifer set xc=?, yc=?, y4326=?, area=?, name=?
        where fid=?
        """
        insert1 = \
        """
        insert into aquifer(fid, xc, yc, y4326, area, name)
        values (?, ?, ?, ?, ?, ?)
        """
        if int(self.file_aquifers.get('upsert')) != 1:
            return

        org = self.file_aquifers.text
        nskip = int(self.file_aquifers.get('nskip'))
        separator = self.file_aquifers.get('separator')

        con = sqlite3.connect(self.db)
        cur = con.cursor()

        with open(org, 'r') as f:
            for i, line in enumerate(f):
                if i < nskip:
                    continue
                words = line.strip().split(separator)
                fid = int(words[0])
                xc = float(words[1])
                yc = float(words[2])
                y4326 = float(words[3])
                area = float(words[4])
                name = words[5]
                row = cur.execute(select1, (fid,)).fetchone()
                if row:
                    cur.execute(update1, (xc, yc, y4326, area, name, fid))
                else:
                    cur.execute(insert1, (fid, xc, yc, y4326, area, name))

        con.commit()
        con.close()


    def outcrop_upsert_from_file(self, org: str,
                                 nskip: int=1, separator: str=';'):
        """
        inserts or update outcrops data from text file -utf8-
        """

        select1 = \
        """
        select fid from outcrop where fid=?
        """
        update1 = \
        """
        update outcrop set aquifer=?, lito=?, period=?, perme=?, area=?,
        ia=?, whc=?, kdirect=?, kuz=?, klateral=?, krunoff=?
        where fid=?
        """
        insert1 = \
        """
        insert into outcrop(fid, aquifer, lito, period, perme, area, ia, whc,
            kdirect, kuz, klateral, krunoff)
        values (?,?,?,?,?,?,?,?,?,?,?,?)
        """
        if int(self.file_outcrops.get('upsert')) != 1:
            return

        org = self.file_outcrops.text
        nskip = int(self.file_outcrops.get('nskip'))
        separator = self.file_outcrops.get('separator')

        con = sqlite3.connect(self.db)
        cur = con.cursor()
        cur.execute('pragma encoding="UTF-8"')

        with open(org, 'r') as f:
            for i, line in enumerate(f):
                if i < nskip:
                    continue
                words = line.strip().split(separator)
                fid = int(words[0])
                acu = int(words[1])
                lito = words[2]
                period = words[3]
                perme = words[4]
                area = float(words[5])
                ia = float(words[6])
                whc = float(words[7])
                kdirect = float(words[8])
                kuz = float(words[9])
                klateral = float(words[10])
                krunoff = float(words[11])
                row = cur.execute(select1, (fid,)).fetchone()
                if row:
                    cur.execute(update1, (acu, lito, period, perme,
                                          area, ia, whc, kdirect, kuz,
                                          klateral, krunoff, fid))
                else:
                    cur.execute(insert1, (fid, acu, lito, period, perme,
                                          area, ia, whc, kdirect, kuz,
                                          klateral, krunoff))

        con.commit()
        con.close()


    def met_upsert_from_file01(self, org: str,
                                 nskip: int=1, separator: str=';'):
        """
        inserts or update meteorological data from text file -utf8-
        with no average temperature
        """

        select1 = \
        """
        select fid, date from met where fid=? and date=?
        """
        update1 = \
        """
        update met set p=?, tmin=?, tmax=?, tavg=?
        where fid=? and date=?
        """
        insert1 = \
        """
        insert into met(fid, date, p, tmin, tmax, tavg)
        values (?, ?, ?, ?, ?, ?)
        """
        if int(self.file_met.get('upsert')) != 1:
            return

        org = self.file_met.text
        nskip = int(self.file_met.get('nskip'))
        separator = self.file_met.get('separator')

        con = sqlite3.connect(self.db)
        cur = con.cursor()

        with open(org, 'r') as f:
            for i, line in enumerate(f):
                if i < nskip:
                    continue
                words = line.strip().split(separator)
                fid = words[0]
                date = words[1]
                p = float(words[2]) / 10.  # dmm -> mm
                tmin = float(words[3]) / 10.  # dºC -> ºC
                tmax = float(words[4]) / 10.
                tavg = (tmin + tmax) / 2.

                row = cur.execute(select1, (fid, date)).fetchone()
                if row:
                    cur.execute(update1, (p, tmin, tmax, tavg, fid, date))
                else:
                    cur.execute(insert1, (fid, date, p, tmin, tmax, tavg))

        con.commit()
        con.close()


    def swb01(self, table_output: str='swb01'):
        """
        soil water balance version 0.01
        """
        import numpy as np

        select_r0 = \
        """
        select r0
        from allen.r0
        where lat = ?
        order by "month"
        """

        con = sqlite3.connect(self.db)
        cur = con.cursor()
        cur.execute('attach database ? as allen', ('r0.db',))

        cur.execute(self.select_aquifers.text)
        aquifers = [row for row in cur.fetchall()]

        for aquifer in aquifers:
            print(f'{aquifer[1]}')
            cur.execute(select_r0, (int(aquifer[2]),))
            sr = [row for row in cur.fetchall()]
            cur.execute(self.select_outcrops.text, (aquifer[0],) )
            outcrops = [row for row in cur.fetchall()]
            cur.execute(self.select_met.text, (aquifer[0],) )
            met = np.array([row for row in cur.fetchall()])
            dates = met[:, 0]
            months = np.array([int(row[5:7])-1 for row in dates], np.int32)
            p = met[:, 1]
            if self.model != 'basic':
                tmin = met[:, 2]
                tmax = met[:, 3]
                tavg = met[:, 4]
                et = np.empty((p.size), np.float32)
                self.hargreaves_samani_01(sr, months, tmax, tmin, tavg, et)
            rch = np.zeros((p.size), np.float32)
            runoff = np.zeros((p.size), np.float32)
            for outcrop in outcrops:
                rch1 = np.zeros((p.size), np.float32)
                runoff1 = np.zeros((p.size), np.float32)
                if self.model == 'basic':
                    swb_basic(outcrop, outcrop, self.et_avg,
                              months, p, rch1, runoff1)
                else:
                    self.swb01_01(self.time_step, self.initial_condition,
                                  self._initial_conditions,
                                  outcrop, p, et, rch1, runoff1)
                rch += rch1
                runoff += runoff1

        con.close()


    @jit(nopython=True)
    def hargreaves_samani_01(r0, im, tmax, tmin, tavg, etp):
        """
        et by hargreaves-samani
        param
        r0: extraterrestrial radiation mm
        im: for a vector of observations dates, im has the month of each
            date -1; ie, for the date 1970-08-22 -> 8 - 1 = 7
        tmax, tmin, tavg: máx, mín, average temperature ºC
        etp (output): et mm
        """
        for i in range(len(im)):
            etp[i] = 0.0023 * (tavg[i] + 17.78) + r0[im[i]] \
                * (tmax[i] - tmin[i])**0.5


    @jit(nopython=True)
    def swb01_01(time_step, wcondition, initial_conditions, outcrop,
                 p, et, rch, runoff):
        """
        soil water balance in a temporal data serie
        parameters:
        time_step, int: number of iterations to work out water balance in
            the soil
        wcondition str: initial vegatation and soil water condition
        initial conditions list str: all possible initial water conditions
        outcrop np array float: outcrops attributes
            outcrop[iia]: initial abstraction mm
            kdirect: k recharge mm
            kuz: k unsatured zone
            klateral: k lateral applied to kuz
            krunof: k from runoff to kuz
        dates np array str: dates as aaaa-mm-dd
        p, np array float: precipitacion mm
        et, np array float: potential evapotranspiration mm
        rch, np array float: recharge mm -output-
        runoff, np array float: runoff mm -output-
        """
        iia = 1
        iwhc = 2

        ia0, whc0 = BHIMES._initial_wstorages(wcondition, initial_conditions)

        ia = outcrop[iia] * ia0
        whc = outcrop[iwhc] * whc0
        for i in range(1, p.size):
            p_ts, ia = BHIMES._storage_01(outcrop[ia], ia, p[i])  # ia
            et_ts, ia = BHIMES._storage_et_01(ia, et[i])

            nts = BHIMES.nstep_set(time_step, p[i], et[i], whc)
            kdirect = outcrop[3] / nts
            kuz = outcrop[4] / nts
            klateral = outcrop[5] / nts
            krunoff = outcrop[6] / nts
            p_ts = p_ts / nts
            et_ts = et_ts / nts

            for its in (range(time_step)):

                if kdirect > 0:  # direct recharge
                    rch[i] += min(kdirect, p_ts)
                    p_ts -= rch[i]

                p_ts, whc = BHIMES._storage_01(outcrop[iwhc], whc, p_ts) # soil

                if abs(whc - outcrop[iwhc]) < 0.001:  # recharge through unsatured zone
                    x = min(kuz, p_ts)
                    rch[i] += x
                    p_ts = p_ts - x

                if p_ts > 0. and krunoff > 0.:
                    x = min(krunoff, p_ts)
                    p_ts = p_ts - x

                if klateral > 0 and rch[i] > 0:
                    x = min(klateral, rch[i])
                    rch[i] -= x

                runoff[i] += p_ts

                BHIMES._storage_et_01(whc, et_ts)


    @jit(nopython=True)
    def _initial_wstorages(wcondition, initial_conditions):
        if wcondition == initial_conditions[0]:
            ia0, whc0 = 0.1, 0.1
        elif wcondition == initial_conditions[1]:
            ia0, whc0 = 0.5, 0.5
        else:
            ia0, whc0 = 0.9, 0.9
        return ia0, whc0


    @jit(nopython=True)
    def _storage_01(smax, scurrent, wi):
        """
        water balance in a storage
        parameters:
            smax: max water storage
            scurrent: current water storage (at the beginning)
            wi: water input
        output:
            wo: water not in the storage
            sfinal: storage at the end of the water balance
        """
        sdry = smax - wi
        if wi > sdry:
            wo = wi - sdry
            sfinal = smax
        else:
            wo = 0.
            sfinal = scurrent + wi
        return wo, sfinal


    @jit(nopython=True)
    def _storage_et_01(scurrent, pwr):
        """
        soil water release -et-
        parameters:
            scurrent: current water storage (at the beginning)
            pwr: potential water release
        output:
            pwr_final: pwr at the end of the balance
            sfinal: storage at the end of the water balance
        """
        if pwr > scurrent:
            pwr_final = pwr - scurrent
            sfinal = 0.
        else:
            pwr_final = 0.
            sfinal = scurrent - pwr
        return pwr_final, sfinal


    @jit(nopython=True)
    def nstep_set(nstep0, p, et, whc):
        """
        sets nstep according to p, et, whc
        """
        if nstep0 == 1:
            return 1

        minp = 2.
        minet = 0.1 * whc
        if minet < 0.25:
            minet = 0.25

        if p < minp:
            nstep = 1
        else:
            for i in range(nstep0):
                nstep = i + 1
                x = et / nstep
                if x < minet:
                    break
        return nstep


    @jit(nopython=True)
    def swb_basic(kuz, whc, et_avg, months, p, rch, runoff):
        """
        basic soil water balance in a temporal data serie
        """
        whc0 = 0.1 * whc
        for i in range(1, p.size):
            p1 = p[i]
            x = whc - whc0
            if p1 > x:
                p1 -= x
                whc0 = whc
            else:
                p1 = 0
                whc0 += p1
            rch[i] = min(kuz, p1)
            p1 -= rch[i]
            if p1 > 0:
                runoff[i] = p1
            whc0 = max(0., whc0-et_avg[months[i]])
