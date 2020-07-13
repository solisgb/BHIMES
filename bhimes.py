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
    name text, aquifer name

outcrop data and parameters
    fid integer,
    aquifer integer, aquifer identifier
    lito text, lithology
    period text, geologic period, age
    perme text, permeability mm/d
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
from numba import jit, void, float32, int32
import numpy as np


class BHIMES():
    """
    Balance HIdroMEteorológico en el Suelo
    Soil water balance
    _initial_conditions: initial water conditions in ia and whc
    _procedures: calculus procedures
    """
    _initial_conditions = ('dry', 'normal', 'wet')  # don't change the order!
    _procedures = ('basic', 'swb01')


    def __init__(self, project: str, procedure: str = 'basic',
                 xml_org: str='bhimes.xml'):
        """
        project: project name in xml_org
        procedure: to calculate soil water balance
        xml_org: name of xml file with execution parameters
        """
        if procedure not in BHIMES._procedures:
            raise ValueError\
            (f'arg procedure not in {",".join(BHIMES._procedures)}')
        self.proc = procedure
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
        self.db = prj.find('db').text.strip()
        self.file_aquifers = prj.find('file_aquifers')
        self.file_outcrops = prj.find('file_outcrops')
        self.file_met = prj.find('file_met')
        self.select_aquifers = prj.find('select_aquifers')
        self.select_outcrops = prj.find('select_outcrops')
        self.select_met = prj.find('select_met')
        self.initial_condition = prj.find('initial_condition').text\
        .strip().lower()
        if self.initial_condition not in BHIMES._initial_conditions:
            raise ValueError(f'initial condition not in ' +\
                             '{",".join(BHIMES.self._initial_condition)}')
        self.time_step = int(prj.find('time_step').text)
        if self.time_step < 1:
            self.time_step = 1
        elif self.time_step > MAX_TIME_STEP:
            self.time_step = MAX_TIME_STEP
        self.et_avg = \
        np.array([item for item in prj.find('et_avg').text.split(',')],
                 np.float32)
        self.table_output = prj.find('table_output').text.strip()
        self.xy_annual_dir = prj.find('xy_annual_dir').text.strip()


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
            name text,
            primary key (fid)
        )
        """

        stm2 = \
        """
        create table if not exists outcrop(
            fid integer primary key,
            aquifer integer,
            lito text,
            period text,
            perme text,
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
        separator = self.file_aquifers.get('sep')

        con = sqlite3.connect(self.db)
        cur = con.cursor()
        cur.execute("PRAGMA auto_vacuum = FULL")

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


    def outcrop_upsert_from_file(self):
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
        separator = self.file_outcrops.get('sep')

        con = sqlite3.connect(self.db)
        cur = con.cursor()
        cur.execute("PRAGMA auto_vacuum = FULL")
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


    def met_upsert_from_file01(self):
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
        separator = self.file_met.get('sep')

        con = sqlite3.connect(self.db)
        cur = con.cursor()
        cur.execute("PRAGMA auto_vacuum = FULL")

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


    def swb01(self):
        """
        soil water balance version 0.01
        """

        select_r0 = \
        """
        select r0
        from allen.r0
        where lat = ?
        order by "month"
        """

        con = sqlite3.connect(self.db)
        cur = con.cursor()
        cur.execute("PRAGMA auto_vacuum = FULL")
        cur.execute('attach database ? as allen', ('r0.db',))
        self._create_output_table(con, cur)

        cur.execute(self.select_aquifers.text)
        aquifers = [Aquifer(row) for row in cur.fetchall()]

        for aquifer in aquifers:
            print(f'{aquifer.name}')

            cur.execute(self.select_outcrops.text, (aquifer.fid,) )
            outcrops = [Outcrop(row) for row in cur.fetchall()]
            cur.execute(self.select_met.text, (aquifer.fid,) )
            met = np.array([row for row in cur.fetchall()])
            dates = met[:, 0]
            imonths = np.array([int(row[5:7])-1 for row in dates], np.int32)
            p = met[:,1].astype(np.float32)
            if self.proc != 'basic':
                # solar radiation
                cur.execute(select_r0, (int(aquifer.y4326),))
                sr = np.array([row[0] for row in cur.fetchall()], np.float32)
                tmin = met[:, 2].astype(np.float32)
                tmax = met[:, 3].astype(np.float32)
                tavg = met[:, 4].astype(np.float32)
                et = np.empty((p.size), np.float32)
                BHIMES._hargreaves_samani(sr, imonths, tmax, tmin, tavg, et)
            rch = np.zeros((p.size), np.float32)
            runoff = np.zeros((p.size), np.float32)
            etr =  np.zeros((p.size), np.float32)
            for outcrop in outcrops:
                c_ia0, c_whc0 = self._coef_initial_wstorages()
                rch1 = np.zeros((p.size), np.float32)
                runoff1 = np.zeros((p.size), np.float32)
                etr1 =  np.zeros((p.size), np.float32)
                if self.proc == 'basic':
                    BHIMES._swb_basic(outcrop.kuz, outcrop.whc*c_whc0,
                                      outcrop.whc,
                                      self.et_avg, imonths, p, rch1, runoff1,
                                      etr1)
                else:
                    storages = np.array([outcrop.ia, outcrop.ia*c_ia0,
                                         outcrop.whc, outcrop.whc*c_whc0],
                                        np.float32)
                    k = np.array([outcrop.kdirect, outcrop.kuz,
                                  outcrop.klateral, outcrop.krunoff],
                                 np.float32)
                    self._swb01_01(self.time_step, storages, k, p, et,
                                   rch1, runoff1, etr1)
                rch += rch1 * (outcrop.area * 0.001)
                runoff += runoff1 * (outcrop.area * 0.001)
                etr += etr1 * (outcrop.area * 0.001)

            self._insert_output(con, cur, aquifer.fid, dates, rch, runoff, etr)

        self._save_metadata(con, cur)

        con.close()


    @jit(void(float32[:], int32[:], float32[:], float32[:], float32[:],
              float32[:]), nopython=True)
    def _hargreaves_samani(r0, im, tmax, tmin, tavg, et):
        """
        et by hargreaves-samani
        param
        r0: extraterrestrial radiation mm
        im: for a vector of observations dates, im has the month of each
            date -1; ie, for the date 1970-08-22 -> 8 - 1 = 7
        tmax, tmin, tavg: máx, mín, average temperature ºC
        etp (output): et mm
        """
        for i in range(et.size):
            et[i] = 0.0023 * (tavg[i] + 17.78) + r0[im[i]] \
            * (tmax[i] - tmin[i])**0.5


#    @jit(void(int32, float32[:], float32[:], float32[:],
#              float32[:], float32[:], float32[:], float32[:]),
#        nopython=True)
    def _swb01_01(ntimestep, storages, k, p, et, rch, runoff, etr):
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
            ia1, p1 = BHIMES._storage_input(iamax, ia1, p[i])  # ia
            ia1, etr[i] = BHIMES._storage_output(ia1, et[i])
            et1 = et[i] - etr[i]

            nts = BHIMES._nstep_set(ntimestep, p, whcmax, whc1)
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

                whc1, pi = BHIMES._storage_input(whcmax, whc1, pi)

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

                whc1, eti = BHIMES._storage_output(whc1, eti)
                x = et1 - eti
                etr[i] -= x


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
    def _storage_output(whc0, pwr):
        """
        soil water release -et-
        parameters:
            whc0: initial water holding content
            pwr: potential water release -max water release-
        output:
            shc1: final water holding content
            wr: water release
        """
        if pwr >= whc0:
            wr = whc0
            whc1 = 0.
        else:
            wr = pwr
            whc1 = whc0 - pwr
        return whc1, wr


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


    @jit(void(float32, float32, float32, float32[:], int32[:],
              float32[:], float32[:], float32[:], float32[:] ), nopython=True)
    def _swb_basic(kuz, whc0, whc, et_avg, im, p, rch, runoff, etr):
        """
        basic soil water balance in a temporal data serie
        parameters:
            kuz: k unsatured zone mm
            whc0: initial whc
            whc: soil water holding content mm
            et_avg: average values of et mm -12 values, Jan. to Dec.-
            im: for a vector of observations dates,
                im contains the month of each date -1;
                i.e., for the date 1970-08-22 -> 8 - 1 = 7
        p: precipitacion mm
        rch, np array float: recharge mm -output-
        runoff: runoff mm -output-
        et: real evapotranspiration mm -output-
        """
        whc1 = whc0
        for i in range(p.size):
            p1 = p[i]
            x = whc - whc1
            if p1 > x:
                p1 -= x
                whc1 = whc
            else:
                p1 = 0
                whc1 += p1
            rch[i] = min(kuz, p1)
            p1 -= rch[i]
            if p1 > 0:
                runoff[i] = p1
            if et_avg[im[i]] >= whc1:
                etr[i] = whc1
                whc1 = 0.
            else:
                etr[i] = et_avg[im[i]]
                whc1 -= et_avg[im[i]]


    def _create_output_table(self, con, cur):
        """
        create the output table
        parameters
        con: connection -already open-
        cur: cursor to conn -already open-
        """

        stm1 = f'drop table if exists {self.table_output}'

        stm2 =\
        f"""
        create table if not exists {self.table_output}(
            aquifer integer,
            date text(10),
            rch real,
            runoff real,
            etr, real,
            primary key (aquifer, date)
        )
        """
        cur.execute(stm1)
        cur.execute(stm2)
        con.commit()


    def _insert_output(self, con, cur, aquifer, dates, rch, runoff, etr):
        """
        insert or update output in output table
        parameters
        con: connection -already open-
        cur: cursos to conn -already open-
        dates, np array str: dates (aaaa-mm-dd)
        rch, np array float: calculated recharge mm
        runoff, np array float: calculated runoff mm
        etr, np array float: calculated et mm -real evapotranspiration-

        There is an elegant solution at
        https://stackoverflow.com/questions/18621513/python-insert-numpy-array
        -into-sqlite3-database
        witch I'm not going to apply because I want to access database outside
        python
        """

        stm1 =\
        f"""
        insert into {self.table_output} (aquifer, date, rch, runoff, etr)
        values (?, ?, ?, ?, ?)
        """

        for i in range(dates.size):
            cur.execute(stm1, (aquifer, dates[i],
                        float(rch[i]), float(runoff[i]), float(etr[i])))
        con.commit()


    def _save_metadata(self, con, cur):
        """
        con: connection to an open database
        cur: cursor to con
        """
        import datetime

        stm1 = f'drop table if exists {self.table_output}_metadata'

        stm2 =\
        f"""
        create table if not exists {self.table_output}_metadata(
            description text
        )
        """

        stm3 =\
        f"""
        insert into {self.table_output}_metadata (description)
        values (?)
        """
        et_avg = [f'{row:0.2f}' for row in self.et_avg]
        et_avg = ','.join(et_avg)
        metadata =\
        (datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
         self.description,
         'calculation procedure: ' + self.proc,
         'aquifers imported from: ' + self.file_aquifers.text,
         'outcrops imported from: ' + self.file_outcrops.text,
         'meteorogical data imported from: ' + self.file_met.text,
         'aquifers select: ' + self.select_aquifers.text,
         'outcrops select: ' + self.select_outcrops.text,
         'met select: ' + self.select_met.text,
         'initial_condition: ' + self.initial_condition,
          f'time_step max: {self.time_step}',
         'et_avg: ' + et_avg
        )
        metadata = '\n'.join(metadata)
        cur.execute(stm1)
        cur.execute(stm2)
        cur.execute(stm3, (metadata,))
        con.commit()


    def save_annual_graphs(self):
        from os.path import join

        stm1 = \
        f"""
        select strftime('%Y', r.date) "year",
            sum(r.rch) rch, sum(r.runoff ) runoff , sum(r.etr) etr
        from aquifer a
        	left join {self.table_output} r on (a.fid = r.aquifer )
        where a.fid = ?
        group by strftime('%Y', r.date)
        order by strftime('%Y', r.date)
        """

        con = sqlite3.connect(self.db)
        cur = con.cursor()
        cur.execute(self.select_aquifers.text)
        aquifers = [Aquifer(row) for row in cur.fetchall()]

        for aquifer in aquifers:
            print(f'{aquifer.name} xy')
            cur.execute(stm1, (aquifer.fid,))
            rows = np.array([row for row in cur.fetchall()])
            years = rows[:, 0].astype(np.int32)
            rch = rows[:, 1].astype(np.float32)
            runoff = rows[:, 2].astype(np.float32)
            etr = rows[:, 3].astype(np.float32)

            title = f'Recarga en el acuífero {aquifer.name}'
            ylabel = 'Rec m3/a'
            dst = join(self.xy_annual_dir, f'{aquifer.name}_annual_rch.png')
            BHIMES._xy_plot_1g(title, years, rch, ylabel, dst)

            title = f'Escorrentía en el acuífero {aquifer.name}'
            ylabel = 'Esc m3/a'
            dst = join(self.xy_annual_dir, f'{aquifer.name}_annual_runoff.png')
            BHIMES._xy_ts_plot_1g(title, years, runoff, ylabel, dst)

            title = f'ETR en el acuífero {aquifer.name}'
            ylabel = 'ETR m3/a'
            dst = join(self.xy_annual_dir, f'{aquifer.name}_annual_etr.png')
            BHIMES._xy_ts_plot_1g(title, years, etr, ylabel, dst)

        con.close()


    def save_annual_data_graphs(self):
        from os.path import join

        stm1 = \
        f"""
        select strftime('%Y', m.date) "year",
            sum(m.p) p, avg(m.tmin ) tmin, avg(m.tmax ) tmax, avg(m.tavg ) tavg
        from aquifer a
        	left join met m on (a.fid = m.fid )
        where a.fid = ?
        group by a.fid, a.name, strftime('%Y', m.date)
        order by a.fid, a.name, strftime('%Y', m.date)
        """

        con = sqlite3.connect(self.db)
        cur = con.cursor()
        cur.execute(self.select_aquifers.text)
        aquifers = [Aquifer(row) for row in cur.fetchall()]

        for aquifer in aquifers:
            print(f'{aquifer.name} xy data')
            cur.execute(stm1, (aquifer.fid,))
            rows = np.array([row for row in cur.fetchall()])
            years = rows[:, 0].astype(np.int32)
            p = rows[:, 1].astype(np.float32)
            tavg = rows[:, 2].astype(np.float32)
            tmin = rows[:, 3].astype(np.float32)
            tmax = rows[:, 4].astype(np.float32)

            title = f'Precipitación en el acuífero {aquifer.name}'
            ylabel = 'P mm/a'
            dst = join(self.xy_annual_dir, f'{aquifer.name}_annual_p.png')
            BHIMES._xy_plot_1g(title, years, p, ylabel, dst)

            title = f'Temperatura media en el acuífero {aquifer.name}'
            ylabel = 'T ºC/a'
            dst = join(self.xy_annual_dir, f'{aquifer.name}_annual_tavg.png')
            BHIMES._xy_ts_plot_1g(title, years, tavg, ylabel, dst)

            title = f'Temperatura mín. media en el acuífero {aquifer.name}'
            ylabel = 'T ºC/a'
            dst = join(self.xy_annual_dir, f'{aquifer.name}_annual_tmin.png')
            BHIMES._xy_ts_plot_1g(title, years, tmin, ylabel, dst)

            title = f'Temperatura máx. media en el acuífero {aquifer.name}'
            ylabel = 'T ºC/a'
            dst = join(self.xy_annual_dir, f'{aquifer.name}_annual_tmax.png')
            BHIMES._xy_ts_plot_1g(title, years, tmax, ylabel, dst)

        con.close()


    def save_annual_eth_graphs(self):
        from os.path import join

        create_table = \
        """
        create table if not exists et(
            date text(10),
            et real,
            primary key(date)
        )
        """
        delete_table = 'delete from et'
        insert_data = 'insert into et (date, et) values (?, ?)'
        drop_table = 'drop table if exists et'

        select_r0 = \
        """
        select r0
        from allen.r0
        where lat = ?
        order by "month"
        """

        stm1 = \
        f"""
        select  m.date, m.tmin, m.tmax, m.tavg
        from met m
        where m.fid = ?
        order by m.date
        """

        stm2 = \
        """
        select strftime('%Y', date) "year", avg(et) etavg
        from et
        group by strftime('%Y', date)
        order by strftime('%Y', date)
        """

        cont = sqlite3.connect('memory')
        curt = cont.cursor()
        curt.execute(create_table)

        con = sqlite3.connect(self.db)
        cur = con.cursor()
        cur.execute('attach database ? as allen', ('r0.db',))
        cur.execute(self.select_aquifers.text)
        aquifers = [Aquifer(row) for row in cur.fetchall()]

        for aquifer in aquifers:
            print(f'{aquifer.name} xy eth')

            # solar radiation -Allen-
            cur.execute(select_r0, (int(aquifer.y4326),))
            r0 = np.array([row[0] for row in cur.fetchall()], np.float32)

            # t diarias
            cur.execute(stm1, (aquifer.fid,))
            rows = np.array([row for row in cur.fetchall()])
            dates = rows[:, 0]
            im = np.array([int(row[5:7])-1 for row in dates], np.int32)
            tmin = rows[:, 1].astype(np.float32)
            tmax = rows[:, 2].astype(np.float32)
            tavg = rows[:, 3].astype(np.float32)

            # eth
            et = np.empty((tmin.size), np.float32)
            BHIMES._hargreaves_samani(r0, im, tmax, tmin, tavg, et)

            # xy
            curt.execute(delete_table)
            for i in range(dates.size):
                curt.execute(insert_data, (dates[i], et[i]))
            cont.commit()
            curt.execute(stm2)
            rows = np.array([row for row in curt.fetchall()])
            years = rows[:, 0].astype(np.int32)
            et = rows[:, 1].astype(np.float32)

            title = f'ET Hargreaves-Samani en el acuífero {aquifer.name}'
            ylabel = 'ET med mm/a'
            dst = join(self.xy_annual_dir, f'{aquifer.name}_annual_eth.png')
            BHIMES._xy_plot_1g(title, years, et, ylabel, dst)

        curt.execute(drop_table)
        con.close()
        cont.close()


    @staticmethod
    def _xy_ts_plot_1g(title: str, x: list, y: list, ylabel: str, dst: str):
        """
        Dibuja una figura con 1 gráfico (axis) xy
        args
        title: título de la figura
        x: lista de objetos date
        y: lista de valores interpolados float
        dst: nombre fichero destino (debe incluir la extensión png)
        """
        import warnings
        import matplotlib.pyplot as plt
        import matplotlib as mpl
        # parámetros específicos
        warnings.simplefilter(action='ignore', category=FutureWarning)
        mpl.rc('font', size=8)
        mpl.rc('axes', labelsize=8, titlesize= 10, grid=True)
        mpl.rc('axes.spines', right=False, top=False)
        mpl.rc('xtick', direction='out', top=False)
        mpl.rc('ytick', direction='out', right=False)
        mpl.rc('lines', linewidth=0.8, linestyle='-', marker='.', markersize=4)
        mpl.rc('legend', fontsize=8, framealpha=0.5, loc='best')

        fig, ax = plt.subplots()

        plt.suptitle(title)
        ax.set_ylabel(ylabel)

        fig.autofmt_xdate()

        ax.plot(x, y)

        fig.savefig(dst)
        plt.close('all')
        plt.rcdefaults()


    @staticmethod
    def _xy_plot_1g(title: str, x: list, y: list, ylabel: str, dst: str):
        """
        Draws an xy grpah
        args
        title: figure title
        x: numpy array
        y: numpy array
        dst: png file name
        """
        import matplotlib.pyplot as plt
        import matplotlib as mpl

        mpl.rc('font', size=8)
        mpl.rc('axes', labelsize=8, titlesize= 10, grid=True)
        mpl.rc('axes.spines', right=False, top=False)
        mpl.rc('xtick', direction='out', top=False)
        mpl.rc('ytick', direction='out', right=False)
        mpl.rc('lines', linewidth=0.8, linestyle='-', marker='.', markersize=4)
        mpl.rc('legend', fontsize=8, framealpha=0.5, loc='best')

        fig, ax = plt.subplots()

        plt.suptitle(title)
        ax.set_ylabel(ylabel)

        ax.plot(x, y)

        fig.savefig(dst)
        plt.close('all')
        plt.rcdefaults()



class Aquifer():
    """
    acuíferos
    """


    def __init__(self, data: list):
        """
        fid integer: identificador
        name text: aquifer name
        y4326 real: latitude epgs 4326
        """
        self.data = tuple([row for row in data])


    @property
    def fid(self):
        return self.data[0]


    @property
    def name(self):
        return self.data[1]


    @property
    def y4326(self):
        return self.data[2]


class Outcrop():
    """
    afloramientos permeables
    """


    def __init__(self, data: list):
        """
        fid: identifier
        area: outcrop area m2
        ia: initial abstraction mm
        kdirect: k recharge mm
        kuz: k unsatured zone mm
        klateral: k lateral applied to kuz mm
        krunof: k from runoff to kuz mm
        """
        self.data = tuple([row for row in data])


    @property
    def fid(self):
        return self.data[0]


    @property
    def area(self):
        return self.data[1]


    @property
    def ia(self):
        return self.data[2]


    @property
    def whc(self):
        return self.data[3]


    @property
    def kdirect(self):
        return self.data[4]


    @property
    def kuz(self):
        return self.data[5]


    @property
    def klateral(self):
        return self.data[6]


    @property
    def runoff(self):
        return self.data[7]
