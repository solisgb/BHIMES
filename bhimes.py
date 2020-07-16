# -*- coding: utf-8 -*-
"""
Created on Sun Jun 21 12:51:28 2020

@author: solis

module to evaluate hidrometeorological soil water balance using
    functions in swb module

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

import swb


class BHIMES():
    """
    Balance HIdroMEteorológico en el Suelo
    Soil water balance
    _initial_conditions: initial water conditions in ia and whc
    _et_procedures: procedures to evaluate et
    """
    _initial_conditions = ('dry', 'normal', 'wet')  # don't change the order!
    _et_procedures = ('basic', 'hargreaves')


    def __init__(self, project: str, et_proc: str = 'basic',
                 xml_org: str='bhimes.xml'):
        """
        project: project name in xml_org
        proc: procedure to calculate et
        xml_org: name of xml file with execution parameters
        """
        if et_proc not in BHIMES._et_procedures:
            raise ValueError\
            (f'arg procedure not in {",".join(BHIMES._et_procedures)}')
        self.proc = et_proc
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
        if self.proc == 'basic':
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
        from os.path import dirname, join

        metadata_file = 'swb01_metadata.csv'

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
        fmeta = open(join(dirname(self.db), metadata_file), 'w')
        self._write_metadata_base(fmeta)

        for aquifer in aquifers:
            print(f'{aquifer.name}')

            cur.execute(self.select_outcrops.text, (aquifer.fid,) )
            outcrops = [Outcrop(row) for row in cur.fetchall()]
            cur.execute(self.select_met.text, (aquifer.fid,) )
            met = np.array([row for row in cur.fetchall()])
            dates = met[:, 0]
            imonths = np.array([int(row[5:7])-1 for row in dates], np.int32)
            p = met[:,1].astype(np.float32)
            et = np.empty((p.size), np.float32)
            if self.proc == 'hargreaves':
                # solar radiation
                cur.execute(select_r0, (int(aquifer.y4326),))
                sr = np.array([row[0] for row in cur.fetchall()], np.float32)
                tmin = met[:, 2].astype(np.float32)
                tmax = met[:, 3].astype(np.float32)
                tavg = met[:, 4].astype(np.float32)
                swb.hargreaves_samani(sr, imonths, tmax, tmin, tavg, et)
            else:
                BHIMES._et_averaged_set(imonths, self.et_avg, et)
            rch = np.zeros((p.size), np.float32)
            runoff = np.zeros((p.size), np.float32)
            etr =  np.zeros((p.size), np.float32)
            for outcrop in outcrops:
                c_ia0, c_whc0 = self._coef_initial_wstorages()
                rch1 = np.zeros((p.size), np.float32)
                runoff1 = np.zeros((p.size), np.float32)
                etr1 =  np.zeros((p.size), np.float32)

                storages = np.array([outcrop.ia, outcrop.ia*c_ia0,
                                     outcrop.whc, outcrop.whc*c_whc0],
                                    np.float32)
                k = np.array([outcrop.kdirect, outcrop.kuz,
                              outcrop.klateral, outcrop.krunoff],
                             np.float32)
                ier, xer = swb.swb01(self.time_step, storages, k, p, et, rch1,
                                     runoff1, etr1)
                if ier >= 0:
                    a = (f'Balance error {xer}',  f'Aquifer: {aquifer.name}',
                         f'Outcrop: {outcrop.fid}',
                         f'Date: {dates[ier]} (i {ier})')
                    a = '\n'.join(a)
                    raise ValueError(a)

                rch += rch1 * (outcrop.area * 0.001)
                runoff += runoff1 * (outcrop.area * 0.001)
                etr += etr1 * (outcrop.area * 0.001)
                self._write_metadata(fmeta, aquifer, outcrop)

            self._insert_output(con, cur, aquifer.fid, dates, rch, runoff, etr)

        con.close()
        fmeta.close()


    @jit(void(int32[:], float32[:], float32[:]), nopython=True)
    def _et_averaged_set(imonths, et_avg, et):
        """
        set et arrays whith month averaged values
        """
        for i in range(imonths.size):
            et[i] = et_avg[imonths[i]]


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


    def _write_metadata_base(self, fmeta):
        """
        save data from selected project
        """
        import datetime

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
          f'time_step max: {self.time_step}'
        )
        metadata = '\n'.join(metadata)
        if self.proc == 'basic':
            et_avg = [f'{row:0.2f}' for row in self.et_avg]
            et_avg = ','.join(et_avg)
            metadata = metadata + f'\n{et_avg}'
        fmeta.write(f'{metadata}\n')
        self._write_header(fmeta)


    def _write_header(self, fmeta):
        """
        writes header of metadata
        """
        a = ('aquifer.name', 'aquifer.y4326',
             'outcrop.fid', 'outcrop.area',
             'outcrop.ia', 'outcrop.whc',
             'outcrop.kdirect', 'outcrop.kuz',
             'outcrop.klateral', 'outcrop.krunoff')
        fmeta.write(f'{",".join(a)}\n')


    def _write_metadata(self, fmeta, aquifer, outcrop):
        """
        writes metada
        """
        a = (f'{aquifer.name}', '{aquifer.y4326}',
             f'{outcrop.fid:n}', f'{outcrop.area:0.1f}',
             f'{outcrop.ia:0.2f}', f'{outcrop.whc:0.2f}',
             f'{outcrop.kdirect:0.2f}', f'{outcrop.kuz:0.2f}',
             f'{outcrop.klateral:0.2f}', f'{outcrop.krunoff:0.2f}')
        fmeta.write(f'{",".join(a)}\n')


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

        if self.proc != 'hargreaves':
            print('et annual graphs are only saved when proc is Hargreaves')
            return

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

            # et
            et = np.empty((tmin.size), np.float32)
            swb.hargreaves_samani(r0, im, tmax, tmin, tavg, et)

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
    def krunoff(self):
        return self.data[7]
