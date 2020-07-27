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

meteorological data related with aquifer -in this version-
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
    _select_r0 = 'select r0 from allen.r0 where lat = ? order by "month"'
    _tname_iters_ts = {'swb01': 'swb01_iters', 'swb24': 'swb24_iters'}
    _tname_ts = {'swb01': 'swb01_ts', 'swb24': 'swb24_ts'}


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
        from os.path import join
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

        self.dir_out = prj.find('dir_db').text
        self.description = prj.find('description').text
        db = prj.find('db').text.strip()
        self.db = join(self.dir_out, db)
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
        sensitivity = prj.find('sensitivity')
        if sensitivity:
            self.par = {'whc': (float(sensitivity.find('whc').get('delta')),
                                int(sensitivity.find('whc').get('neval'))
                               ),
                        'kuz': (float(sensitivity.find('kuz').get('delta')),
                                int(sensitivity.find('kuz').get('neval'))
                               )
                       }
            element = sensitivity.find('bcexp')
            if element:
                self.par['bcexp'] = (float(element.get('delta')),
                                     int(element.get('neval'))
                                    )

    def delta(self, par_name: str):
        return self.par[par_name][0]


    def neval(self, par_name: str):
        return self.par[par_name][1]


    def _create_db(self):
        """
        creates a sqlite db
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
            whcr real,
            exp real,
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
        the column meaning is explained in class Outcrop
        """

        select1 = \
        """
        select fid from outcrop where fid=?
        """
        update1 = \
        """
        update outcrop set aquifer=?, lito=?, period=?, perme=?, area=?,
        ia=?, whc=?, kdirect=?, kuz=?, klateral=?, krunoff=?, whcr=?, exp=?
        where fid=?
        """
        insert1 = \
        """
        insert into outcrop(fid, aquifer, lito, period, perme, area, ia, whc,
            kdirect, kuz, klateral, krunoff, whcr, exp)
        values (?,?,?,?,?,?,?,?,?,?,?,?,?,?)
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
                whcr = float(words[12])
                exp = float(words[13])
                row = cur.execute(select1, (fid,)).fetchone()
                if row:
                    cur.execute(update1, (acu, lito, period, perme,
                                          area, ia, whc, kdirect, kuz,
                                          klateral, krunoff, whcr, exp, fid))
                else:
                    cur.execute(insert1, (fid, acu, lito, period, perme,
                                          area, ia, whc, kdirect, kuz,
                                          klateral, krunoff, whcr, exp))

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
                p = float(words[2])  # dmm -> mm
                tmin = float(words[3])  # dºC -> ºC
                tmax = float(words[4])
                tavg = (tmin + tmax) / 2.

                row = cur.execute(select1, (fid, date)).fetchone()
                if row:
                    cur.execute(update1, (p, tmin, tmax, tavg, fid, date))
                else:
                    cur.execute(insert1, (fid, date, p, tmin, tmax, tavg))

        con.commit()
        con.close()


    def swb(self, procedure):
        """
        function where the call to the specific soil water balance function
            is made
        arg
            procedure: swb function name
        """
        from os.path import dirname, join

        metadata_file = 'metadata.csv'
        procedures = ('swb01', 'swb24')
        if procedure not in procedures:
            raise ValueError(f'arg procedure must be in ',
                             f'{",".join(procedures)}')

        con = sqlite3.connect(self.db)
        cur = con.cursor()
        cur.execute("PRAGMA auto_vacuum = FULL")
        cur.execute('attach database ? as allen', ('r0.db',))
        table_output = self.table_output + '_' + procedure
        self._create_output_table(con, cur, table_output)

        cur.execute(self.select_aquifers.text)
        aquifers = [Aquifer(row) for row in cur.fetchall()]
        fmeta = open(join(dirname(self.db), metadata_file), 'w')
        self._write_metadata_base(fmeta, procedure)

        for aquifer in aquifers:
            print(f'{aquifer.name}')

            cur.execute(self.select_outcrops.text, (aquifer.fid,) )
            outcrops = [Outcrop(row) for row in cur.fetchall()]
            cur.execute(self.select_met.text, (aquifer.fid,) )
            met = np.array([row for row in cur.fetchall()])
            dates = met[:, 0]
            imonths = np.array([int(row[5:7])-1 for row in dates], np.int32)
            p = met[:,1].astype(np.float32)
            if procedure == 'swb24':
                m = np.sum(p > 1.)
                gt1 = np.random.randint(1, high=25, size=(m,))
                nh = swb.nhours_generator01(p, gt1)
            et = np.empty((p.size), np.float32)
            if self.proc == 'hargreaves':
                # solar radiation
                cur.execute(BHIMES._select_r0, (int(aquifer.y4326),))
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
                if procedure == 'swb01':
                    xer, irow = swb.swb01(self.time_step, storages, k, p, et,
                                          rch1, runoff1, etr1)
                else:
                    xer, irow, jm = \
                    swb.swb24(outcrop.whc, outcrop.whcr, outcrop.whc*c_whc0,
                              outcrop.kuz, outcrop.kdirect, outcrop.exp,
                              p, nh, et, rch1, runoff1, etr1)
                if irow >= 0:
                    a = (f'Balance error {xer}',  f'Aquifer: {aquifer.name}',
                         f'Outcrop: {outcrop.fid}',
                         f'Date: {dates[irow]} (i {irow})')
                    a = '\n'.join(a)
                    raise ValueError(a)
                self._insert_output(con, cur, aquifer.fid, outcrop.fid, dates,
                                    rch1, runoff1, etr1, procedure)
                rch += rch1  # * (outcrop.area * 0.001)
                runoff += runoff1  # * (outcrop.area * 0.001)
                etr += etr1  # * (outcrop.area * 0.001)
                self._write_metadata(fmeta, procedure, aquifer, outcrop)

        con.close()
        fmeta.close()


    @jit(void(int32[:], float32[:], float32[:]), nopython=True)
    def _et_averaged_set(imonths, et_avg, et):
        """
        sets et arrays whith month averaged values
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


    def _create_output_table(self, con, cur, procedure):
        """
        creates the output table
        parameters
        con: connection -already open-
        cur: cursor to conn -already open-
        procedure: function name where swb is made
        """

        stm1 = f'drop table if exists {self.table_output}'

        stm2 =\
        f"""
        create table if not exists {self.table_output}(
            aquifer integer,
            outcrop integer,
            date text(10),
            rch real,
            runoff real,
            etr real,
            procedure text,
            primary key (aquifer, outcrop, date)
        )
        """
        cur.execute(stm1)
        cur.execute(stm2)
        con.commit()


    def _insert_output(self, con, cur, aquifer, outcrop, dates,
                       rch, runoff, etr, procedure):
        """
        insert or update output in output table
        parameters
        con: connection -already open-
        cur: cursos to conn -already open-
        procedure: function name where swb is made
        dates, np array str: dates (aaaa-mm-dd)
        rch, np array float: calculated recharge mm
        runoff, np array float: calculated runoff mm
        etr, np array float: calculated et mm -real evapotranspiration-
        """

        stm1 =\
        f"""
        insert into {self.table_output}
        (aquifer, outcrop, date, rch, runoff, etr, procedure)
        values (?, ?, ?, ?, ?, ?, ?)
        """

        for i in range(dates.size):
            cur.execute(stm1, (aquifer, outcrop, dates[i], float(rch[i]),
                               float(runoff[i]), float(etr[i]), procedure))
        con.commit()


    def _write_metadata_base(self, fmeta, procedure):
        """
        saves data from selected project
        """
        import datetime

        metadata =\
        (datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
         self.description,
         'swb produre: ' + procedure,
         'et procedure: ' + self.proc,
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
        self._write_header(fmeta, procedure)


    def _write_header(self, fmeta, procedure):
        """
        writes header of metadata
        """
        if procedure == 'swb01':
            a = ('aquifer.name', 'aquifer.y4326',
                 'outcrop.fid', 'outcrop.area',
                 'outcrop.ia', 'outcrop.whc',
                 'outcrop.kdirect', 'outcrop.kuz',
                 'outcrop.klateral', 'outcrop.krunoff'
                )
        else:
            a = ('aquifer.name', 'aquifer.y4326',
                 'outcrop.fid', 'outcrop.area',
                 'outcrop.whc', 'outcrop.whcr',
                 'outcrop.kuz', 'outcrop.exponent'
                )
        fmeta.write(f'{",".join(a)}\n')


    def _write_metadata(self, fmeta, procedure, aquifer, outcrop):
        """
        writes metada
        """
        if procedure == 'swb01':
            a = (f'{aquifer.name}', '{aquifer.y4326}',
                 f'{outcrop.fid:n}', f'{outcrop.area:0.1f}',
                 f'{outcrop.ia:0.2f}', f'{outcrop.whc:0.2f}',
                 f'{outcrop.kdirect:0.2f}', f'{outcrop.kuz:0.2f}',
                 f'{outcrop.klateral:0.2f}', f'{outcrop.krunoff:0.2f}')
        else:
            a = (f'{aquifer.name}', '{aquifer.y4326}',
                 f'{outcrop.fid:n}', f'{outcrop.area:0.1f}',
                 f'{outcrop.whc:0.2f}', f'{outcrop.whcr:0.2f}',
                 f'{outcrop.kuz:0.2f}', f'{outcrop.exp:0.3f}'
                )
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
            cur.execute(BHIMES._select_r0, (int(aquifer.y4326),))
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


    def swb01_sensitivity(self, dir_out: str):
        """
        calls to swb01 using kuz and whc in a defined set of values ranges
        I test recharge, ret and runoff using whc and kuz in a range of
            values defined in self.params
        As an aquifer can have several outcrops, I average the outcrops
            parameter by area; this way, I just make a call to the function
            swb01 for each aquifer
        The outputs area saved in dir_out in text files: one has the time
            series and the other a data summary of each simulation; additional
            png files are saved (p, runoff and ret) considerig units are
            in mm/day
        """
        from os.path import join

        con = sqlite3.connect(self.db)
        cur = con.cursor()
        cur.execute("PRAGMA auto_vacuum = FULL")
        cur.execute('attach database ? as allen', ('r0.db',))

        cur.execute(self.select_aquifers.text)
        aquifers = [Aquifer(row) for row in cur.fetchall()]

        self._create_tables_ts(cur)

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
                cur.execute(BHIMES._select_r0, (int(aquifer.y4326),))
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

            sum_outcrops_area = Outcrop.sum_area(outcrops)
            coefs = np.array([outcrop.area for outcrop in outcrops],
                             np.float32) / sum_outcrops_area
            pars = {'ia': 0., 'whc': 0., 'kdirect': 0., 'kuz': 0.,
                    'klateral': 0., 'krunoff': 0.}
            for key in pars:
                pars[key] = BHIMES._averaged_parameter(key, coefs, outcrops)

            c_ia0, c_whc0 = self._coef_initial_wstorages()

            storages = np.array([pars['ia'], pars['ia']*c_ia0,
                                 pars['whc'], pars['whc']*c_whc0], np.float32)
            k = np.array([pars['kdirect'], pars['kuz'],
                          pars['klateral'], pars['krunoff']], np.float32)

            n = 0
            xf = p.size / 365.
            contour = np.empty((self.neval('kuz') * self.neval('whc'), 5),
                               np.float32)
            for i in range(self.neval('kuz')):
                xk = np.copy(k)
                xk[1] = pars['kuz'] + (pars['kuz'] * self.delta('kuz') * i)
                for j in range(self.neval('whc')):
                    n += 1
                    print(f'{n:n}')
                    xstorages = np.copy(storages)
                    xstorages[2] = pars['whc'] + (pars['whc'] * \
                                                  self.delta('whc') * j)

                    ier, xer = swb.swb01(self.time_step, xstorages, xk, p, et,
                                         rch, runoff, etr)
                    if ier >= 0:
                        a = (f'Balance error {xer}',
                             f'Aquifer: {aquifer.name}',
                             f'Date: {dates[ier]} (i {ier})')
                        a = '\n'.join(a)
                        raise ValueError(a)

                    contour[n-1][:] = [xstorages[2], xk[1], rch.sum(),
                                 runoff.sum(), etr.sum()]
                    self._insert_iter_ts(cur, n, self.time_step,
                                         xstorages[0], xstorages[1],
                                         xstorages[2], xstorages[3],
                                         xk[0], xk[1], xk[2], xk[3])

                    self._insert_ts(cur, aquifer.fid, n, dates, rch, runoff,
                                    etr)

            for i in range(2, 5):
                contour[:,i] = contour[:,i] * (sum_outcrops_area * 0.001 / xf)
            for i, item in enumerate(('recharge', 'runoff', 'etr')):
                BHIMES._contour(f'{aquifer.name}: {item} m3/yr',
                                contour[:, 0], contour[:, 1],
                                contour[:, i+2], 'whc mm', 'kuz mm',
                                join(dir_out,
                                     f'aq_{aquifer.fid}_sensitivity_{item}'))

        con.commit()
        con.close()


    def swb24_sensitivity(self, dir_out: str):
        """
        calls to swb24 using exp, kuz and whc in a defined set of values ranges
        I test recharge, ret and runoff using exp, whc and kuz in a range of
            values defined in self.params
        As an aquifer can have several outcrops, I average the outcrops
            parameter by area; this way, I just make a call to the function
            swb01 for each aquifer
        The outputs area saved in dir_out in text files: one has the time
            series and the other a data summary of each simulation; additional
            png files are saved (p, runoff and ret) considerig units are
            in mm/day
        """
        from os.path import join

        con = sqlite3.connect(self.db)
        cur = con.cursor()
        cur.execute("PRAGMA auto_vacuum = FULL")
        cur.execute('attach database ? as allen', ('r0.db',))

        cur.execute(self.select_aquifers.text)
        aquifers = [Aquifer(row) for row in cur.fetchall()]

        self._create_tables_ts(cur, 'swb24')

        for aquifer in aquifers:
            print(f'{aquifer.name}')

            cur.execute(self.select_outcrops.text, (aquifer.fid,) )
            outcrops = [Outcrop(row) for row in cur.fetchall()]
            cur.execute(self.select_met.text, (aquifer.fid,) )
            met = np.array([row for row in cur.fetchall()])
            dates = met[:, 0]
            imonths = np.array([int(row[5:7])-1 for row in dates], np.int32)
            p = met[:,1].astype(np.float32)
            m = np.sum(p > 1.)
            gt1 = np.random.randint(1, high=25, size=(m,))
            nh = swb.nhours_generator01(p, gt1)
            nyears = p.size / 365.
            et = np.empty((p.size), np.float32)
            if self.proc == 'hargreaves':
                # solar radiation
                cur.execute(BHIMES._select_r0, (int(aquifer.y4326),))
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

            sum_outcrops_area = Outcrop.sum_area(outcrops)
            coefs = np.array([outcrop.area for outcrop in outcrops],
                             np.float32) / sum_outcrops_area
            pars = {'whc': 0., 'whcr': 0., 'kdirect': 0., 'kuz': 0.,
                    'bcexp': 0.}
            for key in pars:
                pars[key] = BHIMES._averaged_parameter_24(key, coefs, outcrops)

            c_ia0, c_whc0 = self._coef_initial_wstorages()

            n = 0
            contour = np.empty((self.neval('kuz') * self.neval('whc'), 5),
                               np.float32)
            for i in range(self.neval('exp')):
                xexp = pars['bcexp'] + (pars['bcexp'] * self.delta('exp') * i)
                ic = -1
                for j in range(self.neval('kuz')):
                    xkuz = pars['kuz'] + (pars['kuz'] * self.delta('kuz') * j)
                    for k in range(self.neval('whc')):
                        xwhc = \
                        pars['whc'] + (pars['whc'] * self.delta('whc') * k)
                        n += 1
                        print(f'{n:n}')

                        xer, irow, jh = \
                        swb.swb24(xwhc, pars['whcr'], xwhc*c_whc0,
                                  xkuz, pars['kdirect'], xexp, p, nh, et, rch,
                                  runoff, etr)

                        if irow >= 0:
                            a = (f'Balance error {xer}',
                                 f'Aquifer: {aquifer.name}',
                                 f'Date: {dates[irow]} (i {irow:n}, j {jh:n})')
                            a = '\n'.join(a)
                            raise ValueError(a)
                        ic += 1
                        contour[ic][:] = [xwhc, xkuz, rch.sum()/nyears,
                                          runoff.sum()/nyears,
                                          etr.sum()/nyears]
                        self._insert_iter_ts24(cur, ic, xwhc, pars['whcr'],
                                               xwhc*c_whc0, pars['kdirect'],
                                               xkuz, xexp)
                        self._insert_ts(cur, aquifer.fid, n, dates,
                                        rch, runoff, etr)

                for i, item in enumerate(('drainage', 'runoff', 'etr')):
                    title = f'Aquifer {aquifer.name}: {item} mm/yr, ' +\
                    f'B-C exp {xexp:0.3f}'
                    namefile = f'aq_{aquifer.fid}_sensitivity_{item}_{i:n}'
                    BHIMES._contour(title, contour[:, 0], contour[:, 1],
                                    contour[:, i+2], 'whc mm', 'kuz mm',
                                    join(dir_out, namefile))

        con.commit()
        con.close()


    @staticmethod
    def _coef_area_get(outcrops):
        """
        coefficients by area
        """
        x = np.array([outcrop.area for outcrop in outcrops], np.float32)
        y = x[:] / x.sum()
        return y


    @staticmethod
    def _averaged_parameter(param, coefs, outcrops):
        """
        averaged parameter by area
        parameter are different in swb01 and swb24
        """
        if param == 'ia':
            x = np.array([outcrop.ia * coefs[i] for i,
                          outcrop in enumerate(outcrops)], np.float32)
        elif param == 'whc':
            x = np.array([outcrop.whc * coefs[i] for i,
                          outcrop in enumerate(outcrops)], np.float32)
        elif param == 'whcr':
            x = np.array([outcrop.whcr * coefs[i] for i,
                          outcrop in enumerate(outcrops)], np.float32)
        elif param == 'kdirect':
            x = np.array([outcrop.kdirect * coefs[i] for i,
                          outcrop in enumerate(outcrops)],
                        np.float32)
        elif param == 'kuz':
            x = np.array([outcrop.kuz * coefs[i] for i,
                          outcrop in enumerate(outcrops)], np.float32)
        elif param == 'klateral':
            x = np.array([outcrop.klateral * coefs[i] for i,
                          outcrop in enumerate(outcrops)], np.float32)
        elif param == 'krunoff':
            x = np.array([outcrop.krunoff * coefs[i] for i,
                          outcrop in enumerate(outcrops)], np.float32)
        elif param == 'bcexp':
            x = np.array([outcrop.exp * coefs[i] for i,
                          outcrop in enumerate(outcrops)], np.float32)
        else:
            raise ValueError(f'{param} no es un parámetro válido')
        return x.sum()


    @staticmethod
    def _contour(title, x, y, z, xlabel, ylabel, dst, scale: float=1.0):
        """
        3D representarion of the results
        """
        import matplotlib.pyplot as plt

        fig, ax1 = plt.subplots(figsize=(8,6))

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


    def _create_tables_ts(self, cur, proc = 'swb01'):
        """
        saves time series outputs to table swb01_ts
        """

        drop_table0 = f'drop table if exists {BHIMES._tname_iters_ts[proc]}'
        drop_table1 = f'drop table if exists {BHIMES._tname_ts[proc]}'

        if proc == 'swb01':
            create_table0 = \
            f"""
            create table if not exists {BHIMES._tname_iters_ts[proc]} (
                n integer,
                ntimestep integer,
                iamax real,
                ia0 real,
                whcmax real,
                whc0 real,
                kdirect real,
                kuz real,
                klateral real,
                krunoff real,
                primary key (n)
            )
            """
        else:
            create_table0 = \
            f"""
            create table if not exists {BHIMES._tname_iters_ts[proc]} (
                n integer,
                whc real,
                whcr real,
                whc0 real,
                kdirect real,
                kuz real,
                bcexp real,
                primary key (n)
            )
            """
        create_table1 = \
        f"""
        create table if not exists {BHIMES._tname_ts[proc]} (
            aquifer integer,
            n integer,
            date text,
            rch real,
            rnf real,
            etr real,
            primary key(aquifer, n, date)
        )
        """
        for item in (drop_table0, drop_table1, create_table0, create_table1):
            cur.execute(item)


    def _insert_iter_ts(self, cur, n, ntimestep, iamax, ia0, whcmax, whc0,
                        kdirect, kuz, klateral, krunoff):
        proc = 'swb01'
        insert = \
        f"""
        insert into {BHIMES._tname_iters_ts[proc]} (n, ntimestep, iamax, ia0,
            whcmax, whc0, kdirect, kuz, klateral, krunoff)
        values (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """
        cur.execute(insert, (n, ntimestep, float(iamax), float(ia0),
                             float(whcmax), float(whc0),
                             float(kdirect), float(kuz), float(klateral),
                             float(krunoff)))


    def _insert_iter_ts24(self, cur, n, whc, whcr, whc0, kdirect, kuz, bcexp):
        proc = 'swb24'
        insert = \
        f"""
        insert into {BHIMES._tname_iters_ts[proc]} (n, whc, whcr, whc0,
            kdirect, kuz, bcexp)
        values (?, ?, ?, ?, ?, ?, ?)
        """
        cur.execute(insert, (n, float(whc), float(whcr), float(whc0),
                             float(kdirect), float(kuz), float(bcexp)))


    def _insert_ts(self, cur, aquifer_fid, n, dates, rch, rnf, etr):

        insert = f"""
        insert into {BHIMES._tname_ts}
        (aquifer, n, date, rch, rnf, etr)
        values (?, ?, ?, ?, ?, ?)
        """
        for i in range(rch.size):
            cur.execute(insert, (aquifer_fid, n, dates[i], float(rch[i]),
                                 float(rnf[i]), float(etr[i])))


class Aquifer():
    """
    acuíferos
    """


    def __init__(self, data: list):
        """
        fid integer: identifier
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
    permeable outcrops
    currently I can call 2 swb function with sligthly number of outcrops
        parameters. The outcrop class has the parameters of both
        functions
    """


    def __init__(self, data: list):
        """
        fid: identifier
        area: outcrop area m2
        ia: initial abstraction mm
        whc: mas water holding content mm
        kdirect: k recharge mm
        kuz: k unsatured zone mm
        klateral: k lateral applied to kuz mm
        krunof: k from runoff to kuz mm
        whcr: residual whc -whc > whcr- mm
        exp: exponent in the relation ((whc1-whcr)/(whc/whcr))**exp
        """
        self.data = tuple([row for row in data])
        if self.area <= 0:
            raise ValueError(f'outcrop {self.fid}, area <=0')


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


    @property
    def whcr(self):
        return self.data[8]


    @property
    def exp(self):
        return self.data[8]


    @staticmethod
    def sum_area(outcrops: list):
        return np.array([outcrop.area for outcrop in outcrops],
                        np.float32).sum()

