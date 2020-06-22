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

import littleLogging as logging


class BHIMES():
    """
    Balance HIdroMEteorológico en el Suelo
    Soil water balance
    """

    def __init__(self):
        pass


    def create_db(self, dst: str):
        """
        create a sqlite db
        """
        import sqlite3

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
            fid integer primary key,
            date text(10)
            p real,
            tmin real,
            tmax real,
            tavg real,
            foreign key (fid) references aquifer(fid)
        )
        """

        con = sqlite3.connect(dst)
        con.execute("PRAGMA encoding = 'utf-8'")
        con.execute("PRAGMA foreign_keys = 1")
        cur = con.cursor()
        for stm in (stm1, stm2, stm3, stm4):
            cur.execute(stm)
        con.commit()
        con.close()


    def aquifer_upsert_from_file(self, db: str, org: str,
                                 nskip: int=1, separator: str=';'):
        """
        inserts or update aquifer data from text file -utf8-
        """
        import sqlite3

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

        con = sqlite3.connect(db)
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


    def outcrop_upsert_from_file(self, db: str, org: str,
                                 nskip: int=1, separator: str=';'):
        """
        inserts or update outcrops data from text file -utf8-
        """
        import sqlite3

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

        con = sqlite3.connect(db)
        cur = con.cursor()

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


    def met_upsert_from_file01(self, db: str, org: str,
                                 nskip: int=1, separator: str=';'):
        """
        inserts or update meteorological data from text file -utf8-
        with no average temperature
        """
        import sqlite3

        select1 = \
        """
        select fid from met where fid=?
        """
        update1 = \
        """
        update met set date=?, p=?, tmin=?, tmax=?, tavg=?
        where fid=?
        """
        insert1 = \
        """
        insert into aquifer(fid, date, p, tmin, tmax, tavg)
        values (?, ?, ?, ?, ?, ?)
        """

        con = sqlite3.connect(db)
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

                row = cur.execute(select1, (fid,)).fetchone()
                if row:
                    cur.execute(update1, (date, p, tmin, tmax, tavg, fid))
                else:
                    cur.execute(insert1, (fid, date, p, tmin, tmax, tavg))

        con.commit()
        con.close()


    def swb01(self, aquifers: list, table_output: str='swb01'):
        """
        soil water balance version 0.01
        """
        select1 = \
        """

        """