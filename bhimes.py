# -*- coding: utf-8 -*-
"""
Created on Sun Jun 21 12:51:28 2020

@author: solis
"""
import littleLogging as logging


class BHIMES():
    """
    Balance HIdroMEteorológico en el Suelo
    Water balance in soil
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
            name text(100),
            sup_m2 real,
            primary key (fid)
        )
        """

        stm2 = \
        """
        create table if not exists outcrop(
            fid integer,
            aquifer integer,
            lito text(250),
            period text(100),
            perme text(50),
            sup_m2 real,
            ia real,
            whc real,
            kdirect real,
            kuz real,
            klateral real,
            krunoff real,
            primary key(fid)
        );
        references aquifer(fid)
            on update cascade
            on delete cascade;
        create index outcrop_aquifer
            on outcrop(aquifer);
        """

        con = sqlite3.connect(dst)
        cur = con.cursor()
        cur.execute(stm1)
        cur.execute(stm2)
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
        select * from aquifer where fid=?
        """
        update1 = \
        """
        update aquifer set name=?, sup_m2=?
        where fid=?
        """
        insert1 = \
        """
        insert into aquifer(fid, name, sup_m2)
        values (?, ?, ?)
        """

        con = sqlite3.connect(db)
        cur = con.cursor()

        with open(org, 'r') as f:
            for i, line in enumerate(f):
                if i < nskip:
                    continue
                words = line.strip().split(separator)
                fid = int(words[0])
                name = words[1]
                sup = float(words[2])
                row = cur.execute(select1, (fid,)).fetchone()
                if row:
                    cur.execute(update1, (name, sup, fid))
                else:
                    cur.execute(insert1, (fid, name, sup))

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
        select * from outcrop where fid=?
        """
        update1 = \
        """
        update outcrop set aquifer=?, lito=?, period=?, perme=?, sup_m2=?,
        ia=?, whc=?, kdirect=?, kuz=?, klateral=?, krunoff=?
        where fid=?
        """
        insert1 = \
        """
        insert into outcrop(fid, aquifer, lito, period, perme, sup_m2, ia, whc,
            kdirect, kuz, klateral, krunoff)
        values (?,?,?,?,?,?,?,?,?,?,?,?,?)
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
                sup_m2 = float(words[5])
                ia = float(words[6])
                whc = float(words[7])
                kdirect = float(words[8])
                kuz = float(words[9])
                klateral = float(words[10])
                krunoff = float(words[11])
                row = cur.execute(select1, (fid,)).fetchone()
                if row:
                    cur.execute(update1, (acu, lito, period, perme,
                                          sup_m2, ia, whc, kdirect, kuz,
                                          klateral, krunoff, fid))
                else:
                    cur.execute(insert1, (fid, acu, lito, period, perme,
                                          sup_m2, ia, whc, kdirect, kuz,
                                          klateral, krunoff))

        con.commit()
        con.close()
