import os
import sqlite3
import csv

from .log import log
from .utils import load_csv_data, write_csv_dict

__author__ = 'Benjamin Leopold <bleopold@jax.org>'


class SqliteIO():
    """Library of Sqlite DB I/O methods"""

    @staticmethod
    def connect(dbname, row_dict=True):
        """Connect to sqlite db, using Row in factory"""
        try:
            log.info('Connecting to sqlite db: {}'.format(dbname))
            dbkws = {}
            dbkws['detect_types'] = sqlite3.PARSE_DECLTYPES
            con = sqlite3.connect(dbname, **dbkws)
            if row_dict:
                con.row_factory = SqliteIO._dict_row_factory
            else:
                con.row_factory = sqlite3.Row
        except sqlite3.Error as e:
            log.error('Connect to db "{}": {}'.format(dbname, e))
            raise e
        except Exception as e:
            log.error('Connect to db "{}": {}'.format(dbname, e))
            raise e
        else:
            return con


    def _dict_row_factory(cursor, row):
        """Return row from cursor select as dict. {column_name: value}
        Compared to sqlite3.Row (tuple), this dict gives name-access plus mutability.
        """
        d = {}
        for idx, col in enumerate(cursor.description):
            d[col[0]] = row[idx]
        return d


    @staticmethod
    def iter_select(dbname, table, fields=None, where=None, row_dict=True):
        """Iterate on selected rows from table in dbname"""
        try:
            log.info('Selecting data from {}'.format(dbname))
            if isinstance(fields, list):
                field_def = ', '.join(fields)
            elif isinstance(fields, str):
                field_def = fields
            else:
                field_def = '*'

            where_def = where or '1=1'

            select_sql = f'SELECT {field_def} FROM {table} WHERE {where_def};'
            # select_sql = 'SELECT {} FROM {} WHERE {}'.format(field_def, table, where_def) 

            with SqliteIO.connect(dbname, row_dict=row_dict) as db:
                log.debug('Executing: "{}"'.format(select_sql))
                for row in db.execute(select_sql):
                    yield row
        except sqlite3.Error as e:
            log.error('Selecting with "{}" in db: {}\n{}'.format(select_sql, dbname, e))
            raise e
        except Exception as e:
            log.error('Selecting with "{}" in db: {}\n{}'.format(select_sql, dbname, e))
            raise e
        else:
            db.close()


    @staticmethod
    def exec_ddl(dbname, ddl_sql):
        """Create object in db using DDL sql statement"""
        try:
            if not ddl_sql.endswith(';'):
                ddl_sql += ';'
            if sqlite3.complete_statement(ddl_sql):
                with SqliteIO.connect(dbname, row_dict=False) as db:
                    log.debug('Executing: "{}"'.format(ddl_sql))
                    db.execute(ddl_sql)
                    return True
            else:
                log.error('Can''t execute incomplete sql: "{}"'.format(ddl_sql))
                return False
        except sqlite3.Error as e:
            log.error('Defining object "{}" in db: {}\n{}'.format(ddl_sql, dbname, e))
            raise e
        except Exception as e:
            log.error('Defining object "{}" in db: {}\n{}'.format(ddl_sql, dbname, e))
            raise e
        else:
            db.close()


    @staticmethod
    def import_data(data_list, dbname, table):
        """Insert rows from data list of dicts into database table.
                Fieldnames are data_dict[0].keys()
                Values are from data_dict[:].values()
        """
        try: # check data_list is expected format!
            err_msg = "'data_list' must be a list of dict's of field names and values!"
            if type(data_list) != list and type(data_list[0]) != dict:
                log.error(err_msg)
                return err_msg
        except IndexError as e:
            log.error(err_msg+'\n{}'.format(e))
            raise e
        except Exception as e:
            log.error(err_msg+'\n{}'.format(e))
            raise e

        try:
            fields = list(data_list[0].keys())
            with SqliteIO.connect(dbname) as db:
                dbcur = db.cursor()
                sql_cols = ','.join('?' * len(fields))
                fieldnames = ','.join(fields)
                sql_insert = 'INSERT INTO {} ({}) VALUES ({});' \
                             ''.format(table, fieldnames, sql_cols) 

                for row_num, row in enumerate(data_list, start=1):
                    dbcur.execute( sql_insert, ( list(row.values()) ) )
                    # log.debug(' ... Inserted row {} into table "{}"'.format(row_num, table))
                log.info('Inserted {} rows into table "{}"'.format(dbcur.lastrowid, table))
            log.info('Import session complete.')
        except sqlite3.Error as e:
            log.error('Importing into db "{}": {}'.format(dbname, e))
            raise e
        except Exception as e:
            log.error('Importing into db "{}": {}'.format(dbname, e))
            raise e


    #TODO: db export arg: filetype=csv/tsv or plain text!
    @staticmethod
    def export_csv(dbname, table, filepath, fields=None, where=None, delimiter=','):
        """Export from sqlite database to csv file (created or appended).
        Specify list or string of field names for export, defaults to all (select *).
        Specify where clause to select, defaults to all.
        """
        try:
            log.error('Exporting to csv file "{}"'.format(filepath))
            if os.path.isfile(filepath):
                selects = SqliteIO.iter_select(dbname, table, fields, where) 
                records = [ row for row in selects ]
                return write_csv_dict(filepath, values=records, delim=delimiter)
            else:
                log.error('Exporting to csv file "{}" is not file...'.format(filepath))
                return None
        except Exception as e:
            log.error('Exporting from db "{}" to file "{}": {}'.format(dbname, filepath, e))
            raise e


    @staticmethod
    def get_csv_field_datatypes(filename, delim=','):
        """guess datatypes from data in csv file.
        Return dict of {field: typename, ...}
        """
        try:
            if os.path.isfile(filename):
                log.info('Getting types of fields from {}'.format(filename))
                ftypes = {}
                with open(filename, 'r') as csvfh:
                    dreader = csv.DictReader(csvfh, delimiter=delim)
                    fields = dreader.fieldnames
                    data = next(dreader)

                for f in fields:
                    if not data[f]:
                        log.notice('field: {}'.format(f))
                        continue
                    if data[f].isdigit():
                        ftypes[f] = 'INTEGER'
                    elif data[f].replace('.','').isdigit():
                        ftypes[f] = 'REAL'
                    else:
                        ftypes[f] = 'TEXT'
                return ftypes
            else:
                log.warning('Fieldtypes: "{}" not a file?'.format(filename))
                return {}
        except Exception as e:
            log.error('Getting fieldtypes from "{}"; {}'.format(filename, e))
            raise e


    #FIXME: parsing errors in field datatypes and list-length mismatches; skip dataypes?
    @staticmethod
    def import_csv(filename, dbname, table=None, fields=[],
                     delim=',', dialect='unix'):
        """Import csv file to sqlite database (created if not found),
        into a specific table if name is passed, or will create new
        table with same name as file.
        Column names will be taken from first row in file if fields[list] not passed.
        """
        try:
            with SqliteIO.connect(dbname) as db:
                if not table:
                    filebase = filename
                    if filename.rfind('.'): # e.g. '.csv'
                        filebase = filename[0:filename.rfind('.')]
                    table = filebase

                if os.path.isfile(filename):
                    field_datatypes = SqliteIO.get_csv_field_datatypes(filename, delim=delim)
                    log.notice('length fields: {}'.format(len(fields)))
                    if not fields:
                        fields = field_datatypes.keys()
                        log.notice('length fields: {}'.format(len(fields)))

                    if len(field_datatypes):
                        fld_types = {}
                        for fld in fields:
                            fld_types[fld] = field_datatypes.pop(fld, '')
                        # include leftover fields in file:
                        fld_types.update(field_datatypes)

                        col_defs = ', '.join([' '.join(t) for t in fld_types])
                    else:
                        col_defs = ', '.join(fields)
                    log.notice('col_defs: {}'.format(col_defs))

                    try:
                        sql_cols = ','.join('?' * len(fields))
                        fieldnames = ','.join(fields)
                        sql_insert = 'INSERT INTO {} ({}) VALUES ({});' \
                                     ''.format(table, fieldnames, sql_cols) 
                        dbcur = db.cursor()

                        rows = load_csv_data(filename, fields=fields, skip_rows=1,
                                         delim=delim, dialect=dialect)
                        for row in rows:
                            row_values = list(row.values())
                            dbcur.execute( sql_insert, (row_values) )
                            # log.notice('Inserted row {} into table "{}"'.format(dbcur.lastrowid, table))
                        log.notice('Inserted {} rows into table "{}"'.format(dbcur.lastrowid, table))
                    except Exception as e:
                        log.error('Importing records from file "{}": {}'.format(filename, e))
                        raise e

                    return # Success!!
                else:
                    log.error('Reading file "{}" has problems...'.format(filename))
                    return None
        except sqlite3.Error as e:
            log.error('Importing from file "{}" to db: {}'.format(filename, e))
            raise e
        except Exception as e:
            log.error('Importing from file "{}" to db: {}'.format(filename, e))
            raise e
        finally:
            db.close()
