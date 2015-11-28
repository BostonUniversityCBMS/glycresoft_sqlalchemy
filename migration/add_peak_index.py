from glycresoft_sqlalchemy.data_model import *
import sqlite3


def main(conn):
    models = Base._decl_class_registry.values()
    for model in [Peak, Decon2LSPeak]:
            tn = model.__tablename__
            try:
                conn.execute("select scan_peak_index from %s limit 1;" % tn)
            except:
                conn.execute("ALTER TABLE %s ADD COLUMN scan_peak_index INTEGER;" % tn)
                conn.execute("UPDATE %s SET scan_peak_index = rowid;" % tn)
    conn.commit()

if __name__ == '__main__':
    main(sqlite3.connect(sys.argv[1]))
