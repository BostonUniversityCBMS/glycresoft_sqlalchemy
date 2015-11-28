from glycresoft_sqlalchemy.data_model import *
import sqlite3
import uuid


def main(conn):
    models = Base._decl_class_registry.values()
    for model in [SampleRun]:
            tn = model.__tablename__
            try:
                print(list(conn.execute("select uuid from %s;" % tn)))
            except:
                conn.execute("ALTER TABLE %s ADD COLUMN uuid CHAR(32);" % tn)
                for id in conn.execute("SELECT id FROM %s;" % tn):
                    id = id[0]
                    conn.execute("UPDATE %s SET uuid = ? WHERE id = ?" % tn, (uuid.uuid4().hex, id))
    conn.commit()

if __name__ == '__main__':
    main(sqlite3.connect(sys.argv[1]))
