from glycresoft_sqlalchemy.data_model import *
import sqlite3


def main(conn):
    models = Base._decl_class_registry.values()
    for model in models:
        if hasattr(model, "GlycanCompositionAssociation"):
            tn = model.GlycanCompositionAssociation.__tablename__
            try:
                conn.execute("select id from %s limit 1;" % tn)
            except:
                conn.execute("ALTER TABLE %s ADD COLUMN id INTEGER;" % tn)
                conn.execute("UPDATE %s SET id = rowid;" % tn)
    conn.commit()

if __name__ == '__main__':
    main(sqlite3.connect(sys.argv[1]))
