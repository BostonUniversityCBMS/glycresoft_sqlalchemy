from sqlalchemy.orm import make_transient
import sqlitedict

Incomplete = object()


class Migrator(object):
    def __init__(self, source, target):
        self.source = source
        self.target = target
        self.copy_register = sqlitedict.open()

    def find_foreign_keys(self, model):
        fks = []
        inherited_tables = set()
        for layer in model.mro():
            if not hasattr(layer, "__table__"):
                break
            inherited_tables.add(layer.__table__.name)
            for fk in layer.__table__.foreign_keys:
                fks.append(fk)
        fks = [fk for fk in fks if fk.column.table.name not in inherited_tables]
        return fks

    def look_up_reference(self, fk, value):
        if value is None:
            return None
        new_id = self.copy_register.get(str((fk.column.table.name, value)), Incomplete)
        if new_id is Incomplete:
            raise ValueError("Relation not yet copied", fk, value)
        return new_id

    def look_up_reference_for_instance(self, instance):
        for klass in instance.__class__.mro():
            if not hasattr(klass, "__tablename__"):
                continue
            new_id = self.copy_register.get(str((instance.__tablename__, instance.id)), Incomplete)
            if new_id is Incomplete:
                continue
            return new_id
        raise ValueError("Relation not yet copied", instance)

    def set_reference(self, model, old_value, new_value):
        for layer in model.mro():
            if not hasattr(layer, "__table__"):
                break
            self.copy_register[str((layer.__table__.name, old_value))] = new_value

    def copy_model(self, model, filterfunc=lambda q: q, batch_size=10000):
        fks = self.find_foreign_keys(model)
        i = 0
        with self.source.no_autoflush:
            for ref in filterfunc(self.source.query(model)).yield_per(1000):
                make_transient(ref)
                for fk in fks:
                    value = getattr(ref, fk.parent.name)
                    setattr(ref, fk.parent.name, self.look_up_reference(fk, value))
                old_id = ref.id
                ref.id = None
                self.target.add(ref)
                self.target.flush()
                assert ref.id is not None
                self.set_reference(model, old_id, ref.id)
                i += 1
                if i % batch_size == 0:
                    self.target.commit()
        self.target.commit()


class MigrationPipeline(object):
    def __init__(self, source, target, order=tuple(), batch_size=10000):
        self.source = source
        self.target = target
        self.order = list(order)
        self.batch_size = batch_size

    def add(self, model, filterfunc=lambda q: q):
        self.order.append((model, filterfunc))

    def run(self):
        migrator = Migrator(self.source, self.target)
        for model, filterfunc in self.order:
            migrator.copy_model(model, filterfunc, self.batch_size)
