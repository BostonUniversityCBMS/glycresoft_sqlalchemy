CREATE TABLE alembic_version (
    version_num VARCHAR(32) NOT NULL
);

-- Running upgrade  -> 3f08a8fbe5c2

ALTER TABLE "TheoreticalGlycopeptide" ADD COLUMN base_composition_id INTEGER REFERENCES PeakGroupMatch(id);

INSERT INTO alembic_version (version_num) VALUES ('3f08a8fbe5c2');

-- Running upgrade 3f08a8fbe5c2 -> 76205281ff5

ALTER TABLE "GlycopeptideMatch" ADD COLUMN theoretical_glycopeptide_id INTEGER REFERENCES TheoreticalGlycopeptide(id);

UPDATE alembic_version SET version_num='76205281ff5' WHERE alembic_version.version_num = '3f08a8fbe5c2';

