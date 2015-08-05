CREATE TABLE alembic_version (
    version_num VARCHAR(32) NOT NULL
);

-- Running upgrade  -> 3f08a8fbe5c2

ALTER TABLE "TheoreticalGlycopeptide" ADD COLUMN base_composition_id INTEGER References PeakGroupMatch(id);

INSERT INTO alembic_version (version_num) VALUES ('3f08a8fbe5c2');

