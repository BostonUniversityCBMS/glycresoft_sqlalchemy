"""Add association from TheoreticalGlycopeptide to PeakGroupMatch

Revision ID: 3f08a8fbe5c2
Revises: 
Create Date: 2015-08-05 13:01:53.014000

"""

# revision identifiers, used by Alembic.
revision = '3f08a8fbe5c2'
down_revision = None
branch_labels = None
depends_on = None

from alembic import op
import sqlalchemy as sa


def upgrade():
    with op.batch_alter_table("TheoreticalGlycopeptide") as batch_op:
        batch_op.add_column(sa.Column("base_composition_id", sa.Integer))
        batch_op.create_foreign_key(
            "fk_TheoreticalGlycopeptide_base_PeakGroupMatch",
            "PeakGroupMatch", ["base_composition_id"], ["id"])


def downgrade():
    pass
