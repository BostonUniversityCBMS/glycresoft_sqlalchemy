"""Add theoretical_glycopeptide_id to GlycopeptideMatch

Revision ID: 76205281ff5
Revises: 3f08a8fbe5c2
Create Date: 2015-08-05 13:21:48.526000

"""

# revision identifiers, used by Alembic.
revision = '76205281ff5'
down_revision = '3f08a8fbe5c2'
branch_labels = None
depends_on = None

from alembic import op
import sqlalchemy as sa


def upgrade():
    with op.batch_alter_table("GlycopeptideMatch") as batch_op:
        batch_op.add_column(sa.Column("theoretical_glycopeptide_id", sa.Integer))
        batch_op.create_foreign_key(
            "fk_GlycopeptideMatch_base_TheoreticalGlycopeptide",
            "TheoreticalGlycopeptide", ["theoretical_glycopeptide_id"], ["id"])



def downgrade():
    pass
