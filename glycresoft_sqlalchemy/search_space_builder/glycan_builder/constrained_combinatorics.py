from glycresoft_sqlalchemy.data_model import MS1GlycanHypothesis, TheoreticalGlycanComposition, PipelineModule
from glycresoft_sqlalchemy.utils.database_utils import get_or_create
from glypy import GlycanComposition, MonosaccharideResidue, monosaccharides

import logging
from itertools import product

logger = logging.getLogger("glycan_composition_constrained_combinatorics")


from glycresoft_sqlalchemy.search_space_builder.glycan_builder import registry


@registry.composition_source_type.register("constrained_combinatorics")
class ConstrainedCombinatoricsGlycanHypothesisBuilder(PipelineModule):
    HypothesisType = MS1GlycanHypothesis

    def __init__(self, database_path, rules_table, constraints_list, hypothesis_id=None):
        self.manager = self.manager_type(database_path)
        self.rules_table = rules_table
        self.constraints_list = constraints_list
        self.hypothesis_id = hypothesis_id

    def run(self):
        session = self.manager.session()

        hypothesis, _ = get_or_create(session, self.HypothesisType, id=self.hypothesis_id)
        session.add(hypothesis)
        session.commit()

        hypothesis_id = self.hypothesis_id = hypothesis.id

        generator = CombinatoricCompositionGenerator(
            rules_table=self.rules_table, constraints=self.constraints_list)

        acc = []
        for composition in generator:
            mass = composition.mass()
            serialized = composition.serialize()
            rec = TheoreticalGlycanComposition(
                calculated_mass=mass, composition=serialized,
                hypothesis_id=hypothesis_id)
            acc.append(rec)

            if len(acc) > 1000:
                session.add_all(acc)
                session.commit()
                acc = []

        session.close()


def descending_combination_counter(counter):
    keys = counter.keys()
    count_ranges = map(lambda lo_hi: range(lo_hi[0], lo_hi[1] + 1), counter.values())
    for combination in product(*count_ranges):
        yield dict(zip(keys, combination))


class CombinatoricCompositionGenerator(object):
    def __init__(self, residue_list=None, lower_bound=None, upper_bound=None, constraints=None, rules_table=None):
        self.residue_list = residue_list or []
        self.lower_bound = lower_bound or []
        self.upper_bound = upper_bound or []
        self.constraints = constraints or []
        self.rules_table = rules_table

        if len(self.constraints) > 0 and not isinstance(self.constraints[0], CompositionConstraint):
            self.constraints = map(CompositionConstraint.from_list, self.constraints)

        if rules_table is None:
            self._build_rules_table()

    def _build_rules_table(self):
        rules_table = {}
        for i, residue in enumerate(self.residue_list):
            lower = self.lower_bound[i]
            upper = self.upper_bound[i]
            rules_table[residue] = (lower, upper)
        self.rules_table = rules_table
        return rules_table

    def generate(self):
        for combin in descending_combination_counter(self.rules_table):
            passed = True
            for constraint in self.constraints:
                if not constraint(combin):
                    passed = False
                    break
            if passed:
                yield GlycanComposition(**combin)

    __iter__ = generate

    def __repr__(self):
        return repr(self.rules_table) + '\n' + repr(self.constraints)


class CompositionConstraint(object):

    @classmethod
    def from_list(cls, sym_list):
        assert len(sym_list) == 3
        left, op, right = sym_list
        left = SymbolNode.parse(left)
        right = SymbolNode.parse(right)
        op = Operator.operator_map[op]
        return CompositionConstraint(left, op, right)

    def __init__(self, left, operator, right):
        self.left = left
        self.operator = operator
        self.right = right

    def __repr__(self):
        return "{} {} {}".format(self.left, self.operator, self.right)

    def __call__(self, context):
        return self.operator(self.left, self.right, context)


class SymbolNode(object):
    @classmethod
    def parse(cls, string):
        coef = []
        i = 0
        while string[i].isdigit():
            coef.append(string[i])
            i += 1
        coef_val = int(''.join(coef))
        residue_sym = string[i:]
        if string[i] == "(" and string[-1] == ")":
            residue_sym = residue_sym[1:-1]
        return cls(residue_sym, coef_val)

    def __init__(self, symbol, coefficient=1):
        self.symbol = symbol
        self.coefficient = coefficient

    def __hash__(self):
        return hash(self.symbol)

    def __eq__(self, other):
        if isinstance(other, basestring):
            return self.symbol == other
        else:
            return self.symbol == other.symbol and self.coefficient == other.coefficient

    def __ne__(self, other):
        return not self == other

    def __repr__(self):
        return "{}({})".format(self.coefficient, self.symbol)

    def to_monosaccharide(self):
        return MonosaccharideResidue.from_monosaccharide(monosaccharides[self.symbol])


operator_map = {}


def register_operator(cls):
    operator_map[cls.symbol] = cls()
    return cls


@register_operator
class Operator(object):
    symbol = "NoOp"

    def __init__(self, symbol=None):
        if symbol is not None:
            self.symbol = symbol

    def __call__(self, left, right, context):
        return NotImplemented

    def __repr__(self):
        return self.symbol

    operator_map = operator_map


@register_operator
class LessThan(Operator):
    symbol = "<"

    def __call__(self, left, right, context):
        left_val = context[left] * left.coefficient
        right_val = context[right] * right.coefficient
        return left_val < right_val


@register_operator
class LessThanOrEqual(Operator):
    symbol = "<="

    def __call__(self, left, right, context):
        left_val = context[left] * left.coefficient
        right_val = context[right] * right.coefficient
        return left_val <= right_val


@register_operator
class GreaterThan(Operator):
    symbol = ">"

    def __call__(self, left, right, context):
        left_val = context[left] * left.coefficient
        right_val = context[right] * right.coefficient
        return left_val > right_val


@register_operator
class GreaterThanOrEqual(Operator):
    symbol = ">="

    def __call__(self, left, right, context):
        left_val = context[left] * left.coefficient
        right_val = context[right] * right.coefficient
        return left_val >= right_val


@register_operator
class Equal(Operator):
    symbol = "="

    def __call__(self, left, right, context):
        left_val = context[left] * left.coefficient
        right_val = context[right] * right.coefficient
        return left_val == right_val
