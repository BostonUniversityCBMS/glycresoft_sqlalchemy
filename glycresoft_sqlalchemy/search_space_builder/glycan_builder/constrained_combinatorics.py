import datetime
import re

from glycresoft_sqlalchemy.data_model import MS1GlycanHypothesis, TheoreticalGlycanComposition, PipelineModule
from glycresoft_sqlalchemy.utils.database_utils import get_or_create
from glypy import GlycanComposition
from glycresoft_sqlalchemy.search_space_builder.glycan_builder import registry
from glycresoft_sqlalchemy.search_space_builder.glycan_builder.composition_source import make_motif_lookup_index

import logging
from itertools import product

logger = logging.getLogger("glycan_composition_constrained_combinatorics")


motifs = {
    "N-Glycan": "N-Glycan core basic 1",
    "O-Glycan": "O-Glycan core 1"
}


def resolve_motif(name, lookup):
    return lookup[motifs[name]]


@registry.composition_source_type.register("constrained_combinatorics")
class ConstrainedCombinatoricsGlycanHypothesisBuilder(PipelineModule):
    """
    Summary
    
    Attributes
    ----------
    constraints_list : TYPE
        Description
    derivatization : TYPE
        Description
    hypothesis_id : TYPE
        Description
    HypothesisType : TYPE
        Description
    manager : TYPE
        Description
    options : TYPE
        Description
    reduction : TYPE
        Description
    rules_file : TYPE
        Description
    rules_table : TYPE
        Description
    """
    HypothesisType = MS1GlycanHypothesis

    def __init__(self, database_path, rules_file=None, hypothesis_id=None,
                 rules_table=None, constraints_list=None, derivatization=None,
                 reduction=None, *args, **kwargs):
        self.manager = self.manager_type(database_path)
        self.rules_file = rules_file
        self.rules_table = rules_table
        self.constraints_list = constraints_list
        self.hypothesis_id = hypothesis_id
        self.derivatization = derivatization
        self.reduction = reduction
        self.options = kwargs

    def run(self):
        if self.rules_table is None:
            if self.rules_file is None:
                raise Exception("Must provide a text file of glycan composition rules")
            self.rules_table, self.constrains_list = parse_rules_from_file(self.rules_file)

        self.manager.initialize()
        session = self.manager.session()

        motif_lookup = make_motif_lookup_index(session)

        hypothesis, _ = get_or_create(session, self.HypothesisType, id=self.hypothesis_id)
        hypothesis.name = self.options.get(
            "hypothesis_name",
            "glycan-hypothesis-%s" % datetime.datetime.strftime(
                datetime.datetime.now(), "%Y%m%d-%H%M%S"))
        hypothesis.parameters = hypothesis.parameters or {}
        hypothesis.parameters['rules_table'] = self.rules_table
        hypothesis.parameters['constraints'] = self.constraints_list
        session.add(hypothesis)
        session.commit()

        hypothesis_id = self.hypothesis_id = hypothesis.id

        generator = CombinatoricCompositionGenerator(
            rules_table=self.rules_table, constraints=self.constraints_list)

        acc = []
        for composition, classifications in generator:
            mass = composition.mass()
            serialized = composition.serialize()
            matched_motifs = []
            if classifications:
                for cls in classifications:
                    matched_motifs.append(resolve_motif(cls, motif_lookup))
            rec = TheoreticalGlycanComposition(
                calculated_mass=mass, composition=serialized,
                hypothesis_id=hypothesis_id)
            rec.motifs = matched_motifs
            acc.append(rec)

            if len(acc) > 1000:
                session.add_all(acc)
                session.commit()
                acc = []

        session.add_all(acc)
        session.commit()
        acc = []

        session.close()
        return hypothesis_id


def descending_combination_counter(counter):
    keys = counter.keys()
    count_ranges = map(lambda lo_hi: range(lo_hi[0], lo_hi[1] + 1), counter.values())
    for combination in product(*count_ranges):
        yield dict(zip(keys, combination))


class CombinatoricCompositionGenerator(object):
    """
    Summary

    Attributes
    ----------
    constraints : list
        Description
    lower_bound : list
        Description
    residue_list : list
        Description
    rules_table : dict
        Description
    upper_bound : list
        Description
    """
    @staticmethod
    def build_rules_table(residue_list, lower_bound, upper_bound):
        rules_table = {}
        for i, residue in enumerate(residue_list):
            lower = lower_bound[i]
            upper = upper_bound[i]
            rules_table[residue] = (lower, upper)
        return rules_table

    def __init__(self, residue_list=None, lower_bound=None, upper_bound=None, constraints=None, rules_table=None):
        self.residue_list = residue_list or []
        self.lower_bound = lower_bound or []
        self.upper_bound = upper_bound or []
        self.constraints = constraints or []
        self.rules_table = rules_table

        if len(self.constraints) > 0 and not isinstance(self.constraints[0], CompositionConstraint):
            self.constraints = list(map(CompositionConstraint.from_list, self.constraints))

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
            combin = Solution(combin)
            for constraint in self.constraints:
                if not constraint(combin):
                    passed = False
                    break
            classifications = [classifier.classification for classifier in
                               [is_n_glycan_classifier, is_o_glycan_classifier]
                               if classifier(combin)]
            if passed:
                yield GlycanComposition(**combin.context), classifications

    __iter__ = generate

    def __repr__(self):
        return repr(self.rules_table) + '\n' + repr(self.constraints)


def tryopen(obj):
    if hasattr(obj, "read"):
        return obj
    else:
        return open(obj)


def parse_rules_from_file(path):
    """
    Summary
    
    Parameters
    ----------
    path : TYPE
        Description
    
    Raises
    ------
    Exception
        Description
    
    Returns
    -------
    name : TYPE
        Description
    """
    ranges = []
    constraints = []
    stream = tryopen(path)

    def cast(parts):
        return parts[0], int(parts[1]), int(parts[2])

    for line in stream:
        if line.startswith(";"):
            continue
        parts = line.replace("\n", "").split(" ")
        if len(parts) == 3:
            ranges.append(cast(parts))
        elif len(parts) == 1:
            if parts[0] in ["\n", "\r", ""]:
                break
            else:
                raise Exception("Could not interpret line '%r'" % parts)

    for line in stream:
        if line.startswith("#"):
            continue
        line = line.replace("\n", "")
        if line in ["\n", "\r", ""]:
            break
        constraints.append(CompositionConstraint.parse(line))

    rules_table = CombinatoricCompositionGenerator.build_rules_table(*zip(*ranges))
    return rules_table, constraints


class CompositionConstraint(object):
    """
    A wrapper around an ExpressionNode that is callable
    to test for satisfaction.

    Attributes
    ----------
    expression : ExpressionNode
        The symbolic expression that a composition must
        satisfy in order to be retained.
    """
    @classmethod
    def parse(cls, string):
        return CompositionConstraint(ExpressionNode.parse(string))

    @classmethod
    def from_list(cls, sym_list):
        return cls.parse(' '.join(sym_list))

    def __init__(self, expression):
        self.expression = expression

    def __repr__(self):
        return "{}".format(self.expression)

    def __call__(self, context):
        """
        Test for satisfaction of :attr:`expression`

        Parameters
        ----------
        context : Solution
            Context to evaluate :attr:`expression` in

        Returns
        -------
        bool
        """
        return context[self.expression]

    def __and__(self, other):
        """
        Construct an :class:`AndCompoundConstraint` from
        `self` and `other`

        Parameters
        ----------
        other : CompositionConstraint

        Returns
        -------
        AndCompoundConstraint
        """
        return AndCompoundConstraint(self, other)

    def __or__(self, other):
        """
        Construct an :class:`OrCompoundConstraint` from
        `self` and `other`

        Parameters
        ----------
        other : CompositionConstraint

        Returns
        -------
        OrCompoundConstraint
        """
        return OrCompoundConstraint(self, other)


class SymbolNode(object):
    """
    A node representing a single symbol or variable with a scalar multiplication
    coefficient. When evaluated in the context of a :class:`Solution`, it's value
    will be substituted for the value hashed there.

    A SymbolNode is equal to its symbol string and hashes the same way.

    Attributes
    ----------
    coefficient : int
        Multiplicative coefficient to scale the value by
    symbol : str
        A name
    """
    @classmethod
    def parse(cls, string):
        coef = []
        i = 0
        while i < len(string) and string[i].isdigit():
            coef.append(string[i])
            i += 1
        coef_val = int(''.join(coef) or 1)
        residue_sym = string[i:]
        if residue_sym == "":
            residue_sym = None
        elif string[i] == "(" and string[-1] == ")":
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
        if self.symbol is not None:
            return "{}({})".format(self.coefficient, self.symbol)
        else:
            return "{}".format(self.coefficient)


class ValueNode(object):
    """
    Represent a single numeric value in an ExpressionNode

    Attributes
    ----------
    value : int
        The value of this node
    """
    coefficient = 1

    @classmethod
    def parse(cls, string):
        try:
            value = int(string.replace(" ", ""))
        except:
            value = float(string.replace(" ", ""))
        return cls(value)

    def __init__(self, value):
        self.value = value

    def __hash__(self):
        return hash(self.value)

    def __eq__(self, other):
        if isinstance(other, basestring):
            return str(self.value) == other
        else:
            return self.value == other.value

    def __ne__(self, other):
        return not self == other

    def __repr__(self):
        return str(self.value)


def typify(node):
    if isinstance(node, basestring):
        if re.match(r"(\d+(.\d+)?)", node):
            return ValueNode.parse(node)
        else:
            return SymbolNode.parse(node)
    else:
        return node


def parse_expression(string):
    i = 0
    size = len(string)
    paren_level = 0
    current_symbol = ''
    expression_stack = []
    resolver_stack = []

    while i < size:
        c = string[i]
        if c == " ":
            if current_symbol != "":
                expression_stack.append(current_symbol)
            current_symbol = ""
        elif c == "(":
            paren_level += 1
            if current_symbol != "":
                expression_stack.append(current_symbol)
            current_symbol = ""
            resolver_stack.append(expression_stack)
            expression_stack = []
        elif c == ")":
            paren_level -= 1
            if current_symbol != "":
                expression_stack.append(current_symbol)
            current_symbol = ""
            term = collapse_expression_sequence(expression_stack)
            expression_stack = resolver_stack.pop()
            expression_stack.append(term)
        else:
            current_symbol += c
        i += 1

    if current_symbol != "":
        expression_stack.append(current_symbol)

    if len(resolver_stack) > 0:
        raise SyntaxError("Unpaired parenthesis")
    return collapse_expression_sequence(expression_stack)


class ExpressionNode(object):
    """
    Summary
    
    Attributes
    ----------
    coefficient : int
        Description
    left : TYPE
        Description
    op : TYPE
        Description
    parse : TYPE
        Description
    right : TYPE
        Description
    """
    coefficient = 1

    parse = staticmethod(parse_expression)

    def __init__(self, left, op, right):
        self.left = typify(left)
        self.op = operator_map[op]
        self.right = typify(right)

    def __hash__(self):
        return hash((self.left, self.op, self.right))

    def __eq__(self, other):
        if isinstance(other, basestring):
            return str(self) == other
        else:
            return self.left == other.left and self.right == other.right and self.op == other.op

    def __ne__(self, other):
        return not self == other

    def __repr__(self):
        return "{} {} {}".format(self.left, self.op, self.right)

    def evaluate(self, context):
        """
        Execute :attr:`op` on :attr:`left` and :attr:`right` in `context`

        Parameters
        ----------
        context : Solution

        Returns
        -------
        bool or int
        """
        return self.op(self.left, self.right, context)


def collapse_expression_sequence(expression_sequence):
    stack = []
    i = 0
    size = len(expression_sequence)
    if size == 1:
        return typify(expression_sequence[0])

    while i < size:
        next_term = expression_sequence[i]
        stack.append(next_term)
        if len(stack) == 3:
            node = ExpressionNode(*stack)
            if hasattr(node.left, 'op'):
                if node.left.op.precedence < node.op.precedence:
                    node = ExpressionNode(
                        left=node.left.left, op=node.left.op, right=ExpressionNode(
                            left=node.left.right, op=node.op, right=node.right))
            stack = [node]
        i += 1
    if len(stack) != 1:
        raise SyntaxError("Incomplete Expression: %s" % ' '.join(expression_sequence))
    return stack[0]


class Solution(object):
    """
    Summary

    Attributes
    ----------
    context : dict
    """
    def __init__(self, context):
        self.context = context

    def __getitem__(self, node):
        """
        Summary

        Parameters
        ----------
        node : TYPE
            Description

        Returns
        -------
        name : TYPE
            Description
        """
        if isinstance(node, SymbolNode):
            if node.symbol is None:
                return 1
            else:
                return self.context[node]
        elif isinstance(node, ValueNode):
            return node.value
        else:
            return node.evaluate(self)


operator_map = {}


def register_operator(cls):
    """
    Summary
    
    Returns
    -------
    name : TYPE
        Description
    """
    operator_map[cls.symbol] = cls()
    return cls


@register_operator
class Operator(object):
    """
    Summary
    
    Attributes
    ----------
    operator_map : TYPE
        Description
    precedence : int
        Description
    symbol : str
        Description
    """
    symbol = "NoOp"
    precedence = 0

    def __init__(self, symbol=None):
        """
        Summary

        Parameters
        ----------
        symbol : TYPE, optional
            Description
        """
        if symbol is not None:
            self.symbol = symbol

    def __call__(self, left, right, context):
        """
        Summary

        Parameters
        ----------
        left : TYPE
            Description
        right : TYPE
            Description
        context : TYPE
            Description

        Returns
        -------
        name : TYPE
            Description
        """
        return NotImplemented

    def __repr__(self):
        return self.symbol

    def __eq__(self, other):
        """
        Summary

        Parameters
        ----------
        other : TYPE
            Description

        Returns
        -------
        name : TYPE
            Description
        """
        try:
            return self.symbol == other.symbol
        except AttributeError:
            return self.symbol == other

    def __hash__(self):
        return hash(self.symbol)

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


@register_operator
class Subtraction(Operator):
    """
    Summary
    
    Attributes
    ----------
    precedence : int
        Description
    symbol : str
        Description
    """
    symbol = "-"
    precedence = 1

    def __call__(self, left, right, context):
        left_val = context[left] * left.coefficient
        right_val = context[right] * right.coefficient
        return left_val - right_val


@register_operator
class Addition(Operator):
    symbol = "+"
    precedence = 1

    def __call__(self, left, right, context):
        left_val = context[left] * left.coefficient
        right_val = context[right] * right.coefficient
        return left_val + right_val


@register_operator
class Or(Operator):
    symbol = 'or'

    def __call__(self, left, right, context):
        return context[left] or context[right]


@register_operator
class And(Operator):
    symbol = "and"

    def __call__(self, left, right, context):
        return context[left] and context[right]


class OrCompoundConstraint(CompositionConstraint):
    """
    A combination of CompositionConstraint instances which
    is satisfied if either of its components are satisfied.

    Attributes
    ----------
    left : CompositionConstraint
    right : CompositionConstraint
    """
    def __init__(self, left, right):
        self.left = left
        self.right = right

    def __call__(self, context):
        return self.left(context) or self.right(context)

    def __repr__(self):
        return "(({}) or ({}))".format(self.left, self.right)


class AndCompoundConstraint(CompositionConstraint):
    """
    A combination of CompositionConstraint instances which
    is only satisfied if both of its components are satisfied.
    
    Attributes
    ----------
    left : CompositionConstraint
    right : CompositionConstraint
    """
    def __init__(self, left=None, right=None):
        self.left = left
        self.right = right

    def __call__(self, context):
        return self.left(context) and self.right(context)

    def __repr__(self):
        return "(({}) and ({}))".format(self.left, self.right)


class ClassificationConstraint(object):
    """
    Express a classification (label) of a solution by satisfaction
    of a boolean Expression

    Attributes
    ----------
    classification : str
        The label to grant on satisfaction of `constraints`
    constraint : CompositionConstraint
        The constraint that must be satisfied.
    """
    def __init__(self, classification, constraint):
        self.classification = classification
        self.constraint = constraint

    def __call__(self, context):
        """
        Test for satisfaction.

        Parameters
        ----------
        context : Solution

        Returns
        -------
        bool
            Whether the constraint is satisfied
        """
        return self.constraint(context)

    def __repr__(self):
        return "{}\n{}".format(self.classification, self.constraint)


is_n_glycan_classifier = ClassificationConstraint(
    "N-Glycan", (
        CompositionConstraint.from_list(
            ["HexNAc", ">=", "2"]) & CompositionConstraint.from_list(
            ["Hex", ">=", "3"]))
)

is_o_glycan_classifier = ClassificationConstraint(
    "O-Glycan", ((CompositionConstraint.from_list(["HexNAc", ">=", "1"]) &
                  CompositionConstraint.from_list(["Hex", ">=", "1"])) |
                 (CompositionConstraint.from_list(["HexNAc", ">=", "2"])))
)
