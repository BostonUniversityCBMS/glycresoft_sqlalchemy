from contextlib import contextmanager
from collections import defaultdict


def identity(x):
    return x


class TreeReprMixin(object):
    def __repr__(self):
        base = dict(self)
        return repr(base)


class PrefixTree(TreeReprMixin, defaultdict):
    '''
    A hash-based Prefix or Suffix Tree for testing for
    sequence inclusion. This implementation works for any
    slice-able sequence of hashable objects, not just strings.
    '''
    def __init__(self):
        defaultdict.__init__(self, PrefixTree)
        self.labels = set()

    def add(self, sequence, label=None):
        layer = self
        if label is None:
            label = sequence
        if label:
            layer.labels.add(label)
        for i in range(len(sequence)):
            layer = layer[sequence[i]]
            if label:
                layer.labels.add(label)

        return self

    def add_ngram(self, sequence, label=None):
        if label is None:
            label = sequence
        for i in range(1, len(sequence) + 1):
            self.add(sequence[:i], label)

    def __contains__(self, sequence):
        layer = self
        satisfied = True
        for i in sequence:
            if i not in layer.keys():
                satisfied = False
                break
            layer = layer[i]
        return satisfied

    def depth_in(self, sequence):
        layer = self
        count = 0
        for i in sequence:
            if i not in layer.keys():
                break
            count += 1
        return count

    def subsequences_of(self, sequence):
        layer = self
        for i in sequence:
            layer = layer[i]
        return layer.labels

    def __iter__(self):
        return iter(self.labels)


class SuffixTree(PrefixTree):
    '''
    A hash-based Prefix or Suffix Tree for testing for
    sequence inclusion. This implementation works for any
    slice-able sequence of hashable objects, not just strings.
    '''
    def __init__(self):
        defaultdict.__init__(self, SuffixTree)
        self.labels = set()

    def add_ngram(self, sequence, label=None):
        if label is None:
            label = sequence
        for i in range(len(sequence)):
            self.add(sequence[i:], label=label)


class KeyTransformingPrefixTree(PrefixTree):
    def __init__(self, transformer):
        defaultdict.__init__(self, lambda: KeyTransformingPrefixTree(transformer))
        self.transformer = transformer
        self.labels = set()

    def add(self, sequence, label=None):
        layer = self
        t = self.transformer
        if label is None:
            label = sequence
        if label:
            layer.labels.add(label)
        for i in range(len(sequence)):
            layer = layer[t(sequence[i])]
            if label:
                layer.labels.add(label)
            # self.add(sequence[i+1:], label=label)
        return self

    def __contains__(self, sequence):
        layer = self
        satisfied = True
        t = self.transformer
        for i in sequence:
            if t(i) not in layer.keys():
                satisfied = False
                break
            layer = layer[i]
        return satisfied

    def depth_in(self, sequence):
        layer = self
        count = 0
        t = self.transformer
        for i in sequence:
            if t(i) not in layer.keys():
                break
            count += 1
        return count

    def subsequences_of(self, sequence):
        layer = self
        t = self.transformer
        for i in sequence:
            layer = layer[t(i)]
        return layer.labels

    __repr__ = dict.__repr__

    def __iter__(self):
        return iter(self.labels)

    @contextmanager
    def with_identity(self):
        t = self.transformer
        self.transformer = identity
        yield
        self.transformer = t


class KeyTransformingSuffixTree(KeyTransformingPrefixTree):
    def __init__(self, transformer):
        defaultdict.__init__(self, lambda: KeyTransformingSuffixTree(transformer))
        self.transformer = transformer
        self.labels = set()

    def add_ngram(self, sequence, label=None):
        if label is None:
            label = sequence
        for i in range(len(sequence)):
            self.add(sequence[i:], label=label)
