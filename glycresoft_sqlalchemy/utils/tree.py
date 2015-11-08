from collections import defaultdict


class FixTree(defaultdict):
    '''
    A hash-based Prefix or Suffix Tree for testing for
    (sub)sequence inclusion. This implementation works for any
    slice-able iterable of hashable objects, not just strings.
    '''
    def __init__(self):
        defaultdict.__init__(self, FixTree)
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
            self.add(sequence[i+1:], label=label)
        return self

    def __contains__(self, iterable):
        layer = self
        satisfied = True
        for i in iterable:
            if i not in layer.keys():
                satisfied = False
                break
            layer = layer[i]
        return satisfied

    __repr__ = dict.__repr__
