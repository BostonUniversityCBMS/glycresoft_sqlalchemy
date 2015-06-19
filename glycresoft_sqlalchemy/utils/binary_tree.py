

class TreeNode(object):
    __slots__ = ["value", "key"]

    def __init__(self, value, key):
        self.value = value
        self.key = key


def _maxheapify(nodes, i):
    left = (2 * i)
    right = (2 * i) + 1
    largest = i

    if nodes[left].key > nodes[largest].key:
        largest = left
    if nodes[right].key > nodes[largest].key:
        largest = right
    if largest != i:
        nodes[i], nodes[largest] = nodes[largest], nodes[i]
        _maxheapify(nodes, largest)


class Heap(object):

    def __init__(self, nodes):
        self.nodes = nodes

    def maxheapify(self):
        nodes = list(self.nodes)
        _maxheapify(nodes, 1)
