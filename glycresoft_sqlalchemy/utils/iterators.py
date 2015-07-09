from collections import Iterable


class peekable(Iterable):

    def __init__(self, iterable):
        self.iterable = iter(iterable)
        try:
            self._lookahead = self.iterable.next()
        except StopIteration:
            self._lookahead = StopIteration

    def __iter__(self):
        # Counter to determine if we've passed through the
        # inner loop at all
        i = 0
        # Default stopping value for x to make accidental reuse
        # of the iterator safer
        x = StopIteration
        for x in self.iterable:
            val = self._lookahead
            self._lookahead = x
            yield (val, x)
            i += 1
            x = StopIteration
        else:
            val = self._lookahead
            self._lookahead = x
        if i > 0:
            yield val, x
        raise StopIteration()

    @property
    def peek(self):
        return self._lookahead
