from functools import wraps
from copy import deepcopy
import sqlitedict
try:
    import cPickle as pickle
except:
    import pickle

def memoize(maxsize=100):
    """Make a memoization decorator. A negative value of `maxsize` means
    no size limit."""
    def deco(f):
        """Memoization decorator. Items of `kwargs` must be hashable."""
        memo = {}

        @wraps(f)
        def func(*args, **kwargs):
            key = (args, frozenset(kwargs.items()))
            if key not in memo:
                if len(memo) == maxsize:
                    memo.popitem()
                memo[key] = f(*args, **kwargs)
            return memo[key]
        func.memo = memo
        return func
    return deco


def memocpy(maxsize=100):
    def deco(f):
        """Memoization decorator. Items of `kwargs` must be hashable."""
        memo = {}

        @wraps(f)
        def func(*args, **kwargs):
            key = (args, frozenset(kwargs.items()))
            if key not in memo:
                if len(memo) == maxsize:
                    memo.popitem()
                memo[key] = f(*args, **kwargs)
            return deepcopy(memo[key])
        return func
    return deco


def memoclone(maxsize=100):
    def deco(f):
        """Memoization decorator. Items of `kwargs` must be hashable."""
        memo = {}

        @wraps(f)
        def func(*args, **kwargs):
            key = (args, frozenset(kwargs.items()))
            if key not in memo:
                if len(memo) == maxsize:
                    memo.popitem()
                memo[key] = f(*args, **kwargs)
            return (memo[key]).clone()
        return func
    return deco


def sqlitedict_memo():
    def deco(f):
        memo = sqlitedict.open()

        @wraps(f)
        def func(*args, **kwargs):
            key = pickle.dumps((args, frozenset(kwargs.items())))
            if key not in memo:
                value = memo[key] = f(*args, **kwargs)
            else:
                value = memo[key]
            return value
        return func
    return deco
