import itertools
import textwrap
import argparse
import sys

from glycresoft_sqlalchemy.utils import pager
from glycresoft_sqlalchemy.data_model import DatabaseManager, Hypothesis


INDENT = 4
TITLE_SHIFT = 0


class Section(object):
    def __init__(self, title, contents, level=0):
        self.title = title
        self.contents = contents
        self.level = level

    def __iter__(self):
        indent = " " * (self.level * INDENT)
        if len(self.title) > 0:
            yield indent + self.title + ":"
        for line in self.contents:
            if isinstance(line, Section):
                for sub in line:
                    yield sub
            else:
                yield indent + line

    def __repr__(self):
        return "Section(title=%s, level=%d)" % (self.title, self.level)


class InlineSection(Section):
    def __iter__(self):
        indent = " " * ((self.level) * INDENT)
        title_part = indent + self.title + ":"
        titled = len(self.title) < 1
        for line in self.contents:
            if isinstance(line, Section):
                if titled:
                    for sub in line:
                        yield sub
                else:
                    for sub in line:
                        if not titled:
                            chunk = sub.lstrip()
                            yield "%s %s" % (title_part, chunk)
                            titled = True
                        else:
                            yield sub
            else:
                if not titled:
                    yield "%s %s" % (title_part, line)
                else:
                    yield indent + line
            titled = True


class Summary(object):
    def __init__(self, sections):
        self.sections = sections

    def __iter__(self):
        for line in itertools.chain.from_iterable(self.sections):
            yield line


def _summary_dict(data, level=0):
    for k, v in data.items():
        if isinstance(v, dict):
            yield Section(k, _summary_dict(v, level + 1), level + 1)
        elif isinstance(v, (list, tuple)):
            yield Section(k, _summary_collection(v, level + 1), level + 1)
        else:
            yield _summary_scalar("%s: %s" % (k, v), level + 1)


def _summary_scalar(data, level=0):
    return str(data)


def _summary_collection(data, level=0):
    data = list(data)

    results = []
    all_scalar = True
    for i, item in enumerate(data):
        if isinstance(item, (list, tuple)):
            results.append(Section(str(i), _summary_collection(item, level + 1), level + 1))
            all_scalar = False
        elif isinstance(item, dict):
            results.append(Section(str(i), _summary_dict(item, level + 1), level + 1))
            all_scalar = False
        else:
            results.append(_summary_scalar(item, level + 1))
    if all_scalar:
        yield Section("", textwrap.wrap("[%s]" % ', '.join(results), 80), level + 1)
    else:
        for item in results:
            yield item


def stream_hypotheses(session):
    for hypothesis in session.query(Hypothesis).filter(~Hypothesis.is_decoy):
        yield Section("Hypothesis", _summary_dict(hypothesis.summarize()))


def main(database_path, mode='text'):
    session = DatabaseManager(database_path).session()
    if mode == 'text':
        for line in Summary(stream_hypotheses(session)):
            print(line)
    elif mode == "pager":
        pager.page(iter(Summary(stream_hypotheses(session))))
    else:
        print("Unknown Interface - using text.")
        for line in Summary(stream_hypotheses(session)):
            print(line)

    session.close()


app = argparse.ArgumentParser("summarize")
app.add_argument("database_path")


def detect_mode():
    if sys.stdout.isatty():
        return "pager"
    else:
        return "text"


def taskmain():
    args = app.parse_args()
    main(args.database_path, detect_mode())


if __name__ == '__main__':
    taskmain()
