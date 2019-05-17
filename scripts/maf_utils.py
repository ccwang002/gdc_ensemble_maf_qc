from collections import namedtuple
from pathlib import Path
import gzip


class MAF:
    """
    General purpose of MAF reader.

    It assumes the file begins wtih some comments lines (starting with `#`),
    then the column header, and the actual variant records in TSV format.

    Arguments:
        pth (pathlib.Path or str): Path object to the MAF file.
    """
    def __init__(self, pth):
        pth = Path(pth)
        self.pth = pth
        if pth.suffix == '.gz':
            self._file = gzip.open(str(pth), 'rt')
        else:
            self._file = open(str(pth))

        # A reader wrapping the underlying file object
        # which also returns the line number
        self._reader = enumerate(self._file, 1)

        # Header comments appear before column
        self.header_comments = []
        self.raw_columns = self.read_header()

        # Set up columns
        self.columns = self.make_columns(self.raw_columns)
        self._record_cls = self.make_record_class()

    def read_header(self):
        """
        Read the header comments and return the parsed the column header.
        """
        line_no, line = next(self._reader)
        while line.startswith('#'):
            self.header_comments.append(line.rstrip('\n'))
            line_no, line = next(self._reader)

        # Treat the first noncomment line as columns
        return line.rstrip('\n').split('\t')

    def make_columns(self, raw_columns):
        """Define the columns a variant record should store."""
        return [c.lower() for c in raw_columns]

    def make_record_class(self):
        """Define the record class."""
        return namedtuple(f'MAFRecord', self.columns)

    def parse_line(self, line, line_no):
        """Given the MAF record as the whole line, construct the record"""
        cols = line.rstrip('\n').split('\t')
        return cols

    def make_record(self, vals):
        """Given the MAF record values from a row, construct the record"""
        return self._record_cls._make(vals)

    def __iter__(self):
        return self

    def __next__(self):
        line_no, line = next(self._reader)
        cols = self.parse_line(line, line_no)
        record = self.make_record(cols)
        return record


class TrailingTabTrimmedMAF(MAF):
    """
    MAF reader that is tolerable to the rows with trailing tabs removed.
    """
    def parse_line(self, line, line_no):
        cols = line.rstrip('\n').split('\t')

        # If the number of columns are fewer, add empty fields to the end.
        expect_n_cols = len(self.raw_columns)
        n_cols = len(cols)
        if n_cols < expect_n_cols:
            empty_trailing_cols = [''] * (expect_n_cols - n_cols)
            cols += empty_trailing_cols

        return cols


