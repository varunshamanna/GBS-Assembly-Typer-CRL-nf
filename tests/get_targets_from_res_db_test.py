import argparse
import io
import unittest

from bin.get_targets_from_res_db import get_targets, write_line

class TestProcessResults(unittest.TestCase):

    TEST_TARGETS = 'test_data/seqs_of_interest.txt'

    def test_get_targets(self):
        actual = get_targets(self.TEST_TARGETS)
        self.assertEqual(actual, ['7__PARCGBS__PARCGBS-1__7',
            '5__GYRAGBS__GYRAGBS-1__5',
            '11__23S1__23S1-1__11',
            '12__23S3__23S3-3__12',
            '16__RPOBgbs__RPOBgbs-1__16',
            '17__RPOBgbs__RPOBgbs-2__17',
            '18__RPOBgbs__RPOBgbs-3__18',
            '19__RPOBgbs__RPOBgbs-4__19'])

    def setUp(self) -> None:
        self.test_stream = io.StringIO()

    def test_write_line(self):
        line = '>Target\n'
        actual = write_line(line, 'Target', 0, self.test_stream)
        self.assertEqual(actual, 1)

        line = 'ACTG'
        actual = write_line(line, 'Target', 1, self.test_stream)
        self.assertEqual(actual, 1)

        line = '>Seq'
        actual = write_line(line, 'Target', 0, self.test_stream)
        self.assertEqual(actual, 0)
