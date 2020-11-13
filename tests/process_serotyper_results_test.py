import argparse
import io
import unittest

from bin.process_serotyper_results import write_line, make_gene_dict, get_arguments

class TestProcessResults(unittest.TestCase):

    def test_should_make_gene_dict(self):
        actual = make_gene_dict('../test_data/SERO_26237_7#5__fullgenes__GBS_seroT_Gene-DB_Final__results.txt', 10)
        self.assertEqual(actual, {'III': ['III', 'III-1', '100.0', '267.595', '1snp', '', '0.581', '172', '0.011', '4', '4']}
                         )

    def setUp(self) -> None:
        self.test_stream = io.StringIO()

    def test_should_write_line(self):
        write_line('III', {'III': ['III', 'III-1', '100.0', '267.595', '1snp', '', '0.581', '172', '0.011', '4', '4']}, self.test_stream)
        actual = self.extract_written_content()
        self.assertEqual(actual, """III-1\tIII=imperfect\tIII\t267.595\n""")

    def extract_written_content(self):
        self.test_stream.seek(0)
        return "".join(self.test_stream.readlines())

    def test_arguments(self):
        actual = get_arguments().parse_args(['--input', 'infile', '--output', 'outfile', '--depth-threshold', '10'])
        self.assertEqual(actual,
                         argparse.Namespace(input='infile', output='outfile', depth=10))

    def test_arguments_short_options(self):
        actual = get_arguments().parse_args(['-i', 'infile', '-o', 'outfile', '-d', '10'])
        self.assertEqual(actual,
                         argparse.Namespace(input='infile', output='outfile', depth=10))

    def test_arguments_without_depth_threshold(self):
        actual = get_arguments().parse_args(['-i', 'infile', '-o', 'outfile'])
        self.assertEqual(actual,
                        argparse.Namespace(input='infile', output='outfile', depth=10))
