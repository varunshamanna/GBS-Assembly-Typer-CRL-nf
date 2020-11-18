import io
import unittest
import argparse

from bin.combine_results import get_sample_line, get_arguments

class TestCombineResults(unittest.TestCase):

    def test_should_print_line_sample1(self):
        actual = get_sample_line('26077_6#118', '../test_data/26077_6#118_SeroType_Results.txt', '../test_data/26077_6#118_BIN_Res_Results.txt')
        self.assertEqual(actual, '26077_6#118	II	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0')

    def test_should_print_line_sample2(self):
        actual = get_sample_line('26237_7#5', '../test_data/26237_7#5_SeroType_Results.txt', '../test_data/26237_7#5_BIN_Res_Results.txt')
        self.assertEqual(actual, '26237_7#5	III	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0')

    def test_should_print_line_sample3(self):
        actual = get_sample_line('25292_2#85', '../test_data/25292_2#85_SeroType_Results.txt', '../test_data/25292_2#85_BIN_Res_Results.txt')
        self.assertEqual(actual, '25292_2#85	III	0	0	0	1	0	0	0	0	0	0	0	0	0	1	0	0	0	0	1')

    def test_should_print_line_multi_serotype(self):
        actual = get_sample_line('MultiSeroType', '../test_data/MultiSeroType_SeroType_Results.txt', '../test_data/MultiSeroType_BIN_Res_Results.txt')
        self.assertEqual(actual, 'MultiSeroType	III;II	0	0	0	1	0	0	0	0	0	0	0	0	0	1	0	0	0	0	1')

    def test_arguments(self):
        actual = get_arguments().parse_args(['--id', 'id', '--serotyper_results', 'sero_file', '--res_typer_results', 'res_file'])
        self.assertEqual(actual,
                         argparse.Namespace(id='id', sero='sero_file', res='res_file'))

    def test_arguments_short_options(self):
        actual = get_arguments().parse_args(['-i', 'id', '-s', 'sero_file', '-r', 'res_file'])
        self.assertEqual(actual,
                         argparse.Namespace(id='id', sero='sero_file', res='res_file'))
