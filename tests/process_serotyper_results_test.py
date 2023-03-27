import argparse
import io
import unittest
import os

from bin.process_serotyper_results import write_outfile, make_gene_list, get_arguments



class TestProcessResults(unittest.TestCase):

    TEST_SEROTYPE_FULLGENES = 'tests/test_data/input/SERO_26237_7#5__fullgenes__GBS_seroT_Gene-DB_Final__results.txt'
    TEST_SEROTYPE_FULLGENES_MULTI = 'tests/test_data/input/SERO_26237_7#5__fullgenes__GBS_seroT_Gene-DB_Final__results_multi.txt'
    TEST_OUTPUT = 'tests/test_data/output/test_sero_process_results.txt'

    def test_should_make_gene_list(self):
        actual = make_gene_list(self.TEST_SEROTYPE_FULLGENES, 10)
        self.assertEqual(actual, [['GBS-SBG', 'GBS-SBG:III', '100.0', '267.595', '1snp', '', '0.581', '172', '0.011', '4', '4']])

    def test_write_outfile(self):
        gene_dict = make_gene_list(self.TEST_SEROTYPE_FULLGENES, 10)
        write_outfile(gene_dict, self.TEST_OUTPUT)
        f = open(self.TEST_OUTPUT, 'r')
        actual = f.readlines()
        self.assertEqual(actual,
                            ['Matched_Allele\tMatch_Type\tSerotype\tAvgDepth\n','III\tIII=imperfect\tIII\t267.595\n'])
        os.remove(self.TEST_OUTPUT)

    def test_write_outfile_multi(self):
        gene_dict = make_gene_list(self.TEST_SEROTYPE_FULLGENES_MULTI, 10)
        write_outfile(gene_dict, self.TEST_OUTPUT)
        f = open(self.TEST_OUTPUT, 'r')
        actual = f.readlines()
        self.assertEqual(actual,
                            ['Matched_Allele\tMatch_Type\tSerotype\tAvgDepth\n','NT/Ia/III\tNT=imperfect/Ia=identical/III=imperfect\tIa/III\t2.306/137.929/11.098\n'])
        os.remove(self.TEST_OUTPUT)

    def test_arguments(self):
        actual = get_arguments().parse_args(['--srst2_output', 'srst2_output_name', '--sero_db', 'sero_db', '--output', 'outfile', '--min_read_depth', '30.0'])
        self.assertEqual(actual,
                         argparse.Namespace(id='srst2_output_name', db='sero_db', output='outfile', depth=30.0))

    def test_arguments_short_options(self):
        actual = get_arguments().parse_args(['-s', 'srst2_output_name', '-b', 'sero_db', '-o', 'outfile', '-d', '30.0'])
        self.assertEqual(actual,
                         argparse.Namespace(id='srst2_output_name', db='sero_db', output='outfile', depth=30.0))

    def test_arguments_without_depth_threshold(self):
        actual = get_arguments().parse_args(['-s', 'srst2_output_name', '-b', 'sero_db', '-o', 'outfile'])
        self.assertEqual(actual,
                        argparse.Namespace(id='srst2_output_name', db='sero_db', output='outfile', depth=0))
