import io
import unittest
import argparse

from bin.combine_results import get_content_with_id, get_sero_res_contents, get_arguments

class TestCombineResults(unittest.TestCase):

    TEST_LANE = '25292_2#85'
    TEST_DATA_SEROTYPE = 'test_data/' + TEST_LANE + '_SeroType_Results.txt'
    TEST_DATA_RES_INCIDENCE = 'test_data/' + TEST_LANE + '_res_incidence.txt'
    TEST_DATA_RES_ALLELES = 'test_data/' + TEST_LANE + '_res_alleles.txt'
    TEST_DATA_RES_VARIANTS = 'test_data/' + TEST_LANE + '_res_gbs_variants.txt'

    def test_should_get_content_with_id_alleles(self):
        actual = get_content_with_id(self.TEST_LANE, self.TEST_DATA_RES_ALLELES)
        self.assertEqual(actual, "ID\tEC\tFQ\tOTHER\tTET\n25292_2#85\tERMB(ERMB-1):MEF(MEF-1):23S1:23S3\tPARC:GYRA-V1A,M2Q,G3W,K4W\tMsrD_MLS(MsrD_296):Ant6-Ia_AGly(Ant6-Ia_1633):Sat4A_AGly(Sat4A_586):Aph3-III_AGly(Aph3-III_268):msr(D)(msr(D)_2):ant(6)-Ia(ant(6)-Ia):aph(3')-III(aph(3')-III)\tTETO(TETO-1)\n")

    def test_should_get_content_with_id_variants(self):
        actual = get_content_with_id(self.TEST_LANE, self.TEST_DATA_RES_VARIANTS)
        self.assertEqual(actual, "ID\t23S1\t23S3\tGYRA\tPARC\tRPOBGBS-1\tRPOBGBS-2\tRPOBGBS-3\tRPOBGBS-4\n25292_2#85\t23S1\t23S3\tGYRA-V1A,M2Q,G3W,K4W\tPARC\t\t\t\t\n")

    def test_should_get_sero_res_contents(self):
        actual = get_sero_res_contents(self.TEST_LANE, self.TEST_DATA_SEROTYPE, self.TEST_DATA_RES_INCIDENCE)
        self.assertEqual(actual, 'ID\tSerotype\t23S1\t23S3\tCAT\tERM\tGYRA\tLNU\tLSA\tMEF\tPARC\tRPOBGBS-1\tRPOBGBS-2\tRPOBGBS-3\tRPOBGBS-4\tTET\n25292_2#85\tIII\t+\t+\t-\t+\t+\t-\t-\t+\t+\t-\t-\t-\t-\t+\n')

    def test_arguments(self):
        actual = get_arguments().parse_args(['--id', 'id', '--serotyper_results', 'sero_file', '--res_incidence_results', 'res_file',
                                            '--res_alleles_results', 'allele_file', '--res_variants_results', 'variants_file', '--output', 'output_prefix'])
        self.assertEqual(actual,
                         argparse.Namespace(id='id', sero='sero_file', inc='res_file', alleles='allele_file', variants='variants_file', output='output_prefix'))

    def test_arguments_short_options(self):
        actual = get_arguments().parse_args(['-i', 'id', '-s', 'sero_file', '-r', 'res_file', '-a', 'allele_file', '-v', 'variants_file', '-o', 'output_prefix'])
        self.assertEqual(actual,
                         argparse.Namespace(id='id', sero='sero_file', inc='res_file', alleles='allele_file', variants='variants_file', output='output_prefix'))
