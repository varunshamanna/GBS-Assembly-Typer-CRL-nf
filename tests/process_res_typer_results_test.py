import argparse
import unittest
import os
from unittest.mock import patch, call

from bin.process_res_typer_results import get_arguments, codon2aa, extract_seq_by_id, derive_presence_absence_targets, \
    derive_presence_absence_targets_for_arg_res, six_frame_translate, find_mismatches, drugRes_Col, Res_Targets, \
    get_seq_diffs, update_Bin_Res_arr, update_drugRes_Col, EOL_SEP, geneToRef, geneToTargetSeq, Bin_Res_arr, drugToClass


class TestProcessResTyperResults(unittest.TestCase):

    TEST_LANE = "26189_8#5"
    TEST_GBS_FULLGENES_RESULTS_FILE = "test_data/RES_" + TEST_LANE + "__fullgenes__GBS_Res_Gene-DB_Final__results.txt"
    TEST_ARGANNOT_FULLGENES_RESULTS_FILE = "test_data/ARG_" + TEST_LANE + "__fullgenes__ARG-ANNOT__results.txt"
    TEST_RESFINDER_FULLGENES_RESULTS_FILE = "test_data/RESFI_" + TEST_LANE + "__fullgenes__ResFinder__results.txt"
    TEST_FASTA_FILE = "test_data/test-db.fasta"
    TEST_RES_BAM_FILE = "test_data/RES_" + TEST_LANE + "__26189_8#5.GBS_Res_Gene-DB_Final.sorted.bam"
    TEST_RES_DB = "test_data/GBS_Res_Gene-DB_Final_0.0.1.fasta"

    def test_extract_seq_by_id(self):
        self.assertEqual("465__DfrB2_Tmt__DfrB2__1230" +
            EOL_SEP +
            "ATGGGTCAAAGTAGCGATGAAGCCAACGCTCCCGTTGCAGGGCAGTTTGCGCTTCCCCTG" +
            EOL_SEP,
            extract_seq_by_id("465__DfrB2_Tmt__DfrB2__1230", self.TEST_FASTA_FILE))

        self.assertEqual(
            "465__DfrB2_Tmt__DfrB3__1231" +
            EOL_SEP +
            "ATGGACCAACACAACAATGGAGTCAGTACTCTAGTTGCTGGCCAGTTTGCGCTCCCATCGAAG" +
            EOL_SEP,
            extract_seq_by_id("465__DfrB2_Tmt__DfrB3__1231", self.TEST_FASTA_FILE))

        self.assertEqual(
            "466__DfrB4_Tmt__DfrB4__1236" +
            EOL_SEP +
            "ATGAATGAAGGAAAAAATGAGGTCAGTACTTCAGCTGCTGGCCGGTTCGCATTCCCATCAAACGCCACGTTTGCCTTGGGGGATCGCGTACGCAAGAAGTCTGGCGCTGCTTGGCAGG" +
            EOL_SEP,
            extract_seq_by_id("466__DfrB4_Tmt__DfrB4__1236", self.TEST_FASTA_FILE))

        self.assertIsNone(extract_seq_by_id("FRED", self.TEST_FASTA_FILE))

    def test_codon2aa(self):

        self.assertEqual('S', codon2aa('TCA'))
        self.assertEqual('S', codon2aa('TCC'))
        self.assertEqual('S', codon2aa('TCG'))
        self.assertEqual('S', codon2aa('TCT'))
        self.assertEqual('F', codon2aa('TTC'))
        self.assertEqual('F', codon2aa('TTT'))
        self.assertEqual('L', codon2aa('TTA'))
        self.assertEqual('L', codon2aa('TTG'))
        self.assertEqual('Y', codon2aa('TAC'))
        self.assertEqual('Y', codon2aa('TAT'))
        self.assertEqual('*', codon2aa('TAA'))
        self.assertEqual('*', codon2aa('TAG'))
        self.assertEqual('C', codon2aa('TGC'))
        self.assertEqual('C', codon2aa('TGT'))
        self.assertEqual('*', codon2aa('TGA'))
        self.assertEqual('W', codon2aa('TGG'))
        self.assertEqual('L', codon2aa('CTA'))
        self.assertEqual('L', codon2aa('CTC'))
        self.assertEqual('L', codon2aa('CTG'))
        self.assertEqual('L', codon2aa('CTT'))
        self.assertEqual('P', codon2aa('CCA'))
        self.assertEqual('P', codon2aa('CCC'))
        self.assertEqual('P', codon2aa('CCG'))
        self.assertEqual('P', codon2aa('CCT'))
        self.assertEqual('H', codon2aa('CAC'))
        self.assertEqual('H', codon2aa('CAT'))
        self.assertEqual('Q', codon2aa('CAA'))
        self.assertEqual('Q', codon2aa('CAG'))
        self.assertEqual('R', codon2aa('CGA'))
        self.assertEqual('R', codon2aa('CGC'))
        self.assertEqual('R', codon2aa('CGG'))
        self.assertEqual('R', codon2aa('CGT'))
        self.assertEqual('I', codon2aa('ATA'))
        self.assertEqual('I', codon2aa('ATC'))
        self.assertEqual('I', codon2aa('ATT'))
        self.assertEqual('M', codon2aa('ATG'))
        self.assertEqual('T', codon2aa('ACA'))
        self.assertEqual('T', codon2aa('ACC'))
        self.assertEqual('T', codon2aa('ACG'))
        self.assertEqual('T', codon2aa('ACT'))
        self.assertEqual('N', codon2aa('AAC'))
        self.assertEqual('N', codon2aa('AAT'))
        self.assertEqual('K', codon2aa('AAA'))
        self.assertEqual('K', codon2aa('AAG'))
        self.assertEqual('S', codon2aa('AGC'))
        self.assertEqual('S', codon2aa('AGT'))
        self.assertEqual('R', codon2aa('AGA'))
        self.assertEqual('R', codon2aa('AGG'))
        self.assertEqual('V', codon2aa('GTA'))
        self.assertEqual('V', codon2aa('GTC'))
        self.assertEqual('V', codon2aa('GTG'))
        self.assertEqual('V', codon2aa('GTT'))
        self.assertEqual('A', codon2aa('GCA'))
        self.assertEqual('A', codon2aa('GCC'))
        self.assertEqual('A', codon2aa('GCG'))
        self.assertEqual('A', codon2aa('GCT'))
        self.assertEqual('D', codon2aa('GAC'))
        self.assertEqual('D', codon2aa('GAT'))
        self.assertEqual('E', codon2aa('GAA'))
        self.assertEqual('E', codon2aa('GAG'))
        self.assertEqual('G', codon2aa('GGA'))
        self.assertEqual('G', codon2aa('GGC'))
        self.assertEqual('G', codon2aa('GGG'))
        self.assertEqual('G', codon2aa('GGT'))

        self.assertEqual('A', codon2aa('GC?'))
        self.assertEqual('G', codon2aa('GG?'))
        self.assertEqual('P', codon2aa('CC?'))
        self.assertEqual('T', codon2aa('AC?'))
        self.assertEqual('V', codon2aa('GT?'))
        self.assertEqual('R', codon2aa('CG?'))
        self.assertEqual('S', codon2aa('TC?'))

    def test_codon2aa_unknown_codon(self):

        self.assertEqual('x', codon2aa('+++'))

    def do_six_frame_translate(self, opt):
        """ Helper method """
        actual = six_frame_translate(
            "5__GYRAGBS__GYRAGBS-1__5" +
            EOL_SEP +
            "GTTATGGGTAAATACCATCCACATGGTGATTCATCTATTTACGAAGCAATGGTGCGTATG" +
            EOL_SEP +
            "GCACAATGGTGG",
            opt
        )

        self.assertEqual("", actual)

    def test_six_frame_translate(self):

        for i in range(0, 5):
            self.do_six_frame_translate(i)

    def test_arguments(self):
        actual = get_arguments().parse_args(
            ['--srst2_gbs', 'srst2_gbs', '--srst2_argannot', 'srst2_argannot', '--srst2_resfinder', 'srst2_resfinder',
             '--output', 'output', '--output_bin', 'output_bin'])
        self.assertEqual(actual,
                         argparse.Namespace(srst2_gbs_output='srst2_gbs',
                                            srst2_argannot_output='srst2_argannot',
                                            srst2_resfinder_output='srst2_resfinder',
                                            output='output',
                                            output_bin='output_bin'))

    def test_derive_presence_absence_targets(self):

        derive_presence_absence_targets(self.TEST_GBS_FULLGENES_RESULTS_FILE)
        print(drugRes_Col)
        print(Res_Targets)
        derive_presence_absence_targets_for_arg_res(self.TEST_ARGANNOT_FULLGENES_RESULTS_FILE)
        print(drugRes_Col)
        print(Res_Targets)
        derive_presence_absence_targets_for_arg_res(self.TEST_RESFINDER_FULLGENES_RESULTS_FILE)
        print(drugRes_Col)
        print(Res_Targets)

    def test_arguments_short_options(self):
        actual = get_arguments().parse_args(
            ['-g', 'srst2_gbs', '-a', 'srst2_argannot', '-r', 'srst2_resfinder',
             '-o', 'output', '-b', 'output_bin'])
        self.assertEqual(actual,
                         argparse.Namespace(srst2_gbs_output='srst2_gbs',
                                            srst2_argannot_output='srst2_argannot',
                                            srst2_resfinder_output='srst2_resfinder',
                                            output='output',
                                            output_bin='output_bin'))

    def test_find_amino_acid_mismatches(self):
        actual = find_mismatches([], 'HPHGDSSIYDAMVRMSS', geneToRef['PARC'])
        self.assertEqual(actual, ['Q17S'])

        actual = find_mismatches([], 'HHHGDSSIYDAMVRMSS', geneToRef['PARC'])
        self.assertEqual(actual, ['P2H', 'Q17S'])

        actual = find_mismatches([], 'MMGKYHPHGDSSIYEAMVRMAQWW', geneToRef['GYRA'])
        self.assertEqual(actual, ['V1M'])

        actual = find_mismatches([], 'GGSSQLSQFMDQHNPLSELSHKRRLSALGPGGL', geneToRef['RPOB1'])
        self.assertEqual(actual, ['F1G'])

        actual = find_mismatches([], 'SSQLVRSPGV', geneToRef['RPOB2'])
        self.assertEqual(actual, ['V1S'])

        actual = find_mismatches([], 'TTVAQANSKLNEDGTFAEEIVMGRHQGNNQEFPSSI', geneToRef['RPOB3'])
        self.assertEqual(actual, ['F1T'])

        actual = find_mismatches([], 'IIDPKAPYVGT', geneToRef['RPOB4'])
        self.assertEqual(actual, ['L1I'])

    def test_find_nucleotide_mismatches(self):
        actual = find_mismatches([], 'ATTACCCGCGACAGGACGGAAAGACCCCATGGAG', geneToRef['23S1'])
        self.assertEqual(actual, ['G1A'])

        actual = find_mismatches([], 'ATTACCCGCGACAGGACGGAAAGACCCCATGGAT', geneToRef['23S1'])
        self.assertEqual(actual, ['G1A', 'G34T'])

        actual = find_mismatches([], 'GGGCACGCGAGCTGGGTTCAGAACGTCGTGAGACAGTTCGGTCCCTATCCGTCGCGGGCG', geneToRef['23S3'])
        self.assertEqual(actual, ['C1G'])

    @patch('bin.process_res_typer_results.six_frame_translate')
    def test_get_seq_diffs(self, mock_six_frame_translate):
        mock_six_frame_translate.return_value = 'HPHGDSSIYDAMVRMSQ'
        get_seq_diffs('CATCCTCATGGGGATTCCTCTATCTATGACGCGATGGTTCGTATGTCTCAA', geneToTargetSeq['PARC'], geneToRef['PARC'])
        self.assertEqual(mock_six_frame_translate.call_args_list, [call('CATCCTCATGGGGATTCCTCTATCTATGACGCGATGGTTCGTATGTCTCAA', 1)])

    def test_update_Bin_Res_arr(self):
        actual = update_Bin_Res_arr('PARC', ['Q17S'], Bin_Res_arr)
        self.assertEqual(actual, {
            'PARC':'PARC-Q17S',
            'GYRA': '',
            '23S1': '',
            '23S3': '',
            'RPOB1': '',
            'RPOB2': '',
            'RPOB3': '',
            'RPOB4': ''
        })

        actual = update_Bin_Res_arr('PARC', ['Q18S'], Bin_Res_arr)
        self.assertEqual(actual, {
            'PARC':'PARC-Q18S',
            'GYRA': '',
            '23S1': '',
            '23S3': '',
            'RPOB1': '',
            'RPOB2': '',
            'RPOB3': '',
            'RPOB4': ''
        })

        actual = update_Bin_Res_arr('GYRA', [], Bin_Res_arr)
        self.assertEqual(actual, {
            'PARC':'PARC-Q18S',
            'GYRA': 'GYRA',
            '23S1': '',
            '23S3': '',
            'RPOB1': '',
            'RPOB2': '',
            'RPOB3': '',
            'RPOB4': ''
        })

        actual = update_Bin_Res_arr('23S1', ['G1A', 'G34T'], Bin_Res_arr)
        self.assertEqual(actual, {
            'PARC':'PARC-Q18S',
            'GYRA': 'GYRA',
            '23S1': '23S1-G1A,G34T',
            '23S3': '',
            'RPOB1': '',
            'RPOB2': '',
            'RPOB3': '',
            'RPOB4': ''
        })

    def test_update_drugRes_Col(self):
        actual = update_drugRes_Col('PARC', ['Q17S'], drugRes_Col, drugToClass)
        self.assertEqual(actual, {
            'TET': 'TETM',
            'EC': 'neg',
            'FQ': 'PARC-Q17S',
            'OTHER': 'TetM_Tet:tet(M):tet(M):tet(M)'
        })

        actual = update_drugRes_Col('RPOB1', ['F1G'], drugRes_Col, drugToClass)
        self.assertEqual(actual, {
            'TET': 'TETM',
            'EC': 'neg',
            'FQ': 'PARC-Q17S',
            'OTHER': 'TetM_Tet:tet(M):tet(M):tet(M):RPOB1-F1G'
        })
