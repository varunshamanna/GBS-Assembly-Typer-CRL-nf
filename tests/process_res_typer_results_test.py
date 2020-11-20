import argparse
import unittest
from unittest.mock import patch, call, ANY

from bin.process_res_typer_results import get_arguments, codon2aa, extract_seq_by_id, derive_presence_absence_targets, \
    derive_presence_absence_targets_for_arg_res, six_frame_translate, extract_frame_aa, \
    update_presence_absence_target, update_presence_absence_target_for_arg_res, EOL_SEP, MIN_DEPTH


class TestProcessResTyperResults(unittest.TestCase):

    TEST_LANE = "26189_8#5"
    TEST_GBS_FULLGENES_RESULTS_FILE = "test_data/RES_" + TEST_LANE + "__fullgenes__GBS_Res_Gene-DB_Final__results.txt"
    TEST_ARGANNOT_FULLGENES_RESULTS_FILE = "test_data/ARG_" + TEST_LANE + "__fullgenes__ARG-ANNOT__results.txt"
    TEST_RESFINDER_FULLGENES_RESULTS_FILE = "test_data/RESFI_" + TEST_LANE + "__fullgenes__ResFinder__results.txt"
    TEST_FASTA_FILE = "test_data/test-db.fasta"

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
        self.assertEqual('S', codon2aa('tca'))
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

    def test_extract_frame_aa(self):
        self.assertEqual("IHMVIHL", extract_frame_aa("ATCCACATGGTGATTCATCTA", 1))
        self.assertEqual("IHMVIHL", extract_frame_aa("atccacatggtgattcatcta", 1))
        self.assertEqual("STW*FI", extract_frame_aa("ATCCACATGGTGATTCATCTA", 2))
        self.assertEqual("PHGDSS", extract_frame_aa("ATCCACATGGTGATTCATCTA", 3))
        self.assertEqual("*MNHHVD", extract_frame_aa("ATCCACATGGTGATTCATCTA", 4))
        self.assertEqual("*MNHHVD", extract_frame_aa("atccacatggtgattcatcta", 4))
        self.assertEqual("R*ITMW", extract_frame_aa("ATCCACATGGTGATTCATCTA", 5))
        self.assertEqual("DESPCG", extract_frame_aa("ATCCACATGGTGATTCATCTA", 6))

    def test_six_frame_translate__frame1(self):
        self.assertEqual(
            "VMGKYHPHGDSSIYEAMVRMAQWW",
            six_frame_translate(
                ">5__GYRAGBS__GYRAGBS-1__5" +
                EOL_SEP +
                "GTTATGGGTAAATACCATCCACATGGTGATTCATCTATTTACGAAGCAATGGTGCGTATG" +
                EOL_SEP +
                "GCACAATGGTGG",
                1
            )
        )

        self.assertEqual(
            "HPHGDSSIYDAMVRMSQ",
            six_frame_translate(
                ">7__PARCGBS__PARCGBS-1__7" +
                EOL_SEP +
                "CATCCTCATGGGGATTCCTCTATCTATGACGCGATGGTTCGTATGTCTCAA",
                1
            )
        )

    def test_six_frame_translate__frame2(self):
        self.assertEqual(
            "LWVNTIHMVIHLFTKQWCVWHNG",
            six_frame_translate(
                ">5__GYRAGBS__GYRAGBS-1__5" +
                EOL_SEP +
                "TTTATGGGTAAATACCATCCACATGGTGATTCATCTATTTACGAAGCAATGGTGCGTATG" +
                EOL_SEP +
                "GCACAATGGTGG",
                2
            )
        )

        self.assertEqual(
            "ILMGIPLSMTRWFVCL",
            six_frame_translate(
                ">7__PARCGBS__PARCGBS-1__7" +
                EOL_SEP +
                "AATCCTCATGGGGATTCCTCTATCTATGACGCGATGGTTCGTATGTCTCAA",
                2
            )
        )

    def test_six_frame_translate__frame3(self):
        self.assertEqual(
            "YG*IPSTW*FIYLRSNGAYGTMVEKK",
            six_frame_translate(
                ">5__GYRAGBS__GYRAGBS-1__5" +
                EOL_SEP +
                "TTTATGGGTAAATACCATCCACATGGTGATTCATCTATTTACGAAGCAATGGTGCGTATG" +
                EOL_SEP +
                "GCACAATGGTGGAAAAAAAAAA",
                3
            )
        )

        self.assertEqual(
            "SSWGFLYL*RDGSYVS",
            six_frame_translate(
                ">7__PARCGBS__PARCGBS-1__7" +
                EOL_SEP +
                "AATCCTCATGGGGATTCCTCTATCTATGACGCGATGGTTCGTATGTCTCAA",
                3
            )
        )

    def test_six_frame_translate__frame4(self):
        self.assertEqual(
            "FFFSTIVPYAPLLRK*MNHHVDGIYP*",
            six_frame_translate(
                ">5__GYRAGBS__GYRAGBS-1__5" +
                EOL_SEP +
                "TTTATGGGTAAATACCATCCACATGGTGATTCATCTATTTACGAAGCAATGGTGCGTATG" +
                EOL_SEP +
                "GCACAATGGTGGAAAAAAAAAA",
                4
            )
        )

        self.assertEqual(
            "LRHTNHRVIDRGIPMRI",
            six_frame_translate(
                ">7__PARCGBS__PARCGBS-1__7" +
                EOL_SEP +
                "AATCCTCATGGGGATTCCTCTATCTATGACGCGATGGTTCGTATGTCTCAA",
                4
            )
        )

    def test_six_frame_translate__frame5(self):
        self.assertEqual(
            "FFFPPLCHTHHCFVNR*ITMWMVFTHK",
            six_frame_translate(
                ">5__GYRAGBS__GYRAGBS-1__5" +
                EOL_SEP +
                "TTTATGGGTAAATACCATCCACATGGTGATTCATCTATTTACGAAGCAATGGTGCGTATG" +
                EOL_SEP +
                "GCACAATGGTGGAAAAAAAAAA",
                5
            )
        )

        self.assertEqual(
            "*DIRTIAS*IEESP*G",
            six_frame_translate(
                ">7__PARCGBS__PARCGBS-1__7" +
                EOL_SEP +
                "AATCCTCATGGGGATTCCTCTATCTATGACGCGATGGTTCGTATGTCTCAA",
                5
            )
        )

    def test_six_frame_translate__frame6(self):
        self.assertEqual(
            "FFFHHCAIRTIAS*IDESPCGWYLPI",
            six_frame_translate(
                ">5__GYRAGBS__GYRAGBS-1__5" +
                EOL_SEP +
                "TTTATGGGTAAATACCATCCACATGGTGATTCATCTATTTACGAAGCAATGGTGCGTATG" +
                EOL_SEP +
                "GCACAATGGTGGAAAAAAAAAA",
                6
            )
        )

        self.assertEqual(
            "ETYEPSRHR*RNPHED",
            six_frame_translate(
                ">7__PARCGBS__PARCGBS-1__7" +
                EOL_SEP +
                "AATCCTCATGGGGATTCCTCTATCTATGACGCGATGGTTCGTATGTCTCAA",
                6
            )
        )

    def test_six_frame_translate__frame_range_error(self):
        with self.assertRaises(IndexError):
            six_frame_translate(
                ">5__GYRAGBS__GYRAGBS-1__5" +
                EOL_SEP +
                "TTTATGGGTAAATACCATCCACATGGTGATTCATCTATTTACGAAGCAATGGTGCGTATG" +
                EOL_SEP +
                "GCACAATGGTGGAAAAAAAAAA",
                7
            )

        with self.assertRaises(IndexError):
            six_frame_translate(
                ">5__GYRAGBS__GYRAGBS-1__5" +
                EOL_SEP +
                "TTTATGGGTAAATACCATCCACATGGTGATTCATCTATTTACGAAGCAATGGTGCGTATG" +
                EOL_SEP +
                "GCACAATGGTGGAAAAAAAAAA",
                0
            )

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

    def test_update_presence_absence_target(self):
        depth = MIN_DEPTH+1

        # ============== Test ERM ==================
        drug_res_col_dict = {"EC": "neg"}
        res_target_dict = {"ERM": "neg"}
        update_presence_absence_target("GENE1", "***ERM***", depth, drug_res_col_dict, res_target_dict)
        self.assertEqual({"EC": "GENE1"}, drug_res_col_dict)
        self.assertEqual({"ERM": "pos"}, res_target_dict)
        update_presence_absence_target("GENE2", "***ERM***", depth, drug_res_col_dict, res_target_dict)
        self.assertEqual({"EC": "GENE1:GENE2"}, drug_res_col_dict)
        self.assertEqual({"ERM": "pos"}, res_target_dict)

        # Check low depth
        drug_res_col_dict = {"EC": "neg"}
        res_target_dict = {"ERM": "neg"}
        update_presence_absence_target("GENE2", "***ERM***", MIN_DEPTH-1, drug_res_col_dict, res_target_dict)
        self.assertEqual({"EC": "neg"}, drug_res_col_dict)
        self.assertEqual({"ERM": "neg"}, res_target_dict)

        # ============== Test TET ==================
        drug_res_col_dict = {"TET": "neg"}
        res_target_dict = {"TET": "neg"}
        update_presence_absence_target("GENE1", "***TET***", depth, drug_res_col_dict, res_target_dict)
        self.assertEqual({"TET": "GENE1"}, drug_res_col_dict)
        self.assertEqual({"TET": "pos"}, res_target_dict)
        update_presence_absence_target("GENE2", "***TET***", depth, drug_res_col_dict, res_target_dict)
        self.assertEqual({"TET": "GENE1:GENE2"}, drug_res_col_dict)
        self.assertEqual({"TET": "pos"}, res_target_dict)

        # ============== Test CAT ==================
        drug_res_col_dict = {"OTHER": "neg"}
        res_target_dict = {"CAT": "neg"}
        update_presence_absence_target("GENE1", "***CAT***", depth, drug_res_col_dict, res_target_dict)
        self.assertEqual({"OTHER": "GENE1"}, drug_res_col_dict)
        self.assertEqual({"CAT": "pos"}, res_target_dict)
        update_presence_absence_target("GENE2", "***CAT***", depth, drug_res_col_dict, res_target_dict)
        self.assertEqual({"OTHER": "GENE1:GENE2"}, drug_res_col_dict)
        self.assertEqual({"CAT": "pos"}, res_target_dict)

        # ============== Test LNUB ==================
        drug_res_col_dict = {"EC": "neg"}
        res_target_dict = {"LNUB": "neg"}
        update_presence_absence_target("GENE1", "***LNUB***", depth, drug_res_col_dict, res_target_dict)
        self.assertEqual({"EC": "GENE1"}, drug_res_col_dict)
        self.assertEqual({"LNUB": "pos"}, res_target_dict)
        update_presence_absence_target("GENE2", "***LNUB***", depth, drug_res_col_dict, res_target_dict)
        self.assertEqual({"EC": "GENE1:GENE2"}, drug_res_col_dict)
        self.assertEqual({"LNUB": "pos"}, res_target_dict)

        # ============== Test LSA ==================
        drug_res_col_dict = {"EC": "neg"}
        res_target_dict = {"LSA": "neg"}
        update_presence_absence_target("GENE1", "***LSA***", depth, drug_res_col_dict, res_target_dict)
        self.assertEqual({"EC": "GENE1"}, drug_res_col_dict)
        self.assertEqual({"LSA": "pos"}, res_target_dict)
        update_presence_absence_target("GENE2", "***LSA***", depth, drug_res_col_dict, res_target_dict)
        self.assertEqual({"EC": "GENE1:GENE2"}, drug_res_col_dict)
        self.assertEqual({"LSA": "pos"}, res_target_dict)

        # ============== Test MEF ==================
        drug_res_col_dict = {"EC": "neg"}
        res_target_dict = {"MEF": "neg"}
        update_presence_absence_target("GENE1", "***MEF***", depth, drug_res_col_dict, res_target_dict)
        self.assertEqual({"EC": "GENE1"}, drug_res_col_dict)
        self.assertEqual({"MEF": "pos"}, res_target_dict)
        update_presence_absence_target("GENE2", "***MEF***", depth, drug_res_col_dict, res_target_dict)
        self.assertEqual({"EC": "GENE1:GENE2"}, drug_res_col_dict)
        self.assertEqual({"MEF": "pos"}, res_target_dict)

        # ============== Test misc ==================
        misc_list = ["PARC", "GYRA", "23S1", "23S3", "RPOB1"]
        for allele in misc_list:
            drug_res_col_dict = {}
            res_target_dict = {allele: "neg"}

            update_presence_absence_target("GENE1", "***"+allele+"***", depth, drug_res_col_dict, res_target_dict)
            self.assertEqual({}, drug_res_col_dict)
            self.assertEqual({allele: "pos"}, res_target_dict)

        # ============== Test RPOBN ==================
        drug_res_col_dict = {}
        res_target_dict = {"RPOB2": "neg"}
        update_presence_absence_target("GENE1", "***RPOBN***", depth, drug_res_col_dict, res_target_dict)
        self.assertEqual({}, drug_res_col_dict)
        self.assertEqual({"RPOB2": "pos"}, res_target_dict)

        # ============== Test depth ==================
        drug_res_col_dict = {}
        res_target_dict = {}
        update_presence_absence_target("GENE1", "***RPOBN***", MIN_DEPTH-1, drug_res_col_dict, res_target_dict)
        self.assertEqual({}, drug_res_col_dict)
        self.assertEqual({}, res_target_dict)

        # TODO there is a suspected bug in this perl code - see Python module

    def test_update_presence_absence_target_for_arg_res(self):
        depth = MIN_DEPTH+1

        # ============== Test ERM ==================
        drug_res_col_dict = {"EC": "neg"}
        res_target_dict = {"ERM": "neg"}
        update_presence_absence_target_for_arg_res("GENE1", "***ERM***", depth, drug_res_col_dict, res_target_dict)
        self.assertEqual({"EC": "ERM"}, drug_res_col_dict)
        self.assertEqual({"ERM": "pos"}, res_target_dict)

        drug_res_col_dict = {"EC": "LNU"}
        res_target_dict = {"ERM": "neg"}
        update_presence_absence_target_for_arg_res("GENE1", "***ERM***", depth, drug_res_col_dict, res_target_dict)
        self.assertEqual({"EC": "LNU:ERM"}, drug_res_col_dict)
        self.assertEqual({"ERM": "pos"}, res_target_dict)

        drug_res_col_dict = {"EC": "ERM"}
        res_target_dict = {"ERM": "pos"}
        update_presence_absence_target_for_arg_res("GENE1", "***ERM***", depth, drug_res_col_dict, res_target_dict)
        self.assertEqual({"EC": "ERM"}, drug_res_col_dict)
        self.assertEqual({"ERM": "pos"}, res_target_dict)

        # ============== Test LNU ==================
        drug_res_col_dict = {"EC": "neg"}
        res_target_dict = {"LNUB": "neg"}
        update_presence_absence_target_for_arg_res("GENE1", "***LNU***", depth, drug_res_col_dict, res_target_dict)
        self.assertEqual({"EC": "LNU"}, drug_res_col_dict)
        self.assertEqual({"LNUB": "pos"}, res_target_dict)

        drug_res_col_dict = {"EC": "ERM"}
        res_target_dict = {"LNUB": "neg"}
        update_presence_absence_target_for_arg_res("GENE1", "***LNU***", depth, drug_res_col_dict, res_target_dict)
        self.assertEqual({"EC": "ERM:LNU"}, drug_res_col_dict)
        self.assertEqual({"LNUB": "pos"}, res_target_dict)

        drug_res_col_dict = {"EC": "ERM"}
        res_target_dict = {"LNUB": "pos"}
        update_presence_absence_target_for_arg_res("GENE1", "***LNU***", depth, drug_res_col_dict, res_target_dict)
        self.assertEqual({"EC": "ERM"}, drug_res_col_dict)
        self.assertEqual({"LNUB": "pos"}, res_target_dict)

        # ============== Test LSA ==================
        drug_res_col_dict = {"EC": "neg"}
        res_target_dict = {"LSA": "neg"}
        update_presence_absence_target_for_arg_res("GENE1", "***LSA***", depth, drug_res_col_dict, res_target_dict)
        self.assertEqual({"EC": "LSA"}, drug_res_col_dict)
        self.assertEqual({"LSA": "pos"}, res_target_dict)

        drug_res_col_dict = {"EC": "ERM"}
        res_target_dict = {"LSA": "neg"}
        update_presence_absence_target_for_arg_res("GENE1", "***LSA***", depth, drug_res_col_dict, res_target_dict)
        self.assertEqual({"EC": "ERM:LSA"}, drug_res_col_dict)
        self.assertEqual({"LSA": "pos"}, res_target_dict)

        drug_res_col_dict = {"EC": "ERM"}
        res_target_dict = {"LSA": "pos"}
        update_presence_absence_target_for_arg_res("GENE1", "***LSA***", depth, drug_res_col_dict, res_target_dict)
        self.assertEqual({"EC": "ERM"}, drug_res_col_dict)
        self.assertEqual({"LSA": "pos"}, res_target_dict)

        # ============== Test MEF ==================
        drug_res_col_dict = {"EC": "neg"}
        res_target_dict = {"MEF": "neg"}
        update_presence_absence_target_for_arg_res("GENE1", "***MEF***", depth, drug_res_col_dict, res_target_dict)
        self.assertEqual({"EC": "MEF"}, drug_res_col_dict)
        self.assertEqual({"MEF": "pos"}, res_target_dict)

        drug_res_col_dict = {"EC": "ERM"}
        res_target_dict = {"MEF": "neg"}
        update_presence_absence_target_for_arg_res("GENE1", "***MEF***", depth, drug_res_col_dict, res_target_dict)
        self.assertEqual({"EC": "ERM:MEF"}, drug_res_col_dict)
        self.assertEqual({"MEF": "pos"}, res_target_dict)

        drug_res_col_dict = {"EC": "ERM"}
        res_target_dict = {"MEF": "pos"}
        update_presence_absence_target_for_arg_res("GENE1", "***MEF***", depth, drug_res_col_dict, res_target_dict)
        self.assertEqual({"EC": "ERM"}, drug_res_col_dict)
        self.assertEqual({"MEF": "pos"}, res_target_dict)

        # ============== Test TET ==================
        drug_res_col_dict = {"TET": "neg"}
        res_target_dict = {"TET": "neg"}
        update_presence_absence_target_for_arg_res("GENE1", "***TET***", depth, drug_res_col_dict, res_target_dict)
        self.assertEqual({"TET": "TET"}, drug_res_col_dict)
        self.assertEqual({"TET": "pos"}, res_target_dict)

        drug_res_col_dict = {"TET": "ERM"}
        res_target_dict = {"TET": "neg"}
        update_presence_absence_target_for_arg_res("GENE1", "***TET***", depth, drug_res_col_dict, res_target_dict)
        self.assertEqual({"TET": "ERM:TET"}, drug_res_col_dict)
        self.assertEqual({"TET": "pos"}, res_target_dict)

        drug_res_col_dict = {"TET": "ERM"}
        res_target_dict = {"TET": "pos"}
        update_presence_absence_target_for_arg_res("GENE1", "***TET***", depth, drug_res_col_dict, res_target_dict)
        self.assertEqual({"TET": "ERM"}, drug_res_col_dict)
        self.assertEqual({"TET": "pos"}, res_target_dict)

        # ============== Test CAT ==================
        drug_res_col_dict = {"OTHER": "neg"}
        res_target_dict = {"CAT": "neg"}
        update_presence_absence_target_for_arg_res("GENE1", "***CAT***", depth, drug_res_col_dict, res_target_dict)
        self.assertEqual({"OTHER": "CAT"}, drug_res_col_dict)
        self.assertEqual({"CAT": "pos"}, res_target_dict)

        drug_res_col_dict = {"OTHER": "ERM"}
        res_target_dict = {"CAT": "neg"}
        update_presence_absence_target_for_arg_res("GENE1", "***CAT***", depth, drug_res_col_dict, res_target_dict)
        self.assertEqual({"OTHER": "ERM:CAT"}, drug_res_col_dict)
        self.assertEqual({"CAT": "pos"}, res_target_dict)

        drug_res_col_dict = {"OTHER": "CAT"}
        res_target_dict = {"CAT": "pos"}
        update_presence_absence_target_for_arg_res("GENE1", "***CAT***", depth, drug_res_col_dict, res_target_dict)
        self.assertEqual({"OTHER": "CAT"}, drug_res_col_dict)
        self.assertEqual({"CAT": "pos"}, res_target_dict)

        # ============== Test OTHER ==================
        drug_res_col_dict = {"OTHER": "neg"}
        res_target_dict = {}
        update_presence_absence_target_for_arg_res("GENE1", "***FOO***", depth, drug_res_col_dict, res_target_dict)
        self.assertEqual({"OTHER": "GENE1"}, drug_res_col_dict)
        self.assertEqual({}, res_target_dict)
        update_presence_absence_target_for_arg_res("GENE2", "***FOO***", depth, drug_res_col_dict, res_target_dict)
        self.assertEqual({"OTHER": "GENE1:GENE2"}, drug_res_col_dict)
        self.assertEqual({}, res_target_dict)

        # ============== Test depth ==================
        drug_res_col_dict = {}
        res_target_dict = {}
        update_presence_absence_target_for_arg_res("GENE1", "***CAT***", MIN_DEPTH - 1, drug_res_col_dict, res_target_dict)
        self.assertEqual({}, drug_res_col_dict)
        self.assertEqual({}, res_target_dict)

    @patch('bin.process_res_typer_results.update_presence_absence_target')
    def test_derive_presence_absence_targets(self, mock):

        calls = [call("23S1", "23S1-1", 1135.571, ANY, ANY), call("23S3", "23S3-3", 1265.721, ANY, ANY)]
        derive_presence_absence_targets(self.TEST_GBS_FULLGENES_RESULTS_FILE)
        mock.assert_has_calls(calls, any_order=False)

    @patch('bin.process_res_typer_results.update_presence_absence_target_for_arg_res')
    def derive_presence_absence_targets_for_arg_res(self, mock):
        calls = [
            call("tet(M)", "tet(M)_12", 132.04, ANY, ANY),
            call("tet(M)", "tet(M)_4", 185.331, ANY, ANY),
            call("tet(M)", "tet(M)_10", 120.412, ANY, ANY),
        ]
        derive_presence_absence_targets_for_arg_res(self.TEST_RESFINDER_FULLGENES_RESULTS_FILE)
        mock.assert_has_calls(calls, any_order=False)

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
