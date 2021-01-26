import os
import pytest
import argparse

from translate_pbp_genes import SeqData, get_arguments, write_seq_dict


@pytest.fixture(scope="session")
def seq_data_location():
    yield 'test_data/test_blactam_contig_fragments.fasta'


@pytest.fixture(scope="function")
def prep_seq_data(seq_data_location):
    seq_data_to_analyse = SeqData(seq_data_location)
    yield seq_data_to_analyse


def test_read_seq_data(prep_seq_data):
    """
    Test output of sequence data
    """
    seq_data_to_analyse = prep_seq_data
    actual = seq_data_to_analyse.get_data()

    assert actual['GBS1A-1_.26077_6_118.11:39458-40418(+)'] == 'GACATCTACAACAGTGACACTTACATCGCTTATCCAAACAATGAATTACAAATAGCATCTACCATCATGGATGCGACTAATGGTAAAGTCATTGCACAATTAGGCGGGCGTCATCAGAATGAAAATATTTCATTTGGGACAAATCAATCTGTCTTAACAGACCGCGATTGGGGTTCTACAATGAAACCTATCTCAGCTTATGCACCTGCTATTGATAGTGGTGTCTATAATTCAACAGGTCAATCATTAAACGACTCAGTTTACTACTGGCCTGGTACTTCTACTCAACTATATGACTGGGATCGTCAATATATGGGTTGGATGAGTATGCAGACCGCTATTCAACAATCACGTAACGTCCCTGCTGTCAGAGCACTTGAAGCCGCTGGATTAGACGAAGCAAAATCTTTCCTTGAAAAATTAGGCATATACTATCCAGAAATGAACTATTCAAATGCTATTTCAAGTAACAACAGTAGCAGTGATGCAAAATATGGTGCAAGTAGTGAGAAAATGGCAGCGGCTTACTCGGCTTTTGCAAACGGCGGAACTTACTATAAACCGCAATATGTTAATAAAATTGAATTTAGCGATGGAACCAATGATACTTATGCAGCGTCTGGTAGCCGTGCGATGAAAGAGACTACTGCCTACATGATGACGGATATGCTGAAAACAGTACTAACATTTGGTACTGGTACTAAAGCAGCTATCCCTGGTGTTGCACAAGCTGGTAAGACTGGTACTTCCAACTATACGGAAGATGAGTTAGCTAAAATTGAAGCAACTACTGGTATCTACAATAGCGCCGTTGGTACAATGGCTCCTGATGAAAACTTTGTCGGCTATACTTCTAAGTACACAATGGCAATTTGGACTGGTTATAAAAATCGCCTTACACCACTTTATGGTAGCCAACTGGATATTGCTACTGAGGTTTATCGTGCAATGATGTCCTAC'


def test_six_frame_translate(prep_seq_data):
    """
    Test amino acid sequence translation
    """
    seq_data_to_analyse = prep_seq_data
    seq_data_to_analyse.translate_content(1)
    seq_data = seq_data_to_analyse.get_data()

    assert seq_data['GBS1A-1_.26077_6_118.11:39458-40418(+)'] == 'DIYNSDTYIAYPNNELQIASTIMDATNGKVIAQLGGRHQNENISFGTNQSVLTDRDWGSTMKPISAYAPAIDSGVYNSTGQSLNDSVYYWPGTSTQLYDWDRQYMGWMSMQTAIQQSRNVPAVRALEAAGLDEAKSFLEKLGIYYPEMNYSNAISSNNSSSDAKYGASSEKMAAAYSAFANGGTYYKPQYVNKIEFSDGTNDTYAASGSRAMKETTAYMMTDMLKTVLTFGTGTKAAIPGVAQAGKTGTSNYTEDELAKIEATTGIYNSAVGTMAPDENFVGYTSKYTMAIWTGYKNRLTPLYGSQLDIATEVYRAMMSY'


def test_write_seq_file(prep_seq_data):
    """
    Test output of seq file
    """
    seq_data_to_analyse = prep_seq_data
    seq_data_to_analyse.translate_content(1)
    seq_data = seq_data_to_analyse.get_data()
    out_seq_data_location = 'test_data/GBS1A-1.faa'
    write_seq_dict(seq_data, out_seq_data_location)
    fo = open(out_seq_data_location, 'r')

    assert fo.readlines() == ['>GBS1A-1_.26077_6_118.11:39458-40418(+)\n', 'DIYNSDTYIAYPNNELQIASTIMDATNGKVIAQLGGRHQNENISFGTNQSVLTDRDWGSTMKPISAYAPAIDSGVYNSTGQSLNDSVYYWPGTSTQLYDWDRQYMGWMSMQTAIQQSRNVPAVRALEAAGLDEAKSFLEKLGIYYPEMNYSNAISSNNSSSDAKYGASSEKMAAAYSAFANGGTYYKPQYVNKIEFSDGTNDTYAASGSRAMKETTAYMMTDMLKTVLTFGTGTKAAIPGVAQAGKTGTSNYTEDELAKIEATTGIYNSAVGTMAPDENFVGYTSKYTMAIWTGYKNRLTPLYGSQLDIATEVYRAMMSY\n']


def test_arguments():
    actual = get_arguments().parse_args(
        ['--blactam_fasta', 'fasta_file', '--output_file', 'out_file'])
    assert actual == argparse.Namespace(fasta='fasta_file', output='out_file')


def test_arguments_short_options():
        actual = get_arguments().parse_args(
            ['-f', 'fasta_file', '-o', 'out_file'])
        assert actual == argparse.Namespace(fasta='fasta_file', output='out_file')
