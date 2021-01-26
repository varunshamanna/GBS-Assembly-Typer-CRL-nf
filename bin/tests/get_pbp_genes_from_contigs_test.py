import os
import pytest
import argparse

from get_pbp_genes_from_contigs import BlastData, SeqData, FragmentData, get_arguments, check_arguments

@pytest.fixture(scope="session")
def blast_data_location():
    yield 'test_data/test_blast_blactam.out'


@pytest.fixture(scope="session")
def seq_data_location():
    yield 'test_data/test_GBS_bLactam_Ref.fasta'


@pytest.fixture(scope="function")
def prep_blast_data(blast_data_location):
    blast_data_to_process = BlastData(blast_data_location)
    yield blast_data_to_process


@pytest.fixture(scope="function")
def prep_seq_data(seq_data_location):
    seq_data_to_analyse = SeqData(seq_data_location)
    yield seq_data_to_analyse


@pytest.fixture(scope="function")
def prep_fragment_data(prep_blast_data, prep_seq_data):
    blast_data_to_process = prep_blast_data
    seq_data_to_analyse = prep_seq_data
    seq_lengths = seq_data_to_analyse.calculate_seq_length()
    best_blast_hits = blast_data_to_process.get_best_hit()
    fragment_data = FragmentData()
    yield fragment_data, best_blast_hits, seq_lengths


@pytest.mark.parametrize("allele,index,stat", [
    ('GBS1A-1', 0, ['.26077_6_118.11', '100.000', '960', '0', '0', '1', '960', '39459', '40418', '0.0', '1773']),
    ('GBS2B-1', 1, ['.26077_6_118.2', '100.000', '16', '0', '0', '534', '549', '101982', '101967', '1.3', '30.7']),
    ('GBS2X-1', 2, ['.26077_6_118.10', '100.000', '15', '0', '0', '611', '625', '8263', '8277', '4.5', '28.8'])
])
def test_read_blast_out(prep_blast_data, allele, index, stat):
    """
    Test output of blast
    """
    blast_data_to_process = prep_blast_data
    actual = blast_data_to_process.get_data()[allele][index]

    assert actual == stat

@pytest.mark.parametrize("allele,stat", [
    ('GBS1A-1', ['.26077_6_118.11', '100.000', '960', '0', '0', '1', '960', '39459', '40418', '0.0', '1773']),
    ('GBS2B-1', ['.26077_6_118.2', '100.000', '1065', '0', '0', '1', '1065', '185772', '186836', '0.0',	'1967']),
    ('GBS2X-1', ['.26077_6_118.11', '100.000', '1038', '0', '0', '1', '1038', '52297', '51260',	'0.0', '1917'])
])
def test_get_best_blast_hit(prep_blast_data, allele, stat):
    """
    Test the best blast hit of each beta lactam
    """
    blast_data_to_process = prep_blast_data

    assert blast_data_to_process.get_best_hit()[allele] == stat


def test_get_start_end_positions(prep_fragment_data):
    """
    Test getting the start and end blactam positions in the contigs
    """
    fragment_data, best_blast_hits, seq_lengths = prep_fragment_data
    fragment_data.get_start_end_positions(best_blast_hits, seq_lengths, 0.5, 0.5)

    assert fragment_data.get_data() == {'GBS1A-1': ('.26077_6_118.11', '39458', '40418', 'forward', '1', '+'),
                                        'GBS2B-1': ('.26077_6_118.2', '185771', '186836', 'forward', '1', '+'),
                                        'GBS2X-1': ('.26077_6_118.11', '51259', '52297', 'reverse', '1', '-')}


@pytest.mark.parametrize("name,content", [
    ('GBS1A-1', ['.26077_6_118.11\t39458\t40418\tforward\t1\t+\n']),
    ('GBS2B-1', ['.26077_6_118.2\t185771\t186836\tforward\t1\t+\n']),
    ('GBS2X-1', ['.26077_6_118.11\t51259\t52297\treverse\t1\t-\n'])
])
def test_write_start_end_positions(prep_fragment_data, name, content):
    """
    Test output file contents containing start and end positions
    """
    fragment_data, best_blast_hits, seq_lengths = prep_fragment_data
    fragment_data.get_start_end_positions(best_blast_hits, seq_lengths, 0.5, 0.5)
    fragment_data_location = 'test_data/' + name + '.bed'
    fragment_data.write_start_end_positions('test_data/')
    fo = open(fragment_data_location, 'r')

    assert fo.readlines() == content


def test_read_seq_data(prep_seq_data):
    """
    Test output of sequence data
    """
    seq_data_to_analyse = prep_seq_data
    actual = seq_data_to_analyse.get_data()

    assert actual['GBS1A-1'] == 'GACATCTACAACAGTGACACTTACATCGCTTATCCAAACAATGAATTACAAATAGCATCTACCATCATGGATGCGACTAATGGTAAAGTCATTGCACAATTAGGCGGGCGTCATCAGAATGAAAATATTTCATTTGGGACAAATCAATCTGTCTTAACAGACCGCGATTGGGGTTCTACAATGAAACCTATCTCAGCTTATGCACCTGCTATTGATAGTGGTGTCTATAATTCAACAGGTCAATCATTAAACGACTCAGTTTACTACTGGCCTGGTACTTCTACTCAACTATATGACTGGGATCGTCAATATATGGGTTGGATGAGTATGCAGACCGCTATTCAACAATCACGTAACGTCCCTGCTGTCAGAGCACTTGAAGCCGCTGGATTAGACGAAGCAAAATCTTTCCTTGAAAAATTAGGCATATACTATCCAGAAATGAACTATTCAAATGCTATTTCAAGTAACAACAGTAGCAGTGATGCAAAATATGGTGCAAGTAGTGAGAAAATGGCAGCGGCTTACTCGGCTTTTGCAAACGGCGGAACTTACTATAAACCGCAATATGTTAATAAAATTGAATTTAGCGATGGAACCAATGATACTTATGCAGCGTCTGGTAGCCGTGCGATGAAAGAGACTACTGCCTACATGATGACGGATATGCTGAAAACAGTACTAACATTTGGTACTGGTACTAAAGCAGCTATCCCTGGTGTTGCACAAGCTGGTAAGACTGGTACTTCCAACTATACGGAAGATGAGTTAGCTAAAATTGAAGCAACTACTGGTATCTACAATAGCGCCGTTGGTACAATGGCTCCTGATGAAAACTTTGTCGGCTATACTTCTAAGTACACAATGGCAATTTGGACTGGTTATAAAAATCGCCTTACACCACTTTATGGTAGCCAACTGGATATTGCTACTGAGGTTTATCGTGCAATGATGTCCTAC'


@pytest.mark.parametrize("allele,length", [
    ('GBS1A-1', 960),
    ('GBS2B-1', 1065),
    ('GBS2X-1', 1038)
])
def test_calculate_seq_length(prep_seq_data, allele, length):
    """
    Test the sequence lengths output
    """
    seq_data_to_analyse = prep_seq_data
    actual = seq_data_to_analyse.calculate_seq_length()

    assert actual[allele] == length


def test_arguments():
    actual = get_arguments().parse_args(
        ['--blast_out_file', 'blast_out_file', '--query_fasta', 'fasta_file',
        '--frac_align_len_threshold', '0.6', '--frac_identity_threshold', '0.6',
        '--output_prefix', 'out_prefix'])
    assert actual == argparse.Namespace(blast_out='blast_out_file',
    fasta_qu='fasta_file', frac_align=0.6, frac_ident=0.6, output='out_prefix')


def test_arguments_short_options():
    actual = get_arguments().parse_args(
        ['-b', 'blast_out_file', '-f', 'fasta_file',
        '-fa', '0.5', '-fi', '0.5', '-o', 'out_prefix'])
    assert actual ==  argparse.Namespace(blast_out='blast_out_file',
    fasta_qu='fasta_file', frac_align=0.5, frac_ident=0.5, output='out_prefix')


@pytest.mark.parametrize("frac_align,frac_ident", [
    (1.1, 0.5),
    (-0.1, 0.5)
])
def test_check_arguments_frac_align(frac_align, frac_ident):
    args = argparse.Namespace(blast_out='blast_out_file',
    fasta_qu='fasta_file', frac_align=frac_align, frac_ident=frac_ident, output='bed_file')

    with pytest.raises(Exception) as exp:
        check_arguments(args)
    assert str(exp.value) == "Invalid frac_align_len_threshold value. Value must be between 0 and 1."


@pytest.mark.parametrize("frac_align,frac_ident", [
    (0.5, 1.1),
    (0.5, -0.1)
])
def test_check_arguments_frac_ident(frac_align, frac_ident):
    args = argparse.Namespace(blast_out='blast_out_file',
    fasta_qu='fasta_file', frac_align=frac_align, frac_ident=frac_ident, output='bed_file')

    with pytest.raises(Exception) as exp:
        check_arguments(args)
    assert str(exp.value) == "Invalid frac_identity_threshold value. Value must be between 0 and 1."
