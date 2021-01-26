import os
import pytest
import argparse

from get_pbp_alleles import BlastData, SeqData, get_identical_allele, get_imperfect_allele, write_content, get_arguments

@pytest.fixture(scope="session")
def blast_data_location():
    yield 'test_data/test_blast_PBP_alleles.out'


@pytest.fixture(scope="session")
def blast_imperfect_data_location():
    yield 'test_data/test_blast_imperfect_PBP_alleles.out'


@pytest.fixture(scope="session")
def seq_data_location():
    yield 'test_data/test_GBS1A-1.faa'


@pytest.fixture(scope="function")
def prep_blast_data(blast_data_location):
    blast_data_to_process = BlastData(blast_data_location)
    yield blast_data_to_process


@pytest.fixture(scope="function")
def prep_seq_data(seq_data_location):
    seq_data_to_analyse = SeqData(seq_data_location)
    yield seq_data_to_analyse


@pytest.mark.parametrize("fragment,index,stat", [
    ('.26077_6_118.11:39458-40418(+)', 0, ['1||GBS_1A', '100.000', '320', '0', '0', '1', '320', '1', '320', '0.0', '652']),
    ('.26077_6_118.11:39458-40418(+)', 1, ['138||GBS_1A', '99.688', '320', '1', '0', '1', '320', '1', '320', '0.0', '651'])
])
def test_read_blast_out(prep_blast_data, fragment, index, stat):
    """
    Test output of blast
    """
    blast_data_to_process = prep_blast_data
    actual = blast_data_to_process.get_data()[fragment][index]

    assert actual == stat


def test_get_best_blast_hit(prep_blast_data):
    """
    Test the best blast hit of each beta lactam
    """
    blast_data_to_process = prep_blast_data

    assert blast_data_to_process.get_best_hit()['.26077_6_118.11:39458-40418(+)'] == ['1||GBS_1A', '100.000', '320', '0', '0', '1', '320', '1', '320', '0.0', '652']


def test_read_seq_data(prep_seq_data):
    """
    Test output of sequence data
    """
    seq_data_to_analyse = prep_seq_data
    actual = seq_data_to_analyse.get_data()

    assert actual['.26077_6_118.11:39458-40418(+)'] == 'DIYNSDTYIAYPNNELQIASTIMDATNGKVIAQLGGRHQNENISFGTNQSVLTDRDWGSTMKPISAYAPAIDSGVYNSTGQSLNDSVYYWPGTSTQLYDWDRQYMGWMSMQTAIQQSRNVPAVRALEAAGLDEAKSFLEKLGIYYPEMNYSNAISSNNSSSDAKYGASSEKMAAAYSAFANGGTYYKPQYVNKIEFSDGTNDTYAASGSRAMKETTAYMMTDMLKTVLTFGTGTKAAIPGVAQAGKTGTSNYTEDELAKIEATTGIYNSAVGTMAPDENFVGYTSKYTMAIWTGYKNRLTPLYGSQLDIATEVYRAMMSY'


def test_calculate_seq_length(prep_seq_data):
    """
    Test the sequence lengths output
    """
    seq_data_to_analyse = prep_seq_data
    actual = seq_data_to_analyse.calculate_seq_length()

    assert actual['.26077_6_118.11:39458-40418(+)'] == 320


def test_get_identical_allele(prep_blast_data):
    """
    Test detection of identical alleles
    """
    blast_data_to_process = prep_blast_data
    best_blast_hit = blast_data_to_process.get_best_hit()

    assert get_identical_allele(best_blast_hit) == ['.26077_6_118.11:39458-40418(+)\n1||GBS_1A\n']


def test_get_identical_allele_with_no_identical_hits(blast_imperfect_data_location):
    """
    Test detection of identical alleles
    """
    blast_data_to_process = BlastData(blast_imperfect_data_location)
    best_blast_hit = blast_data_to_process.get_best_hit()

    assert get_identical_allele(best_blast_hit) == []


def test_get_imperfect_allele(blast_imperfect_data_location, prep_seq_data):
    """
    Test detection of imperfect alleles
    """
    blast_data_to_process = BlastData(blast_imperfect_data_location)
    best_blast_hit = blast_data_to_process.get_best_hit()

    seq_data_to_analyse = prep_seq_data

    assert get_imperfect_allele(best_blast_hit, seq_data_to_analyse) == ['>.26077_6_118.11:39458-40418(+)\nDIYNSDTYIAYPNNELQIASTIMDATNGKVIAQLGGRHQNENISFGTNQSVLTDRDWGSTMKPISAYAPAIDSGVYNSTGQSLNDSVYYWPGTSTQLYDWDRQYMGWMSMQTAIQQSRNVPAVRALEAAGLDEAKSFLEKLGIYYPEMNYSNAISSNNSSSDAKYGASSEKMAAAYSAFANGGTYYKPQYVNKIEFSDGTNDTYAASGSRAMKETTAYMMTDMLKTVLTFGTGTKAAIPGVAQAGKTGTSNYTEDELAKIEATTGIYNSAVGTMAPDENFVGYTSKYTMAIWTGYKNRLTPLYGSQLDIATEVYRAMMSY\n']


def test_get_imperfect_allele_with_identical_hits(prep_blast_data, prep_seq_data):
    """
    Test detection of imperfect alleles
    """
    blast_data_to_process = prep_blast_data
    best_blast_hit = blast_data_to_process.get_best_hit()

    seq_data_to_analyse = prep_seq_data

    assert get_imperfect_allele(best_blast_hit, seq_data_to_analyse) == []


def test_write_content_of_identical_allele(prep_blast_data):
    """
    Test contents of written file
    """
    blast_data_to_process = prep_blast_data
    best_blast_hit = blast_data_to_process.get_best_hit()
    identical_allele = get_identical_allele(best_blast_hit)

    output_filename = 'test_data/GBS1A-1_identical_PBP_allele.txt'
    write_content(identical_allele, output_filename)
    fo = open(output_filename, 'r')

    assert fo.readlines() == ['.26077_6_118.11:39458-40418(+)\n', '1||GBS_1A\n']


def test_write_content_of_imperfect_allele(blast_imperfect_data_location, prep_seq_data):
    """
    Test contents of written file
    """
    blast_data_to_process = BlastData(blast_imperfect_data_location)
    best_blast_hit = blast_data_to_process.get_best_hit()
    seq_data_to_analyse = prep_seq_data
    imperfect_allele = get_imperfect_allele(best_blast_hit, seq_data_to_analyse)

    output_filename = 'test_data/GBS1A-1_new_PBP_allele.faa'
    write_content(imperfect_allele, output_filename)
    fo = open(output_filename, 'r')

    assert fo.readlines() == ['>.26077_6_118.11:39458-40418(+)\n', 'DIYNSDTYIAYPNNELQIASTIMDATNGKVIAQLGGRHQNENISFGTNQSVLTDRDWGSTMKPISAYAPAIDSGVYNSTGQSLNDSVYYWPGTSTQLYDWDRQYMGWMSMQTAIQQSRNVPAVRALEAAGLDEAKSFLEKLGIYYPEMNYSNAISSNNSSSDAKYGASSEKMAAAYSAFANGGTYYKPQYVNKIEFSDGTNDTYAASGSRAMKETTAYMMTDMLKTVLTFGTGTKAAIPGVAQAGKTGTSNYTEDELAKIEATTGIYNSAVGTMAPDENFVGYTSKYTMAIWTGYKNRLTPLYGSQLDIATEVYRAMMSY\n']


def test_arguments():
    actual = get_arguments().parse_args(
        ['--blast_out_file', 'blast_out_file', '--query_fasta', 'fasta_file',
        '--output_prefix', 'out_prefix'])
    assert actual == argparse.Namespace(blast_out='blast_out_file',
    fasta_qu='fasta_file', output='out_prefix')


def test_arguments_short_options():
    actual = get_arguments().parse_args(
        ['-b', 'blast_out_file', '-f', 'fasta_file', '-o', 'out_prefix'])
    assert actual ==  argparse.Namespace(blast_out='blast_out_file',
    fasta_qu='fasta_file', output='out_prefix')
