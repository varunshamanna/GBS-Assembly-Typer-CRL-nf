import argparse
import io
import unittest
from unittest.mock import patch, call, ANY

from bin.process_surface_typer_results import get_arguments


class TestProcessSurfaceTyperResults(unittest.TestCase):

    def setUp(self) -> None:
        self.test_stream = io.StringIO()

    @patch('bin.process_surface_typer_results.derive_presence_absence')
    @patch('bin.process_surface_typer_results.create_output_contents')
    @patch('bin.process_surface_typer_results.write_output')
    def test_run(self, mock_write_output, mock_create_output_contents, mock_derive_presence_absence):

        # Given

        args = get_arguments().parse_args(
            ['--srst2_gbs_fullgenes', 'srst2_gbs_fullgenes',
             '--surface_db', 'gbs_surface_db',
             '--min_read_depth', '40', '--output_prefix', 'output'])
        create_output_contents.return_value = 'foobar'

        # When

        run(args)

        # Then

        self.assertEqual(
            mock_derive_presence_absence.call_args_list, [call(args.fullgenes_file_id, args.min_depth)])
        mock_create_output_contents.assert_called_once()
        mock_write_output.assert_has_calls([
            call(ANY, args.output + '_surface_protein_incidence.txt'),
            call(ANY, args.output + '_surface_protein_variants.txt')
        ], any_order=False)

"""
    def test_arguments(self):
        actual = get_arguments().parse_args(
            ['--srst2_gbs_fullgenes', 'srst2_gbs_fullgenes', '--srst2_gbs_consensus', 'srst2_gbs_consensus',
             '--srst2_other_fullgenes', 'srst2_argannot_fullgenes', 'srst2_resfinder_fullgenes',
             '--min_read_depth', '30', '--output_prefix', 'output'])
        self.assertEqual(actual,
                         argparse.Namespace(srst2_gbs_fg_output='srst2_gbs_fullgenes',
                                            srst2_gbs_cs_output='srst2_gbs_consensus',
                                            srst2_other_fg_output=['srst2_argannot_fullgenes',
                                                                   'srst2_resfinder_fullgenes'],
                                            min_depth=30,
                                            output='output'))

    @patch('bin.process_res_typer_results.get_arguments')
    @patch('bin.process_res_typer_results.run')
    def test_main(self, mock_run, mock_get_arguments):
        args = mock_get_arguments.return_value.parse_args()
        main()
        self.assertEqual(ANY, [call(args)])
"""
